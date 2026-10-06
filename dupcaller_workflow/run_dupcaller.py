import subprocess
import pandas as pd
import os
from Bio import SeqIO
import time

# paths
# software
dupcallerpath = "/home/tilman/osmia_store/dupcaller-test/DupCaller/build/scripts-3.10/DupCaller.py"

# input data
ref = "/home/tilman/nanoseq_test/data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
raw_data_path = '/home/tilman/osmia_store/20250804_nanoseq/files/P35208/'
germline_vcf = '/media/osmia/taliadoros/raw_data/aligned/sorted/RGed2/merged/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz'
sampleID_path = '/home/tilman/osmia_store/20250804_nanoseq/Beetle_IDs.csv'


# output data
detag_path = "./test-data/detag"
bam_folder = "./test-data/bams"
output_folder = "./test-data/output"

panel= "/home/tilman/nanoseq_test/data/ref/ref.bed"

# NNNXXXX  Barcode+SKIP
UMI_pattern = "NNNXXXX"

#pigz parameters:
PIGZ_N = 1
PIGZ_THREADS = "16"

# BWA mem parameters
bwa_threads = 12        # m: threads per sample (for both bwa and samtools)
bwa_max_jobs = 6       # n: how many samples to run in parallel


# gatk MarkDuplicate parameters
gatk_max_jobs = 4        # N: how many MarkDuplicates jobs to run in parallel
gatk_java_mem_gb = 24    # heap per job, adjust as you like

# samtools index parameters
samtools_threads_per_index = 8   # samtools -@ threads
samtools_max_jobs = 3            # concurrent indexing processes

# dupcaller call parameters
dupcaller_threads_per_job = 16   # adjust between 12–20
dupcaller_max_jobs = 6  

def collect_fastq_pairs_from_tree(
    root_dir: str,
    pattern: Optional[str] = None,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Recursively walk root_dir, find FASTQ files, and group them into R1/R2
    per sample based ONLY on the file name.

    Parameters
    ----------
    root_dir : str
        Top-level directory to search under.
    pattern : str, optional
        Regex pattern for file names, with named groups:
          - 'sample': sample ID
          - 'read':   read indicator, typically 'R1' or 'R2'

        Default pattern matches common Illumina-style names, e.g.:
          SAMPLE_XYZ_R1_001.fastq.gz
          SAMPLE_XYZ_R2.fastq.gz
          SAMPLE-RNAseq.R1.fastq.gz

    Returns
    -------
    Dict[str, Dict[str, List[str]]]
        {
          "SAMPLE1": {
              "R1": ["/path/to/SAMPLE1_R1_001.fastq.gz", ...],
              "R2": ["/path/to/SAMPLE1_R2_001.fastq.gz", ...],
          },
          "SAMPLE2": {
              "R1": [...],
              "R2": [...],
          },
          ...
        }
    """
    # Default: very forgiving pattern for filenames
    # - sample name: anything up to _R1 / _R2 / .R1 / .R2
    # - read: R1 or R2
    # Examples matched:
    #   MySample_R1_001.fastq.gz
    #   MySample_R2.fastq.gz
    #   MySample.R1.fastq.gz
    if pattern is None:
        pattern = r'^(?P<sample>.+?)[._-](?P<read>R[12])(?:[_.\-]\d+)?\.f(ast)?q(\.gz)?$'

    regex = re.compile(pattern)
    root = Path(root_dir)

    samples: Dict[str, Dict[str, List[str]]] = {}

    # Walk the entire tree and look at all files
    for path in root.rglob("*"):
        if not path.is_file():
            continue

        # Quick filter: only look at things that *could* be FASTQ
        if not any(str(path).endswith(ext) for ext in (".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            continue

        m = regex.match(path.name)
        if not m:
            continue  # filename doesn't match our pattern

        sample = m.group("sample")
        read = m.group("read")   # 'R1' or 'R2'

        samples.setdefault(sample, {}).setdefault(read, []).append(str(path.resolve()))

    return samples


# produce nonexistent output folders
subprocess.call(["mkdir", "-p", detag_path])
subprocess.call(["mkdir", "-p", bam_folder])
subprocess.call(["mkdir","-p", output_folder])

samples = collect_fastq_pairs_from_tree(root_dir='/home/tilman/osmia_store/20250804_nanoseq/files/P35208/')

# load sampleIDs from file
sampleIDs = pd.read_csv(sampleIDs)
mn_dict = {}
for i in sampleIDs["Sample"].drop_duplicates():
    s = list(sampleIDs.loc[sampleIDs["Sample"]==i].loc[sampleIDs.Dilution=="diluted"]['NGI Sample ID'])[0]
    caste = list(sampleIDs.loc[sampleIDs["Sample"]==i].loc[sampleIDs.Dilution=="diluted"]['Caste'])[0]
    mn = list(sampleIDs.loc[sampleIDs["Sample"]==i].loc[sampleIDs.Dilution=="undiluted"]['NGI Sample ID'])[0]
    mn_dict[i] = {'sample':s, "matched_normal":mn, 'caste':caste}

# index the reference
cmd = ["python", dupcallerpath, 'index', '-f', ref]
P = subprocess.Popen(cmd)

# run Dupcaller trim
for key, sample in samples.items():
    read1 = sample['R1'][0] # there is only one file, since its not split over multiple lanes
    read2 = sample['R2'][0]
    cmd = ["python", dupcallerpath, "trim", "-i", read1, "-i2", read2, "-p", UMI_pattern, "-o", f"{detag_path}/{key}"]
    subprocess.Popen(cmd)

# use pigz to gzip the dupcaller output
files = os.listdir(detag_path)
for idx in range(0, len(files), PIGZ_N):
    batch = files[idx:idx + PIGZ_N]
    procs = []

    for f in batch:
        full = os.path.join(detag_path, f)
        print("processing:", full)
        cmd = ["pigz", "-p", PIGZ_THREADS, full]
        print("\t" + " ".join(cmd))
        procs.append(subprocess.Popen(cmd))

    # wait for this batch to finish
    for p in procs:
        p.wait()

# run bwa mem to map the samples to the reference
procs = []

for key, sample in samples.items():
    # build paths
    r1 = f"{detag_path}/{key}_1.fastq.gz"
    r2 = f"{detag_path}/{key}_2.fastq.gz"
    samples[key]["R1_detag"] = r1
    samples[key]["R2_detag"] = r2
    bam = f"{bam_folder}/{key}.bam"

    # one shell pipeline per sample:
    # bwa mem -> samtools sort -> samtools index
    cmd = (
        f"bwa mem -C -t {bwa_threads} "
        f"-R '@RG\\tID:{key}\\tSM:{key}\\tPL:ILLUMINA' "
        f"{ref} {r1} {r2} "
        f"| samtools sort -@ {bwa_threads} -o {bam} - "
        f"&& samtools index -@ {bwa_threads} {bam}"
    )

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    procs.append(p)

    # throttle to at most max_jobs running at once
    while len(procs) >= bwa_max_jobs:
        # drop finished processes from the list
        procs = [pr for pr in procs if pr.poll() is None]
        if len(procs) >= bwa_max_jobs:
            time.sleep(5)

# wait for any remaining jobs
for p in procs:
    p.wait()

# run gatk MarkDuplicates
procs = []

for key, sample in smp.items():
    bam_in = f"{bam_folder}/{key}.bam"
    bam_out = f"{bam_folder}/{key}.mkdped.bam"
    metrics = f"{bam_folder}/{key}.mkdp_metrics.txt"

    # optional: store paths back into smp
    smp[key]["bam_plain"] = bam_in
    smp[key]["bam_markdup"] = bam_out
    smp[key]["mkdp_metrics"] = metrics

    cmd = [
        "gatk",
        "--java-options", f"-Xmx{gatk_java_mem_gb}g",
        "MarkDuplicates",
        "-I", bam_in,
        "-O", bam_out,
        "-M", metrics,
        "--READ_NAME_REGEX",
        "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$",
        "--DUPLEX_UMI",
        "--TAGGING_POLICY", "OpticalOnly",
        "--BARCODE_TAG", "DB",
    ]

    print(" ".join(cmd))
    p = subprocess.Popen(cmd)
    procs.append(p)

    # throttle: keep at most gatk_max_jobs running
    while len(procs) >= gatk_max_jobs:
        procs = [pr for pr in procs if pr.poll() is None]
        if len(procs) >= gatk_max_jobs:
            time.sleep(5)

# wait for remaining jobs
for p in procs:
    p.wait()

#index the duplicate-marked bams

procs = []

for f in os.listdir(bam_folder):
    if not f.endswith(".mkdped.bam"):
        continue

    bam = os.path.join(bam_folder, f)
    cmd = ["samtools", "index", "-@", str(samtools_threads_per_index), bam]
    print(" ".join(cmd))

    p = subprocess.Popen(cmd)
    procs.append(p)

    # throttle
    while len(procs) >= samtools_max_jobs:
        procs = [pr for pr in procs if pr.poll() is None]
        if len(procs) >= samtools_max_jobs:
            time.sleep(5)

# wait for remaining jobs
for p in procs:
    p.wait()

# generate panel bed that includes all chromosomes
regions_list = []
with open(panel, "wt") as handle:
    for record in SeqIO.parse(ref, 'fasta'):
        id = record.id.split(" ")[0]
        stop = len(record.seq)
        regions_list.append(id)
        #print("\t".join([id, "0", str(stop)]))
        handle.write("\t".join([id, "0", str(stop)])+'\n')

# use Dupcaller call

procs = []

for sample, item in mn_dict.items():
    sample_bam = f"{bam_folder}/{item['sample']}.mkdped.bam"
    normal_bam = f"{bam_folder}/{item['matched_normal']}.mkdped.bam"
    out_prefix = f"{output_folder}/{sample}.dupcaller"

    cmd = [ "python", dupcallerpath, "call","-b", sample_bam,"-f", ref,"-r"]+regions_list+["-o", out_prefix, "-p", str(dupcaller_threads_per_job),"-n", normal_bam,"-g", germline_vcf]
    #"-m", noise_mask_bed,
    #]
    print(" ".join(cmd))
    p = subprocess.Popen(cmd,)
    procs.append(p)

    while len(procs) >= dupcaller_max_jobs:
        procs = [pr for pr in procs if pr.poll() is None]
        if len(procs) >= dupcaller_max_jobs:
            time.sleep(5)

for p in procs:
    p.wait()
