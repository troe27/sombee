import os
import glob
import pandas as pd
from cyvcf2 import VCF


def load_dupcaller_stats(folder):
    """
    Load all {sample}.dupcaller_stats.txt files in a folder into a DataFrame.
    
    Each line in the file is assumed to be:
        <word1> <word2> ... <wordN-1> <value>
    where the first N-1 tokens form the variable name and the last token is the value.
    
    Returns:
        DataFrame with variable names as rows (index) and samples as columns.
    """
    records = {}  # sample -> {variable_name: value}

    for path in glob.glob(os.path.join(folder, "*.dupcaller_stats.txt")):
        sample = os.path.basename(path).replace(".dupcaller_stats.txt", "")
        metrics = {}

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # skip empty lines
                
                parts = line.split()
                if len(parts) < 2:
                    continue  # not enough info to parse
                
                # Last token is the value, the rest form the variable name
                *name_tokens, value_token = parts
                var_name = " ".join(name_tokens)

                # Try to parse numeric value
                value = value_token
                try:
                    value = int(value_token)
                except ValueError:
                    try:
                        value = float(value_token)
                    except ValueError:
                        pass  # leave as string if not numeric

                metrics[var_name] = value

        records[sample] = metrics

    # Build DataFrame: rows = variable names, columns = samples
    df = pd.DataFrame.from_dict(records, orient="columns")
    return df



def load_duplex_group_tables(folder, pattern="*.dupcaller_duplex_group_stats.txt"):
    """
    Load all per-sample duplex group tables in `folder` into a single DataFrame.

    Assumes each file has columns:
        duplex_group_strand_composition  duplex_group_number  effective_coverage  mutation_count

    The sample name is taken from the filename by stripping the suffix
    '.duplex_group_stats.txt' (or whatever you use).

    Parameters
    ----------
    folder : str
        Path to the folder with the tables.
    pattern : str, optional
        Glob pattern to match the files. Default: '*.duplex_group_stats.txt'

    Returns
    -------
    df : pandas.DataFrame
        Long-format table with an extra 'sample' column.
    """
    dfs = []

    for path in glob.glob(os.path.join(folder, pattern)):
        fname = os.path.basename(path)
        # adjust this line if your suffix is different
        sample = fname.replace(".dupcaller_duplex_group_stats.txt", "")

        # robust to tabs / spaces
        df = pd.read_csv(path, sep=r"\s+", header=0)
        df["sample"] = sample

        dfs.append(df)

    if not dfs:
        return pd.DataFrame(columns=[
            "duplex_group_strand_composition",
            "duplex_group_number",
            "effective_coverage",
            "mutation_count",
            "sample",
        ])

    return pd.concat(dfs, ignore_index=True)


def load_trinuc_by_duplex_group_tables(folder):
    """
    Load all {sample}.dupcaller_trinuc_by_duplex_group.txt files in a folder.

    Returns:
        dict where keys = sample names
              values = pandas DataFrames for each file
    """
    tables = {}

    pattern = os.path.join(folder, "*.dupcaller_trinuc_by_duplex_group.txt")

    for path in glob.glob(pattern):
        fname = os.path.basename(path)

        # extract sample name
        sample = fname.replace(".dupcaller_trinuc_by_duplex_group.txt", "")

        # read file (tab or whitespace separated)
        df = pd.read_csv(path, sep=r"\s+", engine="python")

        # store
        tables[sample] = df

    return tables



import os
import glob
from cyvcf2 import VCF
import pandas as pd


def load_all_dupcaller_vcfs(folder):
    sample_vcfs = {}

    # SNVs
    for path in glob.glob(os.path.join(folder, "*dupcaller_snv.vcf*")):
        fname = os.path.basename(path)
        sample = fname.split(".dupcaller_snv.vcf")[0]
        sample_vcfs.setdefault(sample, []).append(path)

    # INDELs
    for path in glob.glob(os.path.join(folder, "*dupcaller_indel.vcf*")):
        fname = os.path.basename(path)
        sample = fname.split(".dupcaller_indel.vcf")[0]
        sample_vcfs.setdefault(sample, []).append(path)

    return sample_vcfs



def _extract_dp(var):
    """Try to get DP from INFO, then from FORMAT."""
    dp = var.INFO.get("DP")
    if dp is not None:
        return dp
    # Try FORMAT (per-sample DP)
    try:
        dp_fmt = var.format("DP")
        if dp_fmt is not None and len(dp_fmt) > 0:
            # usually shape (n_samples, 1)
            return dp_fmt[0][0] if hasattr(dp_fmt[0], "__len__") else dp_fmt[0]
    except KeyError:
        pass
    return None


def _extract_af(var):
    """Try to get AF from INFO, then from FORMAT."""
    af = var.INFO.get("AF")
    if af is not None:
        return af
    # Try FORMAT (per-sample AF)
    for key in ("AF", "AF1"):
        try:
            af_fmt = var.format(key)
            if af_fmt is not None and len(af_fmt) > 0:
                return af_fmt[0][0] if hasattr(af_fmt[0], "__len__") else af_fmt[0]
        except KeyError:
            continue
    return None


def _info_to_string(var):
    """
    Flatten cyvcf2 INFO into 'key=value;key2=value2' string.
    Works with cyvcf2.INFO, which is not a normal dict.
    """
    info = var.INFO
    try:
        keys = list(info.keys())
    except Exception:
        return ""

    parts = []
    for k in keys:
        v = info.get(k)
        # cyvcf2 can return arrays
        if isinstance(v, (list, tuple)):
            v = ",".join(map(str, v))
        parts.append(f"{k}={v}")
    return ";".join(parts)


def vcfs_to_long(sample_vcfs):
    """
    Convert per-sample VCFs to a long-form table with columns:

      sample, CHROM_POS, CHROM, POS, REF, ALT,
      variant_type, QUAL, FILTER, DP, AF, INFO
    """
    rows = []

    for sample, paths in sample_vcfs.items():
        for vcf_path in paths:
            if vcf_path.endswith(".dupcaller_snv.vcf"):
                variant_type = "snv"
            elif vcf_path.endswith(".dupcaller_indel.vcf"):
                variant_type = "indel"
            else:
                variant_type = "unknown"

            vcf = VCF(vcf_path)

            for var in vcf:
                chrom = var.CHROM
                pos = var.POS
                ref = var.REF
                alts = var.ALT or [None]

                qual = var.QUAL
                flt = var.FILTER or "PASS"
                dp = _extract_dp(var)
                af = _extract_af(var)
                info_str = _info_to_string(var)   # <-- changed

                for alt in alts:  # one row per ALT
                    rows.append({
                        "sample": sample,
                        "CHROM_POS": f"{chrom}:{pos}",
                        "CHROM": chrom,
                        "POS": pos,
                        "REF": ref,
                        "ALT": alt,
                        "variant_type": variant_type,
                        "QUAL": qual,
                        "FILTER": flt,
                        "DP": dp,
                        "AF": af,
                        "INFO": info_str,
                    })

    return pd.DataFrame(rows)







