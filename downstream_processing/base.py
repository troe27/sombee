import pandas as pd
import os


# load csv per sample
def load_csv(sample_folder,file_loc="tmpNanoSeq/post/results.muts.tsv", simple=False):
    """
    Loads a mutation TSV file for a given sample.

    Parameters:
    ----------
    sample_folder : str
        Path to the sample folder containing the mutation file.
    file_loc : str, optional
        Relative path to the mutation file (default: "tmpNanoSeq/post/results.muts.tsv").
    simple : bool, optional
        If True, returns only the list of mutation pyrkeys; otherwise, returns the full DataFrame (default: False).

    Returns:
    -------
    pandas.DataFrame or list
        Returns a DataFrame of mutations if simple=False; otherwise, returns a list of pyrkeys.
    
    Raises:
    ------
    ValueError:
        If the `simple` parameter is not True or False.
    """
    sample_file = os.path.join(sample_folder, file_loc)
    muts = pd.read_csv(sample_file, sep="\t")
    if simple==False:
        return muts
    elif simple==True:
        return list(muts.pyrkey)
    else:
        raise ValueError("simple must be True or False")


def load_dict(folder="./test_conservative_noTGCA", prefix="test_out", simple=False):
    """
    Loads mutation data from multiple samples within a given folder.

    Parameters:
    ----------
    folder : str, optional
        Path to the folder containing mutation files (default: "./test_conservative_noTGCA").
    prefix : str, optional
        Prefix for identifying relevant sample files (default: "test_out").
    simple : bool, optional
        If True, loads only pyrkeys from each sample file; otherwise, loads full mutation data (default: True).

    Returns:
    -------
    dict
        Dictionary mapping sample file names to their corresponding mutation data (either a DataFrame or a list of IDs).
    """
    somvar_dict = {}
    for i in [j for j in  os.listdir(folder) if j.startswith(prefix)]:
        l = os.path.join(folder, i)
        svs = load_csv(sample_folder=l, simple=simple)
        svs['POS'] = svs['chromStart']+1
        svs["ID"] = [str(k.chrom)+':'+str(k['POS']) for i,k in svs.iterrows()]
        
        somvar_dict[i]=svs
    return somvar_dict

# filter by somatic variants from other samples
def get_sv_overlap(somvar_dict, return_overlap_dict=False):
    overlap_dict = {}
    for key, item in somvar_dict.items():
        for k2, i2 in somvar_dict.items():
            if not key == k2:
                
                s = set.intersection(set(item['ID']), set(i2['ID']))
                if not len(list(s))==0:
                   overlap_dict[k2+":"+key] = list(s)
    sites = []
    for key, item in overlap_dict.items():
        sites = sites+item
    sites = dict.fromkeys(list(set(sites)))
    if return_overlap_dict==True:
        return list(sites.keys()), overlap_dict
    else:
        return list(sites.keys())



def filter_sv(somvar_dict, sites, verbose=False):
    """

    """
    filtered_somvar_dict = {}
    if verbose==True:
        print(f"#filtering away {len(sites)} sites in total")
        print('sample\tbefore_after\tn_sites')
    for sample, dataframe in somvar_dict.items():
        if verbose==True:
            print(f"{sample}\tbefore_filtering\t{dataframe.shape[0]}")
        filtered_somvar_dict[sample] = pd.DataFrame([k for i,k in dataframe.iterrows() if not k['ID'] in sites])
        if verbose==True:
            print(f"{sample}\tafter_filtering\t{filtered_somvar_dict[sample].shape[0]}")
    return filtered_somvar_dict