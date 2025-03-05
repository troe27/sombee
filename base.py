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
        svs["ID"] = [str(k.chrom)+':'+str(k.chromStart) for i,k in svs.iterrows()]
        somvar_dict[i]=svs
    return somvar_dict


def filter_somvar(somvar_dict, sites):
    """

    """
    filtered_somvar_dict = {}
    for sample, dataframe in somvar_dict.items():
        filtered_somvar_dict[sample] = pd.DataFrame([k for i,k in dataframe.iterrows() if not k['ID'] in sites])
    return filtered_somvar_dict