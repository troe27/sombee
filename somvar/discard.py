import os
import pandas as pd

def load_discard_csv(sample_folder,
                     file_loc="tmpNanoSeq/post/discardedvariants.csv",
                     simple=True):
    """
    Loads a list of discarded somatic variants from a CSV file and optionally returns them in a simplified format.
    
    Parameters:
    sample_folder (str): Path to the folder containing the discard file.
    file_loc (str, optional): Relative path to the discard file within the sample folder. 
                              Defaults to "tmpNanoSeq/post/discardedvariants.csv".
    simple (bool, optional): If True, returns a list of variant site identifiers in the format "chromosome:position".
                             If False, returns the full DataFrame. Defaults to True.
    
    Returns:
    list or pandas.DataFrame:
        - If `simple=True`, returns a list of discarded variant site identifiers (e.g., ["chr1:12345", "chr2:67890"]).
        - If `simple=False`, returns a pandas DataFrame containing the full discarded variant information.
    """
    discard_file = sample_file = os.path.join(sample_folder, file_loc)
    discard = pd.read_csv(sample_file)
    discard = discard.drop_duplicates(subset=["chrom", "chromStart"])
    discard['POS'] = discard['chromStart']+1
    discard["ID"] = [k['chrom']+':'+str(discard.POS) for _,k in discard.iterrows()]
    
    if simple==False:
        return discard
    else:
        return list(discard["ID"])


def load_discard_list(folder="./test_conservative_noTGCA", 
                      prefix="test_out",
                      simple=True,
                      verbose=True):
    """
    Aggregates discarded somatic variant site identifiers from multiple files within a folder.

    Parameters:
    folder (str, optional): Path to the directory containing the discard files. Defaults to "./test_conservative_noTGCA".
    prefix (str, optional): Prefix used to filter relevant discard files in the folder. Defaults to "test_out".
    simple (bool, optional): Passed to `somvar_load_discards()`. If True, returns a list of variant site identifiers.
                             If False, returns full DataFrames. Defaults to True.

    Returns:
    list: A unique list of discarded somatic variant site identifiers (e.g., ["chr1:12345", "chr2:67890"]).
    """
    sv_discard = []
    for i in [j for j in  os.listdir(folder) if j.startswith(prefix)]:
        l = os.path.join(folder, i)
        if verbose==True:
            print(f'loading data from {l}')
        sv_discard = sv_discard + load_discard_csv(sample_folder=l)
        if verbose==True:
            print('done!')
    sv_discard = list(set(sv_discard))
    return sv_discard