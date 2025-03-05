import os
def filter_discard(somvar_dict,discard_list):
    """
    Filters out specific somatic variant sites from a given dictionary based on a discard list.
    
    Parameters:
    somvar_dict (dict): A dictionary where keys are sample identifiers and values are lists of 
                        variant data strings in the format "chromosome,position,additional_info".
    discard_list (list): A list of variant site identifiers in the format "chromosome:position" 
                         that should be removed.
    
    Returns:
    tuple: 
        - discarded_somvar (list): A list of discarded variants, each represented as a [sample_id, "chromosome:position"] pair.
        - filtered_somvar_dict (dict): A dictionary similar to `somvar_dict` but with the specified 
                                       variant sites removed.
    """
    discardIDs = discard_list
    filtered_somvar_dict = {}
    discarded_somvar = []
    for key, item in somvar_dict.items():
        somID = [i.split(",")[0]+':'+i.split(",")[1] for i in item]
        sites = [i for i in somID if i in discardIDs]
        
        filt_somvar = [i for i in item if i.split(",")[0]+':'+i.split(",")[1] not in discardIDs ]
        filtered_somvar_dict[key] = filt_somvar
        discarded_somvar =  discarded_somvar + [[key, i] for i in sites]
    return discarded_somvar, filtered_somvar_dict


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
    if simple==False:
        return discard
    else:
        return [str(k.chrom)+':'+str(k.chromStart) for i,k in discard.iterrows()]


def load_discard_list(folder="./test_conservative_noTGCA", 
                      prefix="test_out",
                      simple=True):
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
        sv_discard = sv_discard + load_discard_csv(sample_folder=l)
    sv_discard = list(set(sv_discard))
    return sv_discard