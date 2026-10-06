import pandas as pd
import os

def get_double_sites(somvar_dict, verbose=False, return_only_IDs=True):
    """
    Identifies variants that share the same duplex barcode within a sample but occur at different positions.

    Parameters:
    ----------
    somvar_dict : dict
        Dictionary mapping sample names to their mutation data (DataFrame).
    verbose : bool, optional
        If True, prints debugging information for detected double sites (default: True).

    Returns:
    -------
    tuple
        (pandas.DataFrame, pandas.DataFrame)
        - First DataFrame contains unique variant IDs (sample, chromosome, position, barcode).
        - Second DataFrame contains additional information about double sites (sample, chromosome, position 1, position 2, barcode, distance).
    """
    shared_reads_ID = []
    shared_reads_info = []
    for sample, item in somvar_dict.items(): # for each sample
        for i, k in item.iterrows(): #for each variant
            for j,l in item.loc[item.chrom == k.chrom].iterrows(): # ss to chrom
                if l["dplxBarcode"] == k["dplxBarcode"]:
                    if not k.POS == l.POS:
                        if verbose==True:
                            print(sample)
                            print(k.chrom)
                            print(k.POS)
                            print(k.dplxBarcode)
            
                            print(l.chrom)
                            print(l.POS)
                            print(l.dplxBarcode)
                            print('---')
                            print(l.POS-k.POS)
                            print('---')
                        info = [sample,k.chrom, k.POS, l.POS,l.dplxBarcode, l.POS-k.POS ]
                        shared_reads_info.append(info)
                        shared_reads_ID.append([sample, k.chrom, k.POS,k.dplxBarcode])
                        shared_reads_ID.append([sample, l.chrom, l.POS,l.dplxBarcode])
    s_r_IDs = pd.DataFrame(shared_reads_ID, columns=['sample', "CHROM", 'POS', 'UMI']).drop_duplicates()
    s_r_IDs["ID"] = [str(k.CHROM)+':'+str(k.POS) for i,k in s_r_IDs.iterrows()]
    multi_bc = pd.DataFrame([k for i,k in pd.DataFrame(shared_reads_info).iterrows() if k[5]>0])
    multi_bc.columns = ["sample", "CHROM", "POS1", "POS2", "UMI", "DIST"]

    if return_only_IDs==True:
        return list(s_r_IDs['ID'].drop_duplicates())
    else:
        return s_r_IDs, multi_bc    