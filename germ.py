import pandas as pd
import os
from cyvcf2 import VCF

from . import base

def filter_germ(somvar_dict, vcf_file, threads=1):
    """
    Filters somatic variant calls based on a given VCF file of germline variants.

    This function compares somatic variant sites in `somvar_dict` against sites in 
    the provided VCF file (`vcf_file`). Any somatic variants that match a site in 
    the VCF file are classified as germline and filtered out. The function returns 
    the set of germline sites, a list of filtered germline sites per sample, and 
    an updated somatic variant dictionary with germline sites removed.

    Args:
        somvar_dict (dict): A dictionary where keys are sample identifiers and 
                            values are lists of somatic variant strings in the 
                            format "CHROM,POS,...".
        vcf_file (str): Path to the VCF file containing germline variant sites.
        threads (int, optional): Number of threads to use when reading the VCF file. 
                                 Defaults to 1.

    Returns:
        tuple: A tuple containing:
            - IDs (dict): A dictionary of germline variant sites from the VCF file, 
                          with keys in the format "CHROM:POS".
            - filtered_germ_sites (list): A list of lists, where each sublist contains 
                                          a sample identifier and a filtered germline 
                                          site in "CHROM:POS" format.
            - filtered_somvar_dict (dict): A dictionary with the same structure as 
                                           `somvar_dict`, but with germline sites 
                                           removed.

    """
    vcf = VCF(vcf_file, threads=threads)
    IDs = dict.fromkeys([str(v.CHROM)+':'+str(v.POS) for v in vcf])
    filtered_somvar_dict = {}
    filtered_germ_sites = []
    for key, item in somvar_dict.items():
        somID = item["ID"]
        sites = [i for i in somID if i in IDs]
        filtered_germ_sites =  filtered_germ_sites + [[key, i] for i in sites]
    filtered_somvar_dict = base.filter_somvar(somvar_dict=somvr_dict,
                                         sites=filtered_germ_sites)
    return  IDs, filtered_germ_sites, filtered_somvar_dict