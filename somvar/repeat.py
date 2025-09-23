import pandas as pd
import numpy as np
from . import base
def filter_by_bed(bed, somvar_dict,sv_fmt="dict", ret_val='filtered_list'):
    """

    """
    
    if ret_val not in ['filtered_n', 'filtered_list', 'full']:
        raise IOError('ret_val needs to be one of the following: "filtered_n", "filtered_list" or "full".')

        if sv_fmt not in ['dict', 'df',]:
            raise IOError('sv_fmt needs to be one of the following: "dict" or "df".')
    
    if isinstance(bed, str)==True:
        bdata = pd.read_csv(bed, sep=' ', skipinitialspace=True)
    else:
        bdata=bed
        bdata.columns=['CHROM', 'POS1', 'POS2']
    if sv_fmt =="dict":
        svdf = pd.concat([item for _,item in somvar_dict.items()])
        svdf.drop_duplicates(subset='pyrkey',inplace=True)
        svdf.rename(columns={'chrom':'CHROM'}, inplace=True)
    else:
        svdf = somvar_dict
    vars = []
    for chromosome in svdf.CHROM.drop_duplicates():
        intervals_ss = bdata.loc[bdata['CHROM']==chromosome]
        variants_ss = svdf.loc[svdf.CHROM==chromosome]
        for _, variant in variants_ss.iterrows():
            POS_var=variant.POS
            interval_match = intervals_ss.loc[intervals_ss.POS1<=POS_var].loc[intervals_ss.POS2>=POS_var]
            if len(interval_match)>0:
                vars.append(list(variant)+[True]+list(interval_match.iloc[0,:]))
            else:
                vars.append(list(variant)+[False]+[np.nan for i in range(intervals_ss.shape[1])])
    vardf = pd.DataFrame(vars , columns=list(variants_ss.columns)+['filtered']+['interval_'+i for i in intervals_ss.columns])
    
    if ret_val=='filtered_n':
        return vardf['filtered'].sum()
    elif ret_val=='filtered_list':
        return list(vardf.loc[vardf.filtered==True]['ID'])
    else:
        return vardf


def filter_repeats(repeatfile, somvar_dict, return_val='list'):
    """

    """

    if return_val not in ['filtered_n', 'filtered_list', 'full']:
        raise IOError('return_val needs to be one of the following: "filtered_n", "filtered_list" or "full".')
    repeats = pd.read_csv(repeatfile, sep=' ', skipinitialspace=True)
    bed = repeats[['query-sequence', 'position_query_begin','position_query_end']]
    if return_val in {'filtered_n':'', 'filtered_list':'',}:
        return filter_by_bed(bed=bed, somvar_dict=somvar_dict, ret_val=return_val)
    else:
        variants_df =  filter_by_bed(bed=bed, somvar_dict=somvar_dict, ret_val=return_val)
        return pd.merge(left=variants_df, right=repeats, how='left', right_on=['query-sequence', 'position_query_begin','position_query_end'], left_on=['CHROM', 'interval_POS1', 'interval_POS2'])