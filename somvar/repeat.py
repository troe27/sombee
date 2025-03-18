from intervaltree import IntervalTree
import pandas as pd

def filter_by_bed(somvar_dict, bed, buffer=0, return_sv=False):
    if isinstance(bed, str)==True:
        bdata = pd.read_csv(bed, sep=' ', skipinitialspace=True)
    else:
        bdata=bed
        bdata.columns=['CHROM', 'POS1', 'POS2']
    
    # Build the tree
    treedict = {}
    for chromosome in bdata.CHROM.drop_duplicates():
        tree = IntervalTree([it.Interval(int(k.POS1)-buffer, int(k.POS2) + 1+buffer) for i,k in bdata.loc[bdata.CHROM==chromosome].iterrows()])  # +1 to make stop inclusive
        treedict[chromosome]=tree
    # Query positions
    filt_sites = []
    filt_svdict = {}
    for sample, dataframe in svdict.items():
        filt_sv = []
        for i,k in dataframe.iterrows():
            if bool(treedict[k.chrom][int(k.chromStart)])==True: # True if position is in any interval
                filt_sites.append(k['ID'])
                filt_sv.append(k)
        filt_svdict[sample] = pd.DataFrame(filt_sv)
    if return_sv==True:
        return filt_sites, filt_svdict
    else:
        return filt_sites
        
        
def sv_filter_by_repeat(repeatfile, return_sv=False):
    repeatdf = pd.read_csv(bedfile, sep=' ', skipinitialspace=True)
    bed =  repeatdf[['query-sequence', 'position_query_begin','position_query_end']]
    if return_sv==True
        return filter_by_bed(bed=bed, return_sv=True)
    else:
        return filter_by_bed(bed=bed)