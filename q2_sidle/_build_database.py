import copy
import itertools as it
import os
import warnings

import biom
import dask
from dask.delayed import Delayed
import dask.dataframe as dd
import numpy as np
import pandas as pd

from dask.delayed import Delayed
from qiime2 import Metadata, Artifact
from qiime2.plugin import ValidationError
from q2_sidle._utils import (_setup_dask_client, 
                             degen_reps,
                             )                             

def reconstruct_database(
    region: str,
    regional_alignment: pd.DataFrame,
    kmer_map: Delayed,
    blocksize: int=10000,
    debug: bool=False, 
    n_workers: int=1,
    client_address: str=None,
    ) -> (pd.DataFrame):
    """
    Reconstructs a kmer database based on the aligned regions

    Parameters
    ----------
    region
        ...
    regional_alignment
        ...
    kmer_map
        ...
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avalaibel resources.
    client_address: str
        The IP address for an existing dask client/cluster
    """    
    # Sets up the client
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers, address=client_address)
    
    # Adresses the regional assignments
    region_order, region_names, num_regions = _check_regions(region)
    
    # Imports the alignment maps and gets the kmers that were aligned
    align_map = pd.concat(
        axis=0, 
        sort=False, 
        objs=regional_alignment)
    align_map.drop_duplicates(['asv', 'kmer'], inplace=True)
    align_map.replace({'region': region_order}, inplace=True)
    aligned_kmers = _get_unique_kmers(align_map['kmer'])
    print('aligment map loaded')
    
    ### Sets up the kmer database
    # Pulls out the kmers that are of interest based on their aligned sequences
    kmer_alignments = [
        dask.delayed(_filter_to_aligned)(kmer, aligned_kmers) 
        for kmer in kmer_map
    ]
    # Gets the full list of sequences to be aligned and subsets the kmer map 
    # to retain only these sequences
    db_seqs = [_pull_unique(df) for df in kmer_alignments]
    db_seqs = np.unique(np.hstack(dask.compute(*db_seqs)))
    print('database kmers identified')
    kmer_map = dd.from_delayed(
        dfs=[dask.delayed(_build_region_db)(df, db_seqs) 
             for df in kmer_alignments],
        meta=[('db-seq', 'str'), ('region', 'str'), ('kmer', 'str')]
    )
    # Cleans up the kmer dataframe by resizing it, and indexing it by the databasae
    # sequence and region
    npartitons = int(np.ceil(kmer_map.shape[0].compute() / blocksize))
    kmer_map = kmer_map.repartition(npartitons)
    kmer_map = kmer_map.replace({'region': region_order})
    kmer_map['kmer'] = kmer_map['kmer'].apply(
        _clean_kmer_list, 
        meta=('kmer', 'str')
    )
    kmer_map['kmer'] = kmer_map['kmer'].apply(
        _check_db_list, 
        meta=('kmer', 'str'), 
        ref_seqs=set(db_seqs)
    )
    db_seqs = kmer_map['db-seq'].unique().compute()
    print('kmer map built')
    
    ### Combines the sequences to build a single regional database
    # Creates a new label for grouping magic
    kmer_map['composite'] = kmer_map.apply(
        lambda x: '{}-{}'.format(x['db-seq'], x['region']), 
        axis=1,
        meta=(None, 'object')
    )
    kmer_map['kmer'] = kmer_map['kmer'].apply(
        _clean_kmer_list, meta=('kmer', 'str')
    )
    kmer_map['kmer'] = kmer_map['kmer'].apply(
        _check_db_list, meta=('kmer', 'str'), ref_seqs=set(db_seqs)
    )


    region_db = kmer_map.groupby('composite')['kmer'].apply(
        _get_regional_seqs,
        meta=('kmer', 'str')
    ).reset_index()
    region_db.columns = ['composite', 'kmer']
    region_db['region'] = region_db['composite'].apply(
        lambda x: int(x.split("-")[1]), 
        meta=(None, int))
    region_db['db-seq'] = region_db['composite'].apply(
        lambda x: x.split("-")[0], meta=(None, str)
    )
    print('regions grouped')
    
    ### Starts untangling the database
    # Pulls out the intersetion of the existing sequences and finds the current
    # kmer matches
    print('tidying round 1')
    shared_kmers_no1, shared_check_no1 = \
        _check_intersection(region_db, set('X'))
    shared_check_no1 = shared_check_no1.compute()
    # Splits the data into tidy and untidy kmers
    tidy_kmers_no1 = \
        shared_check_no1.loc[shared_check_no1['tidy'], 'kmer'].unique()
    untidy_kmers_no1 = \
        shared_check_no1.loc[~shared_check_no1['tidy'], 'kmer'].unique()
    # Gets tidy and untidy sequences
    tidy_kmer_no1 = shared_check_no1.loc[shared_check_no1['tidy'], 'kmer'].unique()
    tidy_seqs_no1 = \
        shared_kmers_no1.loc[shared_kmers_no1['kmer'].isin(tidy_kmer_no1)]
    tidy_seqs_no1 = tidy_seqs_no1['db-seq'].compute()
    # Starts building the list of tidied database maps
    tidy_maps = [shared_kmers_no1.loc[shared_kmers_no1['db-seq'].isin(tidy_seqs_no1)]]
    print('round 1 tidied')
    
    # Pulls out the second round of umatched kmers and filers againt the
    # already assigned sequences
    print('tidying round 2')
    region_db2 = region_db.loc[~region_db['db-seq'].isin(tidy_seqs_no1)].copy()
    shared_kmers_no2, shared_check_no2 = \
        _check_intersection(region_db2, set('X') | set(tidy_seqs_no1))
    shared_check_no2 = shared_check_no2.compute()
    tidy_kmer_no2 = shared_check_no2.loc[shared_check_no2['tidy'], 'kmer'].unique()
    tidy_seqs_no2 = \
        shared_kmers_no2.loc[shared_kmers_no2['kmer'].isin(tidy_kmer_no1)]
    tidy_seqs_no2 = tidy_seqs_no2['db-seq'].compute()
    tidy_maps.append(shared_kmers_no2.loc[shared_kmers_no2['db-seq'].isin(tidy_seqs_no2)])
    print('round 2 tidied')
    
    # Pulls out the remaining sequences and preps them to be filtered
    print('tidying round 3')
    unkempt = region_db2.loc[~region_db2['db-seq'].isin(tidy_seqs_no2)]
    unkempt = unkempt.groupby('db-seq').apply(
        _define_shared, 
        matched=(set('X') | set(tidy_seqs_no1) | set(tidy_seqs_no2)),
        meta=('kmer', 'str'),
    )
    unkempt = unkempt.compute()
    if len(unkempt) > 0:
        to_map = unkempt.apply(lambda x: pd.Series(x.split("|"))).unstack().dropna()
        to_map = to_map.reset_index()
        to_map.columns = ['counter', 'db-seq', 'clean_name']
        detangled = _detangle_names(to_map.copy())
        detangled = detangled.reset_index()
        detangled.rename(columns={'clean_name': 'kmer'}, inplace=True)
    else:
        detangled = pd.Series([], 
                             index=pd.Index([], name='db-seq'),  
                             name='kmer')
        detangled.reset_index(inplace=True)
    print('round 3 tidied')
        
    # Combines the earlier maps
    print('combining solutions')
    combined = pd.concat(axis=0, objs=[
        dd.concat(dfs=tidy_maps, axis=0).compute(),
        detangled
    ])
    combined.rename(columns={'kmer': 'clean_name'}, inplace=True)
    combined.set_index('db-seq', inplace=True)

    return combined


def _build_region_db(df, kmers):
    df = df.loc[df['db-seq'].isin(kmers)].copy()    
    return df[['db-seq', 'region', 'kmer']]


def _check_db_list(x, ref_seqs):
    """
    Filters out sequences which aren't in the database designation
    """
    f_ = lambda x: x if x in ref_seqs else 'X'
    kmer = '|'.join(np.unique([f_(y) for y in x.split('|')]))
    
    return kmer
    

def _check_intersection_delayed(df, matched=set('X')):
    """
    Identifies the overlap in sequences and checks the quality of updated kmers
    """
    # Gets the combined, cleaned kmers
    shared_kmers = df.groupby('db-seq').apply(
        _define_shared, 
        matched=matched,
        meta=('kmer', 'str')
    )
    shared_kmers = shared_kmers.reset_index()
    shared_kmers.columns = ['db-seq', 'kmer']
 
    # Checks the quality fo the kmer mapping
    shared_check = shared_kmers.groupby('kmer')['db-seq'].apply(
        lambda x: '|'.join(sorted(x)), 
        meta=('db-seq', str)
    ).reset_index()
    shared_check.columns = ['kmer', 'shared-name']
    shared_check['tidy'] = shared_check['kmer'] == shared_check['shared-name']
    
    return shared_kmers, shared_check
    
    
def _check_regions(region):
    """
    Converts the region to a numeric assignment
    """
    region, region_idx = np.unique(region, return_index=True)
    region_order = {region: i for (i, region) in zip(*(region_idx, region))}
    region_names = {i: r for r, i in region_order.items()}
    num_regions = len(region_order)
    
    return region_order, region_names, num_regions


def _clean_kmer_list(x):
    """
    Filters degenerates from the kmer list
    """
    kmer = '|'.join(np.unique([y.split('@')[0] for y in x.split("|")]))
    return kmer


def _define_shared(g, matched):
    regional_sets = g['kmer'].apply(lambda x: set(x.split("|"))).values
    intersect = set.intersection(*regional_sets) - matched
    cleaned = '|'.join(sorted(intersect))
    return cleaned


def _detangle_names(long_):
    """
    Splits up the mapping into clearer clusters

    Parameters
    ----------
    long_: DataFrame
        A mapping between the database sequence (`db-seq`), the matched 
        sequence  (`clean_name`), and a column (`counter`) containing an 
        integer describing how many sequences are mapped to the 
        database sequence in a sorted order.
        It is important that this be a long, defense association: condensed/
        upper triangle matrices will not work here.

    Returns
    -------
    Series
        A mapping between the databset sequence and a unique group name.
    """

    long_.sort_values(['db-seq', 'clean_name'], ascending=True, inplace=True)
    long_['counter'] = 1
    long_['counter'] = long_.groupby(['db-seq'])['counter'].cumsum() - 1
    
    # gets the horizontal mapping that pivots the dataframe into a square
    # distance/mapping matrix
    map_ = \
        long_.set_index(['db-seq', 'clean_name'])['counter'].unstack().notna()
    diagnonal = (np.tril(np.ones(map_.shape)) == np.triu(np.ones(map_.shape)))
    map_ = (map_ | diagnonal)
    
    # Clusters the sequences
    map_.sort_values(list(map_.columns), inplace=True)
    map_ = map_[map_.index]
        
    # Calculates the number of connections within the cluster
    num_connections = (long_.groupby('db-seq')['counter'].max() + 1)
    long_.set_index('clean_name', inplace=True)
    long_['co-connection'] = num_connections[long_.index]
    co_min = long_.groupby('db-seq')['co-connection'].min()
    connection_penalty = (num_connections - co_min + 1)

    # Builds a mask that penalizes for the number of connections that should
    # be present
    penalty_v = (np.ones((len(connection_penalty), 1)) * 
                 connection_penalty.loc[map_.index].values)
    penalty_mask = (penalty_v == penalty_v.T)
    
    overlap = penalty_mask & map_
    overlap = overlap.mask(~overlap, np.nan)
    overlap = overlap.unstack().dropna().reset_index()

    # # In some cases, there are sequences that have the same number of 
    # # connections but don't  connect to the same sequences. So,for pairs 
    # # of unique sequences which are initally mapped, we want to check the
    # # number of sequences that overlap
    shared = \
        overlap.groupby('db-seq')['clean_name'].apply(lambda x: set(x.values))
    def check_shared(x):
        seqs0 = shared.loc[x['clean_name']]
        seqs1 = shared.loc[x['db-seq']]
        if len(seqs0) != len(seqs1):
            return 0
        else:
            return len(seqs0 & seqs1) - 1

    check_diff = overlap['clean_name'] != overlap['db-seq']

    if check_diff.any():
        overlap['seq_count'] = overlap.loc[check_diff].apply(check_shared, axis=1)
    else:
        overlap['seq_count'] = np.nan
    overlap['matches'] = (overlap['seq_count'] == overlap[0]) | ~check_diff

    overlap.sort_values(['db-seq', 'clean_name'], 
                        inplace=True, 
                        ascending=True)
    overlap.drop_duplicates(['db-seq', 'clean_name'], 
                            inplace=True)
    new_name = \
        overlap.groupby('db-seq')['clean_name'].apply(lambda x: "|".join(x))

    return new_name


def _filter_to_aligned(df, aligned_kmers):
    """
    Filters a delayed datafrme to the aligned dataset
    """
    df =  df.loc[df.index.isin(aligned_kmers)].copy()
    return df.reset_index()[['db-seq', 'region', 'kmer']]


def _get_clean(df):
    clean_kmers = \
        df.groupby(['db-seq', 'region'])['value'].apply(
            lambda x: "|".join(x.values))
    return clean_kmers.reset_index()


def _get_regional_seqs(x):
    """
    Gets the sequences within a region
    """
    degenerates = np.unique(np.hstack([y.split("|") for y in x]))
    return '|'.join(degenerates)


def _get_unique_kmers(series):
    """
    Extracts the unique database sequence from a set of kmers
    """
    kmers = np.hstack([[a.split("@")[0] for a in kmer.split('|')] 
                       for kmer in series])
    return np.sort(np.unique(kmers))


@dask.delayed
def _pull_unique(df):
    """
    Helper function for `_get_unique_kmers`
    """
    return _get_unique_kmers(df['kmer'])

