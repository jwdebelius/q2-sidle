# Real mapping

import biom
import numpy as np
import pandas as pd
import skbio
from skbio import DNA

from qiime2 import Artifact, Metadata

coverage = pd.DataFrame(
    data=np.array([
        list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'), 
        list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'),
        list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'),
        list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'),
        list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'), 
        list('000000000000000000000000000000000000000000000000000011111111111111100000000000000000'),
        list('000000000000000000000000000000000000000000000000000011111111111111100000000000000000'),
        list('000000000000000000000000000000000000000000000000000011111111111111100000000000000000'),
        list('111111111111111000000000000000000000000000000000000000000000000000000000000000000000'), 
        list('000000000000000000000000000011111111111111100000000000000000000000000000000000000000'),
    ]).astype(int),
    index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                    'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                    name='feature-id')
    )

extra_alignment = Artifact.import_data('FeatureData[AlignedSequence]', pd.Series(
    data=[
        skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTAAGAGAGGCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq1'}),
        skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTCAGAGTTTCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq2'}),
        skbio.DNA('TATGGTACTCATWTCCGCGTTGGAGTTATGATGATGGGGTGA-CACCTCGTTCCAGTTCCGCGCTTCATGACGTGCTGACC---', metadata={'id': 'seq3'}),
        skbio.DNA('------------------------------AAGGCGGGTGAG-CACCTCGTCCCGGAGACGAGAGGCATGACGTG---ATCCGT', metadata={'id': 'seq4'}),
        skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGTGGATGTAGAGCCACCTGACGTGC--ACCTG-', metadata={'id': 'seq5'}),
        skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGAAAATGTAGAGCCACCTGACGTGC--ACC---', metadata={'id': 'seq6'}),
        skbio.DNA('------------GCGAAGCGGCTCAGG---------------------------------------------------------', metadata={'id': 'asv01'}),
        skbio.DNA('------------ATCCGCGTTGGAGTT---------------------------------------------------------', metadata={'id': 'asv02'}),
        skbio.DNA('------------TTCCGCGTTGGAGTT---------------------------------------------------------', metadata={'id': 'asv03'}),
        skbio.DNA('------------CGTTTATGTATGCCC---------------------------------------------------------', metadata={'id': 'asv04'}),
        skbio.DNA('------------CGTTTATGTATGCCT---------------------------------------------------------', metadata={'id': 'asv05'}),
        skbio.DNA('----------------------------------------------------AGAGAGGCTGAATCC-----------------', metadata={'id': 'asv06'}),
        skbio.DNA('----------------------------------------------------AGAGTTTCTGAATCC-----------------', metadata={'id': 'asv07'}),
        skbio.DNA('----------------------------------------------------CCGGAGACGAGAGGC-----------------', metadata={'id': 'asv08'}),
        skbio.DNA('AAGGCTAGTCATCGT---------------------------------------------------------------------', metadata={'id': 'asv09'}),
        skbio.DNA('----------------------------TGAGTWWGGGAGGGA-----------------------------------------', metadata={'id': 'asv10'}),
        ],
    index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6', 'asv01', 'asv02',
                    'asv03', 'asv04', 'asv05', 'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                    dtype='object')
    )
)

db_alignment = Artifact.import_data('FeatureData[AlignedSequence]', pd.Series({
    'seq1': skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTAAGAGAGGCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq1'}),
    'seq2': skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTCAGAGTTTCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq2'}),
    'seq3': skbio.DNA('TATGGTACTCATWTCCGCGTTGGAGTTATGATGATGGGGTGA-CACCTCGTTCCAGTTCCGCGCTTCATGACGTGCTGACC---', metadata={'id': 'seq3'}),
    'seq4': skbio.DNA('------------------------------AAGGCGGGTGAG-CACCTCGTCCCGGAGACGAGAGGCATGACGTG---ATCCGT', metadata={'id': 'seq4'}),
    'seq5': skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGTGGATGTAGAGCCACCTGACGTGC--ACCTG-', metadata={'id': 'seq5'}),
    'seq6': skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGAAAATGTAGAGCCACCTGACGTGC--ACC---', metadata={'id': 'seq6'}),
    }))
db_sequences =  Artifact.import_data('FeatureData[Sequence]', pd.Series({
    'seq1': skbio.DNA('AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTAAGAGAGGCTGAATCCATGACGTGACCAGC', metadata={'id': 'seq1'}),
    'seq2': skbio.DNA('AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTCAGAGTTTCTGAATCCATGACGTGACCAGC', metadata={'id': 'seq2'}),
    'seq3': skbio.DNA('TATGGTACTCATWTCCGCGTTGGAGTTATGATGATGGGGTGACACCTCGTTCCAGTTCCGCGCTTCATGACGTGCTGACC', metadata={'id': 'seq3'}),
    'seq4': skbio.DNA('AAGGCGGGTGAGCACCTCGTCCCGGAGACGAGAGGCATGACGTGATCCGT', metadata={'id': 'seq4'}),
    'seq5': skbio.DNA('AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGTGGATGTAGAGCCACCTGACGTGCACCTG', metadata={'id': 'seq5'}),
    'seq6': skbio.DNA('AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGAAAATGTAGAGCCACCTGACGTGCACC', metadata={'id': 'seq6'}),
    }))

region1_db_seqs = Artifact.import_data('FeatureData[Sequence]', pd.Series({
    'seq1|seq2': DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq1|seq2'}),
    'seq3@0001': DNA('ATCCGCGTTGGAGTT',  metadata={'id': 'seq3@0001'}),
    'seq3@0002': DNA('TTCCGCGTTGGAGTT', metadata={'id': 'seq3@0002'}),
    'seq5': DNA('CGTTTATGTATGCCC', metadata={'id': 'seq5'}),
    'seq6': DNA('CGTTTATGTATGCCT', metadata={'id': 'seq6'}), 
    }))
region2_db_seqs = Artifact.import_data('FeatureData[Sequence]', pd.Series({
    'seq1': DNA('AGAGAGGCTGAATCC', metadata={'id': 'seq1'}),
    'seq2': DNA('AGAGTTTCTGAATCC', metadata={'id': 'seq2'}),
    'seq3': DNA('CCAGTTCCGCGCTTC', metadata={'id': 'seq3'}),
    'seq4': DNA('CCGGAGACGAGAGGC', metadata={'id': 'seq4'}),
    'seq5': DNA('TGGATGTAGAGCCAC', metadata={'id': 'seq5'}),
    'seq6': DNA('AAAATGTAGAGCCAC', metadata={'id': 'seq6'}),
    }))
region1_db_map = Artifact.import_data('FeatureData[KmerMap]', pd.DataFrame(
    data=np.array([['seq1', 'seq1|seq2', 'Bludhaven', 'WANTCAT', 'ATGATGATG', 13, 28, 15],
                   ['seq2', 'seq1|seq2', 'Bludhaven', 'WANTCAT', 'ATGATGATG', 13, 28, 15],
                   ['seq3@0001', 'seq3@0001', 'Bludhaven', 'WANTCAT', 'ATGATGATG',  13, 28, 15],
                   ['seq3@0002', 'seq3@0002', 'Bludhaven', 'WANTCAT', 'ATGATGATG',  13, 28, 15],
                   ['seq5', 'seq5', 'Bludhaven', 'WANTCAT', 'ATGATGATG',  13, 28, 15],
                   ['seq6', 'seq6', 'Bludhaven', 'WANTCAT', 'ATGATGATG',  13, 28, 15],
                   ]),
    index=pd.Index(['seq1', 'seq2', 'seq3', 'seq3', 'seq5', 'seq6'], 
                   name='db-seq'),
    columns=['seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer', 'fwd-pos', 'rev-pos', 'kmer-length']
    ))
region2_db_map = Artifact.import_data('FeatureData[KmerMap]', pd.DataFrame(
    data=np.array([['seq1', 'seq1', 'Gotham', 'CACCTCGTN', 'MTGACGTG', 53, 67, 15],
                   ['seq2', 'seq2', 'Gotham', 'CACCTCGTN', 'MTGACGTG', 53, 67, 15],
                   ['seq3', 'seq3', 'Gotham', 'CACCTCGTN', 'MTGACGTG', 53, 67, 15],
                   ['seq4', 'seq4', 'Gotham', 'CACCTCGTN', 'MTGACGTG', 53, 67, 15],
                   ['seq5', 'seq5', 'Gotham', 'CACCTCGTN', 'MTGACGTG', 53, 67, 15],
                   ['seq6', 'seq6', 'Gotham', 'CACCTCGTN', 'MTGACGTG', 53, 67, 15],
                   ]),
    index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
                    name='db-seq'),
    columns=['seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer', 'fwd-pos', 'rev-pos',  'kmer-length'],
    ))
region1_rep_seqs = Artifact.import_data('FeatureData[Sequence]', pd.Series({
    'asv01': DNA('GCGAAGCGGCTCAGG', metadata={'id': 'asv01'}),
    'asv02': DNA('ATCCGCGTTGGAGTT', metadata={'id': 'asv02'}),
    'asv03': DNA('TTCCGCGTTGGAGTT', metadata={'id': 'asv03'}),
    'asv04': DNA('CGTTTATGTATGCCC', metadata={'id': 'asv04'}),
    'asv05': DNA('CGTTTATGTATGCCT', metadata={'id': 'asv05'}),
    }))
region1_align = \
    Artifact.import_data('FeatureData[KmerAlignment]', pd.DataFrame(
        data=np.array([['seq1|seq2', 'asv01', 15, 0, 2, 'Bludhaven'],
                       ['seq3@0001', 'asv02', 15, 0, 2, 'Bludhaven'],
                       ['seq3@0001', 'asv03', 15, 1, 2, 'Bludhaven'],
                       ['seq3@0002', 'asv02', 15, 1, 2, 'Bludhaven'],
                       ['seq3@0002', 'asv03', 15, 0, 2, 'Bludhaven'],
                       ['seq5', 'asv04', 15, 0, 2, 'Bludhaven'],
                       ['seq5', 'asv05', 15, 1, 2, 'Bludhaven'],
                       ['seq6', 'asv04', 15, 1, 2, 'Bludhaven'],
                       ['seq6', 'asv05', 15, 0, 2, 'Bludhaven']],
                       dtype=object),
        columns=['kmer', 'asv',  'length', 'mismatch', 'max-mismatch', 'region'],
    ))
region2_align = \
    Artifact.import_data('FeatureData[KmerAlignment]', pd.DataFrame(
        data=np.array([['seq1', 'asv06', 15, 0, 2, 'Gotham'],
                       ['seq2', 'asv07', 15, 0, 2, 'Gotham'],
                       ['seq3', 'asv08', 15, 0, 2, 'Gotham'],
                       ['seq4', 'asv09', 15, 0, 2, 'Gotham'],
                       ['seq5', 'asv10', 15, 0, 2, 'Gotham'],
                       ['seq6', 'asv11', 15, 0, 2, 'Gotham']],
                       dtype=object),
        columns=['kmer', 'asv',  'length', 'mismatch', 'max-mismatch', 'region'],
    ))
region1_counts = Artifact.import_data('FeatureTable[Frequency]', biom.Table(
    np.array([[150,   0,   0,  50, 50],
              [125,  50,  50,  25, 25],
              [100,   0, 100,  50, 50]]).T,
    sample_ids=['sample1', 'sample2', 'sample3'],
    observation_ids=['asv01', 'asv02', 'asv03', 'asv04', 'asv05'],
    ))
region1_alt_counts = Artifact.import_data('FeatureTable[Frequency]', biom.Table(
    np.array([[150,   0,   0,  50, 50, 50],
              [125,  50,  50,  25, 25, 25],
              [100,   0, 100,  50, 50, 0]]).T,
    sample_ids=['sample1', 'sample2', 'sample3'],
    observation_ids=['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 'asv20'],
    ))
region2_counts = Artifact.import_data('FeatureTable[Frequency]', biom.Table(
    data=np.array([[100,  50,    0,  50,  50, 50],
                   [100,  25,  100,  25,  25, 25],
                   [  0, 100,  100,   0,  50, 50]]).T,
    sample_ids=['sample1', 'sample2', 'sample3'],
    observation_ids=['asv06', 'asv07', 'asv08', 'asv09', 'asv10', 'asv11'],
    ))

taxonomy_gg = pd.Series(
    data=np.array(['k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batgirl; g__Gordon; s__Barbara',
                   'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batgirl; g__Gordon; s__Barbara',
                   'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batman; g__Wayne; s__',
                   'k__DCU; p__Superhero; c__Gotham; o__Civillian; f__Doctor; g__Thompson; s__Leslie',
                   'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Grayson; s__Dick',
                   'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Todd; s__Jason',
                   ]),
    index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], name='Feature ID'),
    name='Taxon',
    )
taxonomy =  Artifact.import_data('FeatureData[Taxonomy]', taxonomy_gg)

seq_map = pd.DataFrame(
            data=np.array([['seq1', 0, 'WANTCAT', 13,  1, 'CACCTCGTN', 53, 15],
                           ['seq2', 0, 'WANTCAT', 13, 1, 'CACCTCGTN', 53, 15],
                           ['seq3', 0, 'WANTCAT', 13, 1, 'CACCTCGTN', 53, 15],
                           ['seq4', 1, 'CACCTCGTN', 53, 1, 'CACCTCGTN', 53, 15],
                           ['seq5', 0, 'WANTCAT', 13, 1, 'CACCTCGTN', 53, 15],
                           ['seq6', 0, 'WANTCAT', 13, 1, 'CACCTCGTN', 53, 15],
                           ]),
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], name='db-seq'),
            columns=['clean_name', 'first-region', 'first-fwd-primer', 'first-fwd-pos',
                     'last-region', 'last-fwd-primer', 'last-fwd-pos', 'last-kmer-length']
            )
seq_map = Artifact.import_data('FeatureData[SidleReconstruction]', seq_map)
    
db_summary = pd.DataFrame.from_dict(orient='index', data={
    'seq1': {'num-regions': 2, 
             'total-kmers-mapped': 2, 
             'mean-kmer-per-region': 1.,
             'stdv-kmer-per-region': 0.,
             'mapped-asvs': 'asv01|asv06'
            },
    'seq2': {'num-regions': 2, 
             'total-kmers-mapped': 2, 
             'mean-kmer-per-region': 1,
             'stdv-kmer-per-region': 0,
             'mapped-asvs': 'asv01|asv07',
            },
    'seq3': {'num-regions': 2, 
             'total-kmers-mapped': 3, 
             'mean-kmer-per-region': 1.5,
             'stdv-kmer-per-region': np.std([1, 2], ddof=1),
             'mapped-asvs': 'asv02|asv03|asv08'
            },
    'seq4': {'num-regions': 1, 
             'total-kmers-mapped': 1, 
             'mean-kmer-per-region': 1,
             'stdv-kmer-per-region': 0,
             'mapped-asvs': 'asv09'
            },
    'seq5': {'num-regions': 2, 
             'total-kmers-mapped': 2, 
             'mean-kmer-per-region': 1,
             'stdv-kmer-per-region': 0,
             'mapped-asvs': 'asv04|asv05|asv10',
            },
    'seq6': {'num-regions': 2, 
             'total-kmers-mapped': 2, 
             'mean-kmer-per-region': 1,
             'stdv-kmer-per-region': 0,
             'mapped-asvs': 'asv04|asv05|asv11',
            },
    })
db_summary.index.set_names('feature-id', inplace=True)   
db_summary = Artifact.import_data('FeatureData[ReconstructionSummary]', Metadata(db_summary))
