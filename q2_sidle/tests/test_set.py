# Real mapping

import biom
import numpy as np
import pandas as pd
from skbio import DNA

from qiime2 import Artifact

region1_db_seqs = Artifact.import_data('FeatureData[Sequence]', pd.Series({
    'seq1|seq2': DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq1|seq2'}),
    'seq3@0001': DNA('ATCCGCGTTGGAGTT',  metadata={'id': 'seq3@0001'}), 
    'seq3@0002': DNA('TTCCGCGTTGGAGTT', metadata={'id': 'seq3@0002'}), 
    'seq5': DNA('CGTTTATGTATGCCC', metadata={'id': 'seq5'}), 
    'seq6': DNA('CGTTTATGTATGCCT', metadata={'id': 'seq6'})
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
    data=np.array([['seq1', 'seq1|seq2', 'WANTCAT-CATCATCAT'],
                   ['seq2', 'seq1|seq2', 'WANTCAT-CATCATCAT'],
                   ['seq3@0001', 'seq3@0001', 'WANTCAT-CATCATCAT'],
                   ['seq3@0002', 'seq3@0002', 'WANTCAT-CATCATCAT'],
                   ['seq5', 'seq5', 'WANTCAT-CATCATCAT'],
                   ['seq6', 'seq6', 'WANTCAT-CATCATCAT'],
                   ]),
    index=pd.Index(['seq1', 'seq2', 'seq3', 'seq3', 'seq5', 'seq6'], 
                   name='db_seq'),
    columns=['seq_name', 'kmer', 'region'],
    ))
region2_db_map = Artifact.import_data('FeatureData[KmerMap]', pd.DataFrame(
    data=np.array([['seq1', 'seq1', 'Gotham'],
                   ['seq2', 'seq2', 'Gotham'],
                   ['seq3', 'seq3', 'Gotham'],
                   ['seq4', 'seq4', 'Gotham'],
                   ['seq5', 'seq5', 'Gotham'],
                   ['seq6', 'seq6', 'Gotham'],
                   ]),
    index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
                    name='db_seq'),
    columns=['seq_name', 'kmer', 'region']
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
        data=np.array([['seq1|seq2', 'asv01', 0, 15, 'WANTCAT-CATCATCAT'],
                       ['seq3@0001', 'asv02', 0, 15, 'WANTCAT-CATCATCAT'],
                       ['seq3@0001', 'asv03', 1, 15, 'WANTCAT-CATCATCAT'],
                       ['seq3@0002', 'asv02', 1, 15, 'WANTCAT-CATCATCAT'],
                       ['seq3@0002', 'asv03', 0, 15, 'WANTCAT-CATCATCAT'],
                       ['seq5', 'asv04', 0, 15, 'WANTCAT-CATCATCAT'],
                       ['seq5', 'asv05', 1, 15, 'WANTCAT-CATCATCAT'],
                       ['seq6', 'asv04', 1, 15, 'WANTCAT-CATCATCAT'],
                       ['seq6', 'asv05', 0, 15, 'WANTCAT-CATCATCAT']], 
                       dtype=object),
        columns=['kmer', 'asv', 'mismatch', 'length', 'region'],
    ).set_index('kmer'))
region2_align = \
    Artifact.import_data('FeatureData[KmerAlignment]', pd.DataFrame(
        data=np.array([['seq1', 'asv06', 0, 15, 'Gotham'],
                       ['seq2', 'asv07', 0, 15, 'Gotham'],
                       ['seq3', 'asv08', 0, 15, 'Gotham'],
                       ['seq4', 'asv09', 0, 15, 'Gotham'],
                       ['seq5', 'asv10', 0, 15, 'Gotham'],
                       ['seq6', 'asv11', 0, 15, 'Gotham']],
                       dtype=object),
        columns=['kmer', 'asv', 'mismatch', 'length', 'region']
    ).set_index('kmer'))
region1_counts = Artifact.import_data('FeatureTable[Frequency]', biom.Table(
    np.array([[150,   0,   0,  50, 50],
              [125,  50,  50,  25, 25],
              [100,   0, 100,  50, 50]]).T,
    sample_ids=['sample1', 'sample2', 'sample3'],
    observation_ids=['asv01', 'asv02', 'asv03', 'asv04', 'asv05'],
    ))
region2_counts = Artifact.import_data('FeatureTable[Frequency]', biom.Table(
    data=np.array([[100,  50,    0,  50,  50, 50],
                   [100,  25,  100,  25,  25, 25],
                   [  0, 100,  100,   0,  50, 50]]).T,
    sample_ids=['sample1', 'sample2', 'sample3'],
    observation_ids=['asv06', 'asv07', 'asv08', 'asv09', 'asv10', 'asv11']
    ))

# ### ASV Sequences and fastas
# rep_seq1_fasta = ('>asv01\nGCGAAGCGGCTCAGG\n\n'
#                   '>asv02\nATCCGCGTTGGAGTT\n\n'
#                   '>asv03\nTTCCGCGTTGGAGTT\n\n'
#                   '>asv04\nCGTTTATGTATGCCC\n\n'
#                   '>asv05\nCGTTTATGTATGCCT\n\n'
#                   )

# seqs_1 = pd.DataFrame(
#     data=np.array([list('GCGAAGCGGCTCAGG'), # seq1/2
#                    list('ATCCGCGTTGGAGTT'), # seq3
#                    list('TTCCGCGTTGGAGTT'), # seq3
#                    list('CGTTTATGTATGCCC'), # seq5
#                    list('CGTTTATGTATGCCT'), # seq6
#                    ]),
#     index=['asv01', 'asv02', 'asv03', 'asv04', 'asv05']
#     )

# rep_seq2_fasta = ('>asv06\nAGAGAGGCTGAATCC\n\n'
#                   '>asv07\nAGAGTTTCTGAATCC\n\n'
#                   '>asv08\nCCAGTTCCGCGCTTC\n\n'
#                   '>asv09\nCCGGAGACGAGAGGC\n\n'
#                   '>asv10\nTGGATGTAGAGCCAC\n\n'
#                   '>asv11\nAAAATGTAGAGCCAC\n\n'
#                   )

# seqs_2 = pd.DataFrame(
#     data=np.array([list('AGAGAGGCTGAATCC'), # seq1
#                    list('AGAGTTTCTGAATCC'), # seq2
#                    list('CCAGTTCCGCGCTTC'), # seq3
#                    list('CCGGAGACGAGAGGC'), # seq4
#                    list('TGGATGTAGAGCCAC'), # seq5
#                    list('AAAATGTAGAGCCAC'), # seq6
#         ]),
#     index=['asv06', 'asv07', 'asv08', 'asv09', 'asv10', 'asv11'],
#     )

# ### Alignment betweewn asvs and database
# match1 = pd.DataFrame(
#     data=np.array([['seq1 | seq2', 'asv01', 0, 15],
#                    ['seq3@0001', 'asv02', 0, 15],
#                    ['seq3@0001', 'asv03', 1, 15],
#                    ['seq3@0002', 'asv02', 1, 15],
#                    ['seq3@0002', 'asv03', 0, 15],
#                    ['seq5', 'asv04', 0, 15],
#                    ['seq5', 'asv05', 1, 15],
#                    ['seq6', 'asv04', 1, 15],
#                    ['seq6', 'asv05', 0, 15]], 
#                    dtype=object),
#     columns=['kmer', 'asv', 'mismatch', 'length'],
#     )
# match1['length'] = match1['length'].astype(int)
# match1['mismatch'] = match1['mismatch'].astype(int)
# match1['region'] = 'Bludhaven'

# match2 = pd.DataFrame(
#     data=np.array([['seq1', 'asv06', 0, 15],
#                    ['seq2', 'asv07', 0, 15],
#                    ['seq3', 'asv08', 0, 15],
#                    ['seq4', 'asv09', 0, 15],
#                    ['seq5', 'asv10', 0, 15],
#                    ['seq6', 'asv11', 0, 15]],
#                    dtype=object),
#     columns=['kmer', 'asv', 'mismatch', 'length']
#     )
# match2['length'] = match2['length'].astype(int)
# match2['mismatch'] = match2['mismatch'].astype(int)
# match2['region'] = 'Gotham'

# seq_map = pd.Series({'seq1': 'seq1',
#                      'seq2': 'seq2',
#                      'seq3': 'seq3',
#                      'seq4': 'seq4',
#                      'seq5': 'seq5',
#                      'seq6': 'seq6'}, name='clean_name')
# seq_map.index.set_names('db_seq', inplace=True)
# seq_summary = pd.DataFrame.from_dict(orient='index', data={
#     'seq1': {'num_regions': 2, 
#              'total_kmers_mapped': 2, 
#              'mean_kmer_per_region': 1,
#              'stdv_kmer_per_region': 0,
#              'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batgirl;Gordon;Barbara',
#              'mapped_asvs': 'asv01 | asv06'
#             },
#     'seq2': {'num_regions': 2, 
#              'total_kmers_mapped': 2, 
#              'mean_kmer_per_region': 1,
#             'stdv_kmer_per_region': 0,
#             'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batgirl;Gordon;Barbara',
#             'mapped_asvs': 'asv01 | asv07',
#             },
#     'seq3': {'num_regions': 2, 
#              'total_kmers_mapped': 3, 
#              'mean_kmer_per_region': 1.5,
#              'stdv_kmer_per_region': np.std([1, 2], ddof=1),
#              'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batman;Wayne;',
#              'mapped_asvs': 'asv02 | asv03 | asv08'
#             },
#     'seq4': {'num_regions': 1, 
#              'total_kmers_mapped': 1, 
#              'mean_kmer_per_region': 1,
#              'stdv_kmer_per_region': 0,
#              'taxonomy': 'DCU;Superhero;Gotham;Civillian;Doctor;Thompson;Leslie',
#              'mapped_asvs': 'asv09'
#             },
#     'seq5': {'num_regions': 2, 
#              'total_kmers_mapped': 2, 
#              'mean_kmer_per_region': 1,
#              'stdv_kmer_per_region': 0,
#              'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Robin;Grayson;Dick',
#              'mapped_asvs': 'asv04 | asv05 | asv10',
#             },
#     'seq6': {'num_regions': 2, 
#              'total_kmers_mapped': 2, 
#              'mean_kmer_per_region': 1,
#              'stdv_kmer_per_region': 0,
#              'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Robin;Todd;Jason',
#              'mapped_asvs': 'asv04 | asv05 | asv11',
#             },
#     })
# seq_summary.index.set_names('clean_name', inplace=True)

# seq_align = pd.DataFrame(
#     data=np.array([['seq1', 'asv01', 0],
#                    ['seq1', 'asv06', 0],
#                    ['seq2', 'asv01', 0],
#                    ['seq2', 'asv07', 0],
#                    ['seq3', 'asv02', 0],
#                    ['seq3', 'asv03', 0], # This is zero because we take the minimum of the degenerate match
#                    ['seq3', 'asv08', 0],
#                    ['seq4', 'asv09', 0],
#                    ['seq5', 'asv04', 0],
#                    ['seq5', 'asv05', 1], # This is not degenrate, but a single differente, so...
#                    ['seq5', 'asv10', 0],
#                    ['seq6', 'asv04', 1],
#                    ['seq6', 'asv05', 0],
#                    ['seq6', 'asv11', 0],
#                    ], dtype=object),
#     columns=['clean_name', 'asv', 'mismatch'])
# seq_align.set_index('clean_name', inplace=True)

# # True table for samples
# true = pd.DataFrame(
#     data=np.array([[100,  50,   0,  50,  50, 50],
#                    [100,  25, 100,  25,  25, 25],
#                    [  0, 100, 100,   0,  50, 50],
#                    ]) * 1.,
#     index=['sample1', 'sample2', 'sample3'],
#     columns=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], name='clean_name'),
#     ).T
# # true = true.divide(true.sum(axis=0))

# counts_1 = pd.DataFrame(
#     data=np.array([[150,   0,   0,  50, 50],
#                    [125,  50,  50,  25, 25],
#                    [100,   0, 100,  50, 50],
#                   ]),
#     index=['sample1', 'sample2', 'sample3'],
#     columns=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05'], name='asv_id'),
#     ).T
# counts_2 = pd.DataFrame(
#     data=np.array([[100,  50,    0,  50,  50, 50],
#                    [100,  25,  100,  25,  25, 25],
#                    [  0, 100,  100,   0,  50, 50]]),
#     index=['sample1', 'sample2', 'sample3'],
#     columns=pd.Index(['asv06', 'asv07', 'asv08', 'asv09', 'asv10', 'asv11'], name='asv_id'),
#     ).T


# # Amplified sequences


# ### Reference database
# seq1 = list('ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGACACCTCGTAAGAGAGGCTGAATCCTGACGTGAAC')
# # #seq2: same first region as seq1, 3 nt difference in second
# seq2 = list('ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGACACCTCGTAAGAGTTTCTGAATCCTGACGTGAAC')
# # seq3: degeneracy in first region
# seq3 = list('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAAACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCAC')
# # Skipped in 1st region
# seq4 = list('GATTTTTTTATGTTCGCATGCGGAATGATGATGCCGACACCTCGTCCCGGAGACGAGAGGCTGACGTGAGC')
# seq5 = list('CATAGTCATCGTTTATGTATGCCCATGATGATGCGAGCACCTCGTATGGATGTAGAGCCACTGACGTGCGG')
# # # seq6: 1 nt diff in region 1, 3 nt in region 2
# seq6 = list('CATAGTCATCGTTTATGTATGCCTATGATGATGCGAGCACCTCGTAAAAATGTAGAGCCACTGACGTGCGG')

# primers = [['WANTCAT', 'CATCATCAT'], ['CACCTCGTN', 'CACGTCAK']]

# full_db = pd.DataFrame(
#     data=np.array([seq1, seq2, seq3, seq4, seq5, seq6]),
#     index=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
#     )
# # Somtimes a girl has to have fun
# taxonomy = pd.DataFrame(
#     data=np.array([['DCU', 'Superhero', 'Gotham', 'Batfamily', 'Batgirl', 'Gordon', "Barbara"],
#                    ['DCU', 'Superhero', 'Gotham', 'Batfamily', 'Batgirl', 'Gordon', 'Barbara'],
#                    ['DCU', 'Superhero', 'Gotham', 'Batfamily', 'Batman', 'Wayne', ''],
#                    ['DCU', 'Superhero', 'Gotham', 'Civillian', 'Doctor', 'Thompson', "Leslie"],
#                    ['DCU', 'Superhero', 'Gotham', 'Batfamily', 'Robin', 'Grayson', "Dick"],
#                    ['DCU', 'Superhero', 'Gotham', 'Batfamily', 'Robin', 'Todd', "Jason"],
#                    ]),
#     columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
#     index=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
#     )