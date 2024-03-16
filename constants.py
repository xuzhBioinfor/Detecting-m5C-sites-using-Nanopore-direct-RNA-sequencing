r"""
This module is a collection of builder functions used during training and inference for building Python object from configuration files
"""
from itertools import product
import numpy as np
import pandas as pd
import os
NUM_NEIGHBORING_FEATURES = 1
print(os.getcwd())
kmer_file = os.getcwd() + 'model/kmer_list/KMER_LIST'
CENTER_MOTIFS = pd.read_csv(kmer_file)
CENTER_MOTIFS = list(CENTER_MOTIFS['KMER'])
#print(len(CENTER_MOTIFS))
FLANKING_MOTIFS = [['G', 'A', 'C','T'] for i in range(NUM_NEIGHBORING_FEATURES)]
FLANKING_MOTIFS = list(["".join(x) for x in product(*(FLANKING_MOTIFS))])
ALL_KMERS_seq = []
for flank_motif in FLANKING_MOTIFS:
    for center_motif in CENTER_MOTIFS:
        ALL_KMERS_seq.append(flank_motif+center_motif+flank_motif)
ALL_KMERS_seq = np.unique(ALL_KMERS_seq)

ALL_KMERS = []
for seq in ALL_KMERS_seq:
    seq_list = []
    for index in range(len(seq)-4):
        seq_list.append(seq[index:index+5])
        ALL_KMERS = ALL_KMERS + seq_list

ALL_KMERS = np.unique(ALL_KMERS)
print(len(ALL_KMERS))
KMER_TO_INT = {ALL_KMERS[i]: i for i in range(len(ALL_KMERS))}
INT_TO_KMER = {i: ALL_KMERS[i] for i in range(len(ALL_KMERS))}
M5C_KMERS = CENTER_MOTIFS

