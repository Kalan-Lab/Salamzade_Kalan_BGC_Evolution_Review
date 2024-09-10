import os
import sys
from Bio import SeqIO
import math
from collections import defaultdict

lsabgc_pan_resdir = os.path.abspath(sys.argv[1]) + "/"

og_matrix_file = lsabgc_pan_resdir + 'Sample_by_Orthogroups_Matrix.tsv'
reconcile_file = lsabgc_pan_resdir + 'Final_Results/Visualizations/lsaBGC_Reconcile_Results/Orthogroup_Summary_Info.txt'

bgc_ogs = set([])
with open(reconcile_file) as orf:
    for i, line in enumerate(orf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        bgc_ogs.add(ls[0])

#print('og\tfrequency\tfrequency_bin\tbgc_related')
tot_samples = None
freq_bin_counts = defaultdict(lambda: defaultdict(int))
with open(og_matrix_file) as oomf:
    for i, line in enumerate(oomf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            tot_samples = len(ls[1:])
        else:
            og = ls[0]
            samp_with = 0
            for lts in ls[1:]:
                if lts.strip() != '':
                    samp_with += 1
            freq = samp_with/tot_samples
            is_bgc_og = 'False'
            if og in bgc_ogs:
                is_bgc_og = 'True'
            freq_bin = math.floor(freq/0.1)
            #print(og + '\t' + str(freq) + '\t' + str(freq_bin) + '\t' + is_bgc_og)
            freq_bin_counts[freq_bin][is_bgc_og] += 1

print('freq_bin\tbgc_related\tvalue')
for fb in freq_bin_counts:
    for tf in ['True', 'False']:
        print(str(fb) + '\t' + tf + '\t' + str(freq_bin_counts[fb][tf]))
