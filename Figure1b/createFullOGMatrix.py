import os
import sys
from collections import defaultdict

og_mat = 'OrthoFinder_Results/Results_Aug12/Orthogroups/Orthogroups.tsv'
sog_mat = 'OrthoFinder_Results/Results_Aug12/Orthogroups/Orthogroups_UnassignedGenes.tsv'

with open(og_mat) as oom:
    for i, line in enumerate(oom):
        line = line.strip('\n')
        ls = line.split('\t') 
        if i == 0:
            print('\t'.join(['.fna'.join(x.split('.fna')[:-1]) for x in ls]))
        else:
            og = ls[0]
            printlist = [og]
            for lts in ls[1:]:
                if lts.strip() != '':
                    printlist.append('1')
                else:
                    printlist.append('0')
            print('\t'.join(printlist))

with open(sog_mat) as som:
    for i, line in enumerate(som):
        line = line.strip('\n')
        ls = line.split('\t')
        if i != 0:
            og = ls[0]
            printlist = [og]
            for lts in ls[1:]:
                if lts.strip() != '':
                    printlist.append('1')
                else:
                    printlist.append('0')
            print('\t'.join(printlist))

