import os
import sys

"""
0       genus
1       gcf id
2       gcf annotation
3       homolog group
4       annotation
5       hg order index
6       hg consensus direction
7       hg proportion multi-copy genome-wide
8       median gene length
9       is core to bgc
10      num of hg instances
11      samples with hg
12      proportion of samples with hg
13      ambiguous sites proporition
14      Tajimas D
15      proportion variable sites
16      proportion nondominant major allele
17      median beta rd
18      median dn ds
19      mad dn/ds
20      populations with hg
21      proportion of total populations with hg
22      most significant Fisher exact pvalues presence absence
23      median Tajimas D per population
24      mad Tajimas D per population
25      most negative population Tajimas_D
26      most positive population Tajimas D
27      population entropy
28      median fst-like estimate
29      population proportion of members with hg
30      all domains
"""

print('genus\tgcf\thg\tis_core\ttajimas_d')
with open('S10_Data.txt') as osd:
    for i, line in enumerate(osd):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        genus = ls[0]
        gcf = ls[1]
        hg = ls[3]
        prop = float(ls[12])
        samps = int(ls[11])
        core = ls[9]
        tajd = ls[14]
        sc = float(ls[7])
        if tajd == 'NA' or sc > 0.0 or samps < 10 or prop < 0.75: continue
        print('\t'.join([genus, gcf, hg, core, tajd]))
