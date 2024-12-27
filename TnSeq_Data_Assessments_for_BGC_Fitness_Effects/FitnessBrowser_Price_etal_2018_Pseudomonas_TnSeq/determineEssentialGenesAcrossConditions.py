import os
import sys
from collections import defaultdict

fit_file = 'fit_organism_pseudo3_N2E3.tsv' # not included, can be downloaded from: https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=pseudo3_N2E3
t_file = 't_organism_pseudo3_N2E3.tsv' # not included, can be downloaded from: https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=pseudo3_N2E3
bgc_genes_file = 'BGC_Protocore_Genes.txt'

bgc_genes = set([])
with open(bgc_genes_file) as obgf:
    for line in obgf:
        line = line.strip()
        bgc_genes.add(line)

t_values = defaultdict(dict)
conds = []
with open(t_file) as otf:
    for i, line in enumerate(otf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            conds = ls[5:]
        else:
            gid = ls[1]
            for j, val in enumerate(ls[5:]):
                val = float(val)
                c = conds[j]
                t_values[gid][c] = val

conds = []
with open(fit_file) as off:
    for i, line in enumerate(off):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            conds = ls[5:]
        else:
            gid = ls[1]
            for j, val in enumerate(ls[5:]):
                fit = float(val)
                c = conds[j]
                t = t_values[gid][c]
                if abs(t) <= 4: continue
                if fit <= -1:
                    bgc_gene = 'False'
                    if gid in bgc_genes:
                        bgc_gene = 'True'
                    print(gid + '\t' + c + '\t' + str(fit) + '\t' + str(t) + '\t' + str(bgc_gene))
