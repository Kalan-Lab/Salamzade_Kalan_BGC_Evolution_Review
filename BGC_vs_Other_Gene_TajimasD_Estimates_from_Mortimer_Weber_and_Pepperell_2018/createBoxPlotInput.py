import os
import sys
from Bio import SeqIO

as_dir = 'Mtb_AntiSMASH7/'

bgc_ogs = set([])
for f in os.listdir(as_dir):
    if f.endswith('.gbk') and '.region' in f:
        with open(as_dir + f) as of:
            for rec in SeqIO.parse(of, 'genbank'):
                for feat in rec.features:
                    if not feat.type == 'CDS': continue
                    lt = feat.qualifiers['locus_tag'][0]
                    bgc_ogs.add(lt)

gene_to_len = {}
with open('geneLengths.txt') as oglf: # was downloaded from https://raw.githubusercontent.com/tatumdmortimer/mtbDrugResistance/refs/heads/master/data/geneLengths.txt
    for line in oglf:
        line = line.strip()
        gene, length = line.split('\t')
        gene_to_len[gene] = length

name_to_num = {}
with open('rv_numbers.txt') as orn:
    for line in orn:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split('\t')
        name_to_num[ls[0]] = ls[1]

print('rv\ttajimas_d\tbgc_related\tgene_length')
with open('genewiseStats.txt') as ogs: # was downloaded from https://raw.githubusercontent.com/tatumdmortimer/mtbDrugResistance/refs/heads/master/data/genewiseStats.txt
    for i, line in enumerate(ogs):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        if not ls[0] in name_to_num: continue
        rv = name_to_num[ls[0]]
        if not rv.startswith('Rv') or ls[-1] == 'None': continue
        bgc_related = 'False'
        if rv in bgc_ogs:
            bgc_related = 'True'
        print(rv + '\t' + str(ls[-1]) + '\t' + bgc_related + '\t' + gene_to_len[ls[0]])
