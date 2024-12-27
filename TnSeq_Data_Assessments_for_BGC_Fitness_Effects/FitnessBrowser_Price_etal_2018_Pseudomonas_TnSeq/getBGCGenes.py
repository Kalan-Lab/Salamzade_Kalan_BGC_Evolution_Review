import os
import sys
from collections import defaultdict
from Bio import SeqIO

antismash_dir = 'AntiSMASH_Results/'
genes_file = 'organism_pseudo3_N2E3_genes.tab'

bgc_coords = defaultdict(set)
for f in os.listdir(antismash_dir):
    if '.region' in f and f.endswith('.gbk'):
        start = None
        end = None
        scaffold = None
        protocore_coords = set([])
        with open(antismash_dir + f) as osf:
            for line in osf:
                line = line.strip()
                if 'ACCESSION' in line:
                    scaffold = line.split()[1]
                elif 'Orig. start  ::' in line:
                    start = int(line.split('::')[1])
                elif 'Orig. end    ::' in line:
                    end = int(line.split('::')[1])
                elif 'proto_core' in line:
                    s, e = [int(x) for x in line.split()[-1].split('..')]
                    for p in range(s, e+1):
                        protocore_coords.add(p)
            
        for i, pos in enumerate(range(start, end+1)):
            if (i+1) in protocore_coords:
                bgc_coords[scaffold].add(pos)

with open(genes_file) as ogf:
    for i, line in enumerate(ogf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        scaffold = ls[2]
        gid = ls[0]
        bgc_protocore = False
        for pos in range(int(ls[3]), int(ls[4])+1):
            if pos in bgc_coords[scaffold]:
                bgc_protocore = True
        if bgc_protocore:
            print(gid)
