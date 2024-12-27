import os
import sys
from collections import defaultdict
from Bio import SeqIO

antismash_dir = 'AntiSMASH_Results/'
st5_file = 'Supplementary_Table_S5.txt' # Not included here can be downloaded as Rosconi et al. 2022 supplementary table S5 
st7_file = 'Supplementary_Table_S7.txt' # Not included here can be downloaded as Rosconi et al. 2022 supplementary table S7

bgc_coords = defaultdict(set)
for s in os.listdir(antismash_dir):
    samp_dir = antismash_dir + s + '/'
    for f in os.listdir(samp_dir):
        if '.region' in f and f.endswith('.gbk'):
            start = None
            end = None
            scaffold = None
            protocore_coords = set([])
            with open(samp_dir + f) as osf:
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

bgc_ogs = set([])
with open(st5_file) as osf:
    for i, line in enumerate(osf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clust_id = ls[1]
        scaffold = ls[3]
        start = int(ls[6])
        end = int(ls[7])
        for pos in range(start, end+1):
            if pos in bgc_coords[scaffold]:
                bgc_ogs.add(clust_id)

with open(st7_file) as osf:
    for i, line in enumerate(osf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clust_id = ls[1]
        eclass = ls[7]
        in_bgc = 'False'
        if clust_id in bgc_ogs:
            in_bgc = 'True'
        print(clust_id + '\t' + in_bgc + '\t' + eclass)
