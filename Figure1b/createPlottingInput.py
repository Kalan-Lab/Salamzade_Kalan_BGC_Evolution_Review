import os
import sys
from collections import defaultdict
import statistics
from ete3 import Tree
from Bio import SeqIO

cs_lists = 'Clades/'
cs_dir = 'Clade_Subtrees/'
as_dir = 'AntiSMASH_Results/'
og_file = 'Full_OrthoGroup_Matrix.txt'

genome_to_clade = {}
for f in os.listdir(cs_lists):
    clade = f.split('_')[0]
    with open(cs_lists + f) as ocf:
        for line in ocf:
            line = line.strip()
            ls = line.split('\t')
            genome_to_clade[ls[0]] = clade
    
bgcome_props = defaultdict(list)
bgcome_sizes = defaultdict(list)
for s in os.listdir(as_dir): 
    samp_dir = as_dir + s + '/'
    genome_size = 0
    bgcome_size = 0
    for f in os.listdir(samp_dir):
        if f.endswith('.gbk') and not '.region' in f:
            with open(samp_dir + f) as ogf:
                for rec in SeqIO.parse(ogf, 'genbank'):
                    genome_size += len(str(rec.seq))
        elif f.endswith('.gbk') and '.region' in f:
            with open(samp_dir + f) as obf:
                for rec in SeqIO.parse(obf, 'genbank'):
                    bgcome_size += len(str(rec.seq))
    assert(genome_size != 0)
    if not s in genome_to_clade: continue
    bgcome_prop = bgcome_size/float(genome_size)
    clade = genome_to_clade[s]
    bgcome_props[clade].append(bgcome_prop)
    bgcome_sizes[clade].append(bgcome_size)

clade_ogs = defaultdict(set)
samples = []
with open(og_file) as oof:
    for i, line in enumerate(oof):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            samples = ls[1:]
        else:
            og = ls[0]
            for j, op in enumerate(ls[1:]):
                if op == '1':
                    s = samples[j]
                    if not s in genome_to_clade: continue
                    c = genome_to_clade[s]
                    clade_ogs[c].add(og)

print('\t'.join(['clade', 'avg_bgcome_prop', 'avg_bgcome_size', 'tot_ogs', 'phylo_breadth']))
for f in os.listdir(cs_dir):
    clade = f.split('_')[0]
    st = Tree(cs_dir + f)
    phylo_breadth = 0.0
    for n in st.traverse('postorder'):
        phylo_breadth += n.dist
    
    avg_bgcome_prop = statistics.mean(bgcome_props[clade])
    avg_bgcome_size = statistics.mean(bgcome_sizes[clade])
    tot_ogs = len(clade_ogs[clade])
    print('\t'.join([str(x) for x in [clade, avg_bgcome_prop, avg_bgcome_size, tot_ogs, phylo_breadth]]))
