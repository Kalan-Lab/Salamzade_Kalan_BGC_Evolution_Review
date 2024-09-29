import os
import sys
from collections import defaultdict
import statistics
from ete3 import Tree
from Bio import SeqIO
import traceback

def determineBranchSumForGroup(species_tree_file, group_genomes):
    try:
        t = Tree(species_tree_file) 
        t.prune(group_genomes, preserve_branch_length=True)
        phylo_breadth = 0.0
        for n in t.traverse('postorder'):
            phylo_breadth += n.dist
        return(phylo_breadth)
    except:
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def determineOgCount(genome_ogs, group_genomes, core_genome=80.0):
    try:
        og_counts = defaultdict(int)
        for g in group_genomes:
            for og in genome_ogs[g]:
                og_counts[og] += 1
        tot_genomes = len(group_genomes)

        core_ogs = set([])
        for og in og_counts:
            og_prop = og_counts[og]/float(tot_genomes)
            if (og_prop*100.0) >= core_genome:
                core_ogs.add(og)

        aux_ogs = set([])
        for g in group_genomes:
            aux_ogs = aux_ogs.union(genome_ogs[g].difference(core_ogs))
		
        tot_aux_ogs = len(aux_ogs)

        return(tot_aux_ogs)
		
    except:
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

cs_lists = 'Clades/'
as_dir = 'AntiSMASH_Results/'
og_file = 'Full_OrthoGroup_Matrix.txt'
species_tree = 'Species_Phylogeny.tre'

genome_to_clade = {}
clade_genomes = defaultdict(list)
for f in os.listdir(cs_lists):
    clade = f.split('_')[0]
    with open(cs_lists + f) as ocf:
        for line in ocf:
            line = line.strip()
            ls = line.split('\t')
            genome_to_clade[ls[0]] = clade
            clade_genomes[clade].append(ls[0])

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

genome_ogs = defaultdict(set)
genomes = []
with open(og_file) as oof:
    for i, line in enumerate(oof):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            genomes = ls[1:]
        else:
            og = ls[0]
            for j, op in enumerate(ls[1:]):
                if op == '1':
                    s = genomes[j]
                    if not s in genome_to_clade: continue
                    genome_ogs[s].add(og)

print('\t'.join(['clade', 'avg_bgcome_prop', 'avg_bgcome_size', 'aux_ogs', 'phylo_breadth']))
for clade in clade_genomes:
    cgs = clade_genomes[clade]
    clade_bs = determineBranchSumForGroup(species_tree, cgs)
    aux_ogs = determineOgCount(genome_ogs, cgs, core_genome=80.0)
    
    avg_bgcome_prop = statistics.mean(bgcome_props[clade])
    avg_bgcome_size = statistics.mean(bgcome_sizes[clade])
    print('\t'.join([str(x) for x in [clade, avg_bgcome_prop, avg_bgcome_size, aux_ogs, clade_bs]]))
