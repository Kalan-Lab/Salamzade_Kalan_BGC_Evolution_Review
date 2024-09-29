import os
import sys
from scipy import stats
import statistics
import numpy

gene_lens = []
with open('Boxplot_Input.txt') as obif:
    for i, line in enumerate(obif):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gene_lens.append(int(ls[3]))

per25 = numpy.percentile(gene_lens, 25)
per75 = numpy.percentile(gene_lens, 75)

bgc = []
other = []
for i, line in enumerate(open('Boxplot_Input.txt')):
    line = line.strip()
    ls = line.split('\t')
    if i == 0:
        print('\t'.join(ls + ['length_cat']))
        continue
    #if int(ls[3]) > per25 and int(ls[3]) < per75:
    #    print(line + '\tmiddle_50_percentile') 
    if ls[2] == 'True':
        bgc.append(float(ls[1]))
    else:
        other.append(float(ls[1]))
    #print(line + '\tall')
stat, pval = stats.ranksums(bgc, other, alternative='less')
print(pval)
print(stat)
