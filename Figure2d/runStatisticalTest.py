import os
import sys
from scipy import stats

bgc = []
other = []
for i, line in enumerate(open('Boxplot_Input.txt')):
    if i == 0: continue
    line = line.strip()
    ls = line.split('\t')
    if ls[-1] == 'True':
        bgc.append(float(ls[-2]))
    else:
        other.append(float(ls[-2]))

stat, pval = stats.ranksums(bgc, other, alternative='less')
print(pval)
print(stat)
