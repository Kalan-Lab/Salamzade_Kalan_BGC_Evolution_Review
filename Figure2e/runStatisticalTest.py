import os
import sys
from scipy import stats

core = []
noncore = []
for i, line in enumerate(open('Plot_Input.txt')):
    if i == 0: continue
    line = line.strip()
    ls = line.split('\t')
    if ls[-2] == 'TRUE':
        core.append(float(ls[-1]))
    else:
        noncore.append(float(ls[-1]))

stat, pval = stats.mannwhitneyu(core, noncore, alternative='less')
print(pval)
print(stat)
