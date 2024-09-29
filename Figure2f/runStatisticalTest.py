import os
import sys
from scipy import stats

core = []
noncore = []
for i, line in enumerate(open('TajD_Stats.txt')):
    if i == 0: continue
    line = line.strip()
    ls = line.split('\t')
    if ls[2] == 'True':
        core.append(float(ls[3]))
    else:
        noncore.append(float(ls[3]))

stat, pval = stats.mannwhitneyu(core, noncore, alternative='greater')
print(pval)
print(stat)
