import os
import sys
import scipy 
from collections import defaultdict
import statistics

# crt gene IDs
foc_ids = ['SAOUHSC_02877', 'SAOUHSC_02879'] #'SAOUHSC_02882', 'SAOUHSC_02881', 'SAOUHSC_02880', 'SAOUHSC_02879']
fit_file = 'full_fitvals.txt' # full file of fitness values from Santiago et al. 2018

corr_stats = defaultdict(list)
for foc_id in foc_ids:
    foc_vals = []
    with open(fit_file) as odf:
        for i, line in enumerate(odf):
            if i == 0: continue
            line = line.strip("\n")
            ls = line.split('\t')
            if ls[0] == foc_id:
                foc_vals = [float(x) for x in ls[1:]]
                  
    with open(fit_file) as odf:
        for i, line in enumerate(odf):
            if i == 0: continue
            line = line.strip('\n')
            ls = line.split('\t')
            missing = False
            for val in ls[1:]:
                if val.strip() == '':
                    missing = True
            if missing: continue
            vals = [float(x) for x in ls[1:]]
            stat, pval = scipy.stats.spearmanr(foc_vals, vals)
            corr_stats[ls[0]].append(stat)

for g in corr_stats:
    mean_cor = statistics.mean(corr_stats[g])
    print(g + '\t' + str(mean_cor))
