import os
import sys
import scipy 
from collections import defaultdict
import statistics

foc_ids = ['SAOUHSC_02877', 'SAOUHSC_02879']
fit_file = 'full_fitvals.txt'

print('foc_id\tfoc_val\tcom_val')
for foc_id in foc_ids:
    foc_vals = []
    with open(fit_file) as odf:
        for i, line in enumerate(odf):
            if i == 0: continue
            line = line.strip("\n")
            ls = line.split('\t')
            if ls[0] == foc_id:
                foc_vals = ls[1:]
            
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
            if ls[0] != 'SAOUHSC_02962': continue
            vals = ls[1:]
            for j, val in enumerate(vals):
                print(foc_id + '\t' + foc_vals[j] + '\t' + val)
