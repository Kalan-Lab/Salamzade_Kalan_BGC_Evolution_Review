## Overview

- population genetic stats for *Mycobacterium tuberculosis* were downloaded from GitHub repos associated with the study Mortimer et al. 2018. Specifically the genewiseStats.txt file was gathered from: https://github.com/pepperell-lab/mtbDrugResistance/blob/master/data/genewiseStats.txt
- BGC region coordinates were determined by running antiSMASH (v7.0.0) on the refernece genome Mortimer et al. 2018 had used.
- Custom scripts (`createBoxPlotInput.py` and `runStatisticalTest.py`) were used to bin population genetic stats depending on whether genes are within or outside BGC regions, to determine whether gene lengths were in the 25th to 75th percentile, and perform statistical testing.
