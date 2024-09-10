## Overview 

- For each of the three species, representative genomes were gathered from previous dereplication efforts for entire genera/families via skDER - that were uploaded to Zenodo (https://zenodo.org/records/10041203). skDER had been run in greedy dereplication mode using cutoffs of 99% ANI and 90% AF.
- antiSMASH (v7.0.0) was run on the selected genomes and, subsequently, lsaBGC-Pan (v1.1.0) was run with Panaroo orthology inference.
- The custom script computePangenomeFrequency.py was used to create plotting input for U-shaped gene frequency plots. Plotting was performed via ggplot2 in R (see script plot.R).
