## Overview of Analysis

- skDER (v1.2.3) was used to select a distinct set of 258 Streptomyces genomes from across the genus. The dereplication cutoffs used were 85% ANI, 25% AF in greedy mode. The input provided was a set of distinct Streptomyces representative genomes previously selected from GTDB R214 using skDER at 99% ANI, 90% AF in greedy dereplication mode and uploaded on Zenodo (https://zenodo.org/records/10041203).
- Pyrodigal was used to perform gene calling and predict coding sequences, which were then provided as input to OrthoFinder (v2.5.5) for orthology inference and GToTree (1.8.6) for phylogeny inference. OrthoFinder was run only until coarse orthology inference (hierarchical orthogroups were not computed). FastTree2 (v2.1.11) was used within GToTree for phylogeny construction.
- AntiSMASH (v7.0.0) was used for BGC prediction with largely default parameters (see antismash.cmds). 
- The web application iTol was used to visualize the phylogeny and coarse clades manually defined.
- Custom script (included here, createPlottingInput.py) was written to determine the number of distinct orthogroups (including singletons), prune the overall species tree for individual clades and aggregate phylogenetic branch length for them (via ete3 v3.1.3), and the average BGC-ome size across genomes for each of the nine clades.
- Visualization of the scatterplot was performed using ggplot2 in R (see script plot.R).
