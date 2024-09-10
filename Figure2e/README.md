## Overview

- Population genetics data was gathered from Supplementary Data table S10 from Salamzade et al. 2023, Microbial Genomics for BGC information
- The script filterDataForPlotting.py was used to filter for specific genes that met stringent criteria for inclusion in the analysis. Namely, genes included must be present in more than 10 samples, be single-copy across genomes, and must be present in at least 75% of BGC instances for the GCF they are associated with.
- Plotting was performed via ggplot2 in R.
