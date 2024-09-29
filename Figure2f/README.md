## Overview

- Input for lsaBGC-Pan (v1.1.1) was run with antiSMASH (v7.0.0) BGC predictions for 10 Aspergillus flavus genomes (see Genome_Listing.txt) using gene calling information available in GFF format from NCBI's GenBank database.
- The script determineWhichGenesAreKey.py was used to process zol (v1.4.12) results from lsaBGC-Pan analysis of Streptomyces albidoflavus antiSMASH (v7.0.0) results. This lsaBGC-Pan analysis was also used for Figure2a (reference its repo for more info). Specifically Tajima's D values were extract for all BGC genes. The script also processes the antiSMASH BGC region GenBank files to determine which genes are key biosynthesis ones and which are auxiliary.
- Only ortholog groups which were single-copy in the context of the focal GCF, were found in >= 80% of GCF instances, and were found in >= 4 GCF instances were considered.
- Comparison between BGC auxiliary and key/protocore genes was performed using the one-sided Wilcoxon rank-sum test using scipy.
- Plotting was performed via ggplot2 in R.
