## Overview

- TnSeq fitness values were gathered as a matrix from the spreadsheet "Supplementary Data" in Santiago et al. 2018 (https://www.nature.com/articles/s41589-018-0041-4) - is referred to in scripts as the file: "full_fitvals.txt".
- Sequence alignment was used to determine locus tag identifiers of crt genes in the genome of the strain used for TnSeq analysis.
- Custom scripts were written (included here; performCorrelationsAndTakeAverageAcrossFocalGenes.py and createPlotInput.py) to first determine which genes are correlated in fitness response with *BOTH* crtM and crtN and to then create the input for plotting Figure 2B. Plotting itself was performed using ggplot2 in R.
