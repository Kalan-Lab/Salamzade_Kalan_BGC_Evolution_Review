## Overview

- The genome, coding gene coordinates, fitness data (both values and t scores) for strain: *Pseudomonas fluorescens* FW300-N2E3 were downloaded from LBL's fitness browser (accessed on 12/23/2024): https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=pseudo3_N2E3 . This browser is associated with the publiation: Price et al. 2018: https://doi.org/10.1038/s41586-018-0124-0
- AntiSMASH (v7.0.0) was used for BGC predictions with largely default parameters (see `antismash.cmd`).
- The custom script `getBGCGenes.py` was used to identify gene IDs from `organism_pseudo3_N2E3_genes.tab` which overlapped with protocore regions for BGC predictions. This resulted in the identification of 26 such genes listed in "BGC_Protocore_Genes.txt".
- Finally, the script `determineEssentialGenesAcrossConditions.py` was used to determine essential or significantly growth disadvantageous genes across different conditions using fitness values and t scores for genes. Based on information provided on the "Help" page of the fitness browser, cutoffs of abs(t) > 4 and fitness value < -1 were employed. Results of the script are stored in the file: `Essential_Genes_across_Conditions.txt`.
