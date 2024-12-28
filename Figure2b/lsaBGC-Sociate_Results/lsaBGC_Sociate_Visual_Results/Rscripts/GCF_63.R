library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(dplyr)
library(grid)
library(gridExtra)

phylo.tree <- read.tree("/workspace/lab/kalanlab/salamzar/lsaBGC/minireview/Strep_albido/lsaBGC-Pan_Results_v1.0.9/Create_Species_Phylogeny/Species_Phylogeny_FastTree2.treefile")
phylo.tree <- midpoint.root(phylo.tree)
tree.labels <- phylo.tree$tip.label

heatmap.data <- read.table("/workspace/lab/kalanlab/salamzar/lsaBGC/minireview/Strep_albido/lsaBGC-Pan_Results_v1.0.9/lsaBGC_Sociate_Results/Tracks_for_Visual/GCF_63.txt", header=F, sep="\t")
colnames(heatmap.data) <- c("allele", "label", "order", "value")
colors <- c("#242424", "#6f7070", "#4b7dc9", "#db6588")
names(colors) <- c("focal", "-", "positive", "negative")

pdf("/workspace/lab/kalanlab/salamzar/lsaBGC/minireview/Strep_albido/lsaBGC-Pan_Results_v1.0.9/Final_Results/Visualizations/lsaBGC_Sociate_Visual_Results/Plots/GCF_63.pdf", height=10, width=20)
gg_tr <- ggtree(phylo.tree)
gg_hm <- ggplot(heatmap.data, aes(x = reorder(allele, order), y = label, fill = value)) +
theme_classic() + scale_fill_manual(values=colors, na.value="white") +
xlab("") + ylab("") + geom_tile(color="white", show.legend=F) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg_hm %>% insert_left(gg_tr, width=0.4)
dev.off()

