library(ggplot2)

dat <- read.table('Plotting_Input.txt', sep='\t', header=T)
# clade   avg_bgcome_prop tot_ogs phylo_breadth

png('Fig1B_Base.png', height=5, width=10, units='in', res=1200)
ggplot(dat, aes(y=tot_ogs, x=phylo_breadth, color=avg_bgcome_size_rank, size=avg_bgcome_size_rank, label=clade)) + theme_bw() + geom_smooth(method='lm', color='black', linetype=2, fill='grey')  + geom_point( show.legend=F)  + scale_color_gradient(low='#a9aaab',  high='#121212') + xlab("Core genome phylogenetic breadth") + ylab("Total orthogroups") + theme(text = element_text(size = 20), axis.text.x = element_text())
dev.off()
