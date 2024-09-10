library(ggplot2)

dat <- read.table('Plot_Input.txt', header=T, sep='\t')

pdf('Plot.pdf', height=3, width=3)
#ggplot(dat, aes(x=is_core, y=tajimas_d)) + geom_boxplot() + theme_classic()
ggplot(dat, aes(x=tajimas_d, fill=is_core)) + geom_histogram(color='black', show.legend=F) + theme_classic() + facet_wrap(~is_core, ncol=1) + scale_fill_manual(values=c('#303030', '#0f0f0f')) + xlab("") + ylab("") + scale_y_log10()
dev.off()

