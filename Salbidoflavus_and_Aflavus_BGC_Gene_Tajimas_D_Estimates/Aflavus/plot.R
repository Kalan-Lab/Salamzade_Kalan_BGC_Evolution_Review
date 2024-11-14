library(ggplot2)

dat <- read.table('TajD_Stats.txt', header=T, sep='\t')

pdf('Aflavus.pdf', height=3, width=5)
#ggplot(dat, aes(x=is_core, y=tajimas_d)) + geom_boxplot() + theme_classic()
ggplot(dat, aes(x=tajimas_d, fill=protocore)) + geom_histogram(color='black', show.legend=F) + theme_classic() + facet_wrap(~protocore, ncol=1) + scale_fill_manual(values=c('#303030', '#0f0f0f')) + xlab("") + ylab("") + scale_y_log10()
dev.off()

