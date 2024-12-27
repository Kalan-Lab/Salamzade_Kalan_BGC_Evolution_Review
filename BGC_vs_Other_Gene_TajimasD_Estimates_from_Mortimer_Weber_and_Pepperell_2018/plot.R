library(ggplot2)


dat <- read.table('Boxplot_Input_with_GeneLengthCategory.txt', header=T, sep='\t')

pdf('Figure_2D.pdf', height=3, width=4)
ggplot(dat, aes(x=length_cat, y=tajimas_d, fill=bgc_related)) + geom_violin(draw_quantiles = c(0.5), show.legend=F) + theme_classic() + scale_fill_manual(values=c('#e0dede', '#8c8c8c')) + xlab("") + ylab("") 
#ggplot(dat, aes(x=tajimas_d, fill=bgc_related)) + geom_histogram(show.legend=F) + theme_classic()  + scale_fill_manual(values=c('#e0dede', '#8c8c8c')) + xlab("") + ylab("") + scale_y_log10() + facet_wrap(~bgc_related, ncol=1)
dev.off()
