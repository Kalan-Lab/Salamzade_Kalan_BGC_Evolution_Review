library(ggplot2)

dat <- read.table('pangenome_uplot.txt', header=T, sep='\t')

pdf('Frequency_Plot.pdf', height=3, width=4)
ggplot(dat, aes(x=(10*freq_bin), y=value, fill=bgc_related)) + geom_bar(stat='identity', color='black', bins=10, show.legend=F) + theme_classic() + scale_fill_manual(values=c('#a8a5a5', '#303030')) + xlab("") + ylab("")  +  theme(text = element_text(size = 15))
dev.off()
