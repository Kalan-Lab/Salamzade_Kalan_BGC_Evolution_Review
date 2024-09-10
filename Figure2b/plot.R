library(ggplot2)

dat <- read.table("Plot_Input.txt", header=T, sep='\t')

png('fitness.png', height=3, width=4, res=600, units='in')
ggplot(dat, aes(x=foc_val, y=com_val, fill=foc_id, color=foc_id)) + geom_point(size=3, show.legend=F, alpha=0.7) + theme_bw() + scale_color_manual(values=c('#ba233a', '#d49442'))+xlab("") + ylab("") + geom_smooth(method='lm', show.legend=F) + scale_fill_manual(values=c('#ba233a', '#d49442'))
dev.off()
