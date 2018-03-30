


library(ggplot2)


correlate_gs_ngs_vaf <- function(file="vaf_concordance_overlap.tsv"){
    df <- read.table(file, header=TRUE)
    print(cor.test(df$genescan, df$ngs))
    q <- ggplot(data=df, aes(x=genescan, y=ngs)) + geom_point(size=5, alpha=0.5) + geom_smooth(method='lm',formula=y~x, se=FALSE) + theme_classic() + xlab("GeneScan VAF estimate (%)") + ylab("NGS VAF estimate (%)") + ggtitle("GeneScan - NGS VAF estimate", subtitle="ITDs of MOLM14 & dx and rl evolution samples detected by both methods") 
    pdf("mrd_genescan_vs_ngs_vaf_dotplot.pdf")
    print(q)
    dev.off()
    
}





correlate_gs_ngs_vaf()
