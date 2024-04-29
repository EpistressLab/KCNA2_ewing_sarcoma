########## packages ########

library(ggplot2)
library(tidyverse)

########## GSEA plot ########

GSEA_siCT <- read.table("gsea_report_for_A673_siCT_1710168507845.tsv",sep="\t",h=T)
GSEA_siCT <- separate(GSEA_siCT,"LEADING.EDGE",into=c('a','b','c','d'),sep="=") %>%
             separate("b",into=c("Tags","b2"),sep="%") %>%
             separate("c",into=c("List","c2"),sep="%") %>%
             separate("d",into=c("Signal","d2"),sep="%")
GSEA_siCT <- GSEA_siCT[,-c(11,13,15,17,18)]
GSEA_siCT$Signal <- as.numeric(GSEA_siCT$Signal)*0.01


pdf("dotplot_GSEA_siCT_enriched.pdf",width=11.81, height=9.84)
ggplot(GSEA_siCT[1:30,],aes(x=NES, y=reorder(NAME, -NES))) +
    geom_segment(aes(xend=0, yend = NAME)) +
    geom_label(aes(label=NAME,x=0,y=NAME),hjust=0)+
    geom_point(aes(color=FDR.q.val, size = Signal),alpha=0.5) +
    scale_color_gradient(low = "red",high = "blue",guide = guide_colorbar(reverse=TRUE))+
    scale_size_continuous(range=c(2, 10)) +
    theme_bw() +
    xlab("Normalized Enrichment Score") +
    ylab('Gene set Description') +
    ggtitle("KEGG")+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 10)) +
    scale_x_reverse()
dev.off()