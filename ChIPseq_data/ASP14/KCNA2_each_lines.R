
library(ggplot2)
library(gggenes)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)

##################### Functions ###############

make_cov_plot <- function(data_cov,start,end,max_y){
  plot <- ggplot(data_cov, aes(x=position)) +
  theme_minimal() +
  theme(legend.position="none") %+replace% theme(
    axis.ticks.y = element_blank(),
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9)) +
  xlim(start,end) + 
  ylim(0,max_y) +
  xlab("Position on chr1") +
  geom_line(aes(y = reads_all), linewidth=0.5) +
  geom_area(aes(y = reads_all), linewidth=0.05,linetype = 0, fill="#68bd8e", alpha=0.5, outline.type="full")
  return(plot)
}

make_macs2_plot <- function(macs2_data,start,end,min_gradient,max_gradient){
    plot_macs2 <- ggplot(macs2_data) +
        geom_rect(aes(xmin=V2,xmax=V3,ymin=0,ymax=1,fill=V9),alpha=0.7) +
        xlim(start,end) +
        theme_minimal() +
        theme(legend.position="none") %+replace% theme(
            axis.ticks.y = element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        scale_fill_gradient(
            low = "#68bd8e",
            high = "#00451f",
            limits=c(min_gradient,max_gradient))
    return(plot_macs2)
}

make_macs2_legend <- function(macs2_data,start,end,min_gradient,max_gradient){
    #macs2_data <- read.table(macs2_file, h=F, sep="\t")
    plot_macs2 <- as_ggplot(get_legend(ggplot(macs2_data) +
        geom_rect(aes(xmin=V2,xmax=V3,ymin=0,ymax=1,fill=V9),alpha=0.7) +
        xlim(start,end) +
        theme_minimal() + 
        theme(legend.title = element_text(hjust = 0.5, size=9)) +
        scale_fill_gradient(
            low = "#68bd8e",
            high = "#00451f",
            limits=c(min_gradient,max_gradient)) +
        labs(fill=bquote(atop(bold("Macs2"),"\n-log10(qvalue)")))))
    return(plot_macs2)
}

###################### H3K4me3 KCNA2 plots ############

####### GTF file
gtf <- read.table("~/bird/bone_epigenetics/anais/genomes/hg38/GCF_000001405.40_GRCh38.p14_genomic.gtf", sep="\t", h=F)
gtf <- gtf %>% separate('V9',into=c('a','b','c','d','e','f','g','h','i'),sep="; ",remove=T) %>% separate('a',into=c("gene_id","gene"),sep=" ") %>% separate('b',into=c("transcript_id","transcript"),sep=" ")
gtf <- gtf[,-c(2,9,11,13:19)]

gtf_KCNA2 <- subset(gtf,gene=="KCNA2")
gtf_KCNA2_transcripts <- subset(gtf_KCNA2, V3=='transcript')
gtf_KCNA2_transcripts <- gtf_KCNA2_transcripts[,-c(2,5,7)]
colnames(gtf_KCNA2_transcripts) <- c("chr","start_transcript","end_transcript","strand","gene","transcript")


gtf_KCNA2 <- subset(gtf_KCNA2, V3 %in% c("exon","CDS"))
gtf_KCNA2 <- merge(gtf_KCNA2,gtf_KCNA2_transcripts)
gtf_KCNA2$orientation <- 1
gtf_KCNA2$orientation[which(gtf_KCNA2$V7 == '-')] <- 0
gtf_KCNA2$transcript <- factor(gtf_KCNA2$transcript, levels = c("NM_004974.4","XM_011541396.3","XM_017001213.2","XM_011541400.3","NM_001204269.2","XM_011541398.3"))
transcripts_nicknames <- data.frame('transcript'=c("NM_004974.4","XM_011541396.3","XM_017001213.2","XM_011541400.3","NM_001204269.2","XM_011541398.3"),'nickname'=c('F','E','D','C','B','A'))
gtf_KCNA2 <- merge(gtf_KCNA2,transcripts_nicknames)
gtf_KCNA2$nickname <- factor(gtf_KCNA2$nickname, levels=c('F','E','D','C','B','A'))
gtf_KCNA2$nickname_2 <- paste(gtf_KCNA2$nickname,gtf_KCNA2$transcript,sep=" : ")
gtf_KCNA2$nickname_2 <- factor(gtf_KCNA2$nickname_2, levels=c("F : NM_004974.4","E : XM_011541396.3","D : XM_017001213.2","C : XM_011541400.3","B : NM_001204269.2","A : XM_011541398.3"))

plot_annot <- ggplot(gtf_KCNA2,aes(xmin=start_transcript, xmax=end_transcript, y=nickname, forward = orientation)) +
    geom_gene_arrow(fill="white") +
    geom_subgene_arrow(data=subset(gtf_KCNA2, V3=="exon"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=0.6, fill='coral1') +
    geom_subgene_arrow(data=subset(gtf_KCNA2, V3=="CDS"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=1, fill='coral4') +
    theme_genes() +
    theme(legend.position="none") %+replace% theme(
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
    xlim(110583000,110641000) +
    ylab(expression(italic(KCNA2)~" transcripts"))

######## Coverage plots

#GGAA_coords <- data.frame('start'=c(110611020),'end'=c(110611071))
GGAA_coords <- data.frame('start'=c(110610970),'end'=c(110611121))

# H3K27ac
ASP14_d0_H3K27ac_cov <- read.table("./Samples/ASP14_d0_H3K27ac/ASP14_d0_H3K27ac_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d0_H3K27ac_macs2 <- subset(read.table("./Samples/ASP14_d0_H3K27ac/macs2/ASP14_d0_H3K27ac_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d7_H3K27ac_cov <- read.table("./Samples/ASP14_d7_H3K27ac/ASP14_d7_H3K27ac_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d7_H3K27ac_macs2 <- subset(read.table("./Samples/ASP14_d7_H3K27ac/macs2/ASP14_d7_H3K27ac_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")

max_cov_plot <- max(c(ASP14_d0_H3K27ac_cov$reads_all,ASP14_d7_H3K27ac_cov$reads_all))
plot_ASP14_d0_H3K27ac <- make_cov_plot(ASP14_d0_H3K27ac_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 0"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)
plot_ASP14_d7_H3K27ac <- make_cov_plot(ASP14_d7_H3K27ac_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 7"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)

all_peaks_H3K27ac <- rbind(ASP14_d0_H3K27ac_macs2,ASP14_d7_H3K27ac_macs2)
all_peaks_H3K27ac <- subset(all_peaks_H3K27ac, V2>=110583000 & V3 <= 110641000)
min_gradient <- floor(min(all_peaks_H3K27ac$V9))
max_gradient <- ceiling(max(all_peaks_H3K27ac$V9))

legend_macs2_H3K27ac <- make_macs2_legend(ASP14_d0_H3K27ac_macs2,110583000,110641000,min_gradient,max_gradient)
plot_macs2_ASP14_d0_H3K27ac <- make_macs2_plot(ASP14_d0_H3K27ac_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 0"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)
plot_macs2_ASP14_d7_H3K27ac <- make_macs2_plot(ASP14_d7_H3K27ac_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 7"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)

png("plots/KCNA2_ASP14_H3K27ac.png", width = 20, height = 21, units = "cm", res = 300)
plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label="ASP14 H3K27ac ChIP sequencing") + theme_void(),
            ggplot() + xlim(110583000,110641000) + annotate("text",x=max(GGAA_coords$start), y = 1, size=3, label="(GGAA)[n]",parse=TRUE) + theme_void(),
            plot_ASP14_d0_H3K27ac,plot_macs2_ASP14_d0_H3K27ac,
            plot_ASP14_d7_H3K27ac,plot_macs2_ASP14_d7_H3K27ac,
            plot_annot, align = "v", axis="tb", ncol=1, nrow=7, rel_heights = c(0.04,0.04,rep(c(0.25, 0.06),2),0.3)),
        ggdraw(),legend_macs2_H3K27ac, ncol=3, rel_widths=c(0.85,0.05,0.15))
dev.off()

# FLI1
ASP14_d7_FLI1_cov <- read.table("./Samples/ASP14_d7_FLI1/ASP14_d7_FLI1_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d7_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d7_FLI1/macs2/ASP14_d7_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d10_FLI1_cov <- read.table("./Samples/ASP14_d10_FLI1/ASP14_d10_FLI1_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d10_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d10_FLI1/macs2/ASP14_d10_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d11_FLI1_cov <- read.table("./Samples/ASP14_d11_FLI1/ASP14_d11_FLI1_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d11_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d11_FLI1/macs2/ASP14_d11_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d14_FLI1_cov <- read.table("./Samples/ASP14_d14_FLI1/ASP14_d14_FLI1_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d14_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d14_FLI1/macs2/ASP14_d14_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d17_FLI1_cov <- read.table("./Samples/ASP14_d17_FLI1/ASP14_d17_FLI1_KCNA2_coverage.tsv", h=T, sep="\t")
ASP14_d17_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d17_FLI1/macs2/ASP14_d17_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")

max_cov_plot <- max(c(ASP14_d7_FLI1_cov$reads_all,ASP14_d10_FLI1_cov$reads_all,ASP14_d11_FLI1_cov$reads_all,ASP14_d14_FLI1_cov$reads_all,ASP14_d17_FLI1_cov$reads_all))

plot_ASP14_d7_FLI1 <- make_cov_plot(ASP14_d7_FLI1_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 7"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)
plot_ASP14_d10_FLI1 <- make_cov_plot(ASP14_d10_FLI1_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 7 + 3"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)
plot_ASP14_d11_FLI1 <- make_cov_plot(ASP14_d11_FLI1_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 7 + 4"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)
plot_ASP14_d14_FLI1 <- make_cov_plot(ASP14_d14_FLI1_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 7 + 7"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)
plot_ASP14_d17_FLI1 <- make_cov_plot(ASP14_d17_FLI1_cov,110583000,110641000,max_cov_plot) + 
    ylab(bquote(atop(bold("Day 7 + 10"),"Depth coverage"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="red", alpha=1)


all_peaks_FLI1 <- rbind(rbind(rbind(rbind(ASP14_d7_FLI1_macs2,ASP14_d10_FLI1_macs2),ASP14_d11_FLI1_macs2),ASP14_d14_FLI1_macs2),ASP14_d17_FLI1_macs2)
all_peaks_FLI1 <- subset(all_peaks_FLI1, V2>=110583000 & V3 <= 110641000)
min_gradient <- floor(min(all_peaks_FLI1$V9))
max_gradient <- ceiling(max(all_peaks_FLI1$V9))

legend_macs2_FLI1 <- make_macs2_legend(ASP14_d7_FLI1_macs2,110583000,110641000,min_gradient,max_gradient)
plot_macs2_ASP14_d7_FLI1 <- make_macs2_plot(ASP14_d7_FLI1_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 7"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)
plot_macs2_ASP14_d10_FLI1 <- make_macs2_plot(ASP14_d10_FLI1_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 7 + 3"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)
plot_macs2_ASP14_d11_FLI1 <- make_macs2_plot(ASP14_d11_FLI1_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 7 + 4"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)
plot_macs2_ASP14_d14_FLI1 <- make_macs2_plot(ASP14_d14_FLI1_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 7 + 7"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)
plot_macs2_ASP14_d17_FLI1 <- make_macs2_plot(ASP14_d17_FLI1_macs2,110583000,110641000,min_gradient,max_gradient) + 
    ylab(bquote(atop(bold("Day 7 + 10"),"Macs2"))) +
    geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="red", alpha=1)

png("plots/KCNA2_ASP14_FLI1.png", width = 20, height = 28, units = "cm", res = 300)
plot_grid(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label="ASP14 FLI1 ChIP sequencing") + theme_void(),
            ggplot() + xlim(110583000,110641000) + annotate("text",x=max(GGAA_coords$start), y = 1, size=3, label="(GGAA)[n]",parse=TRUE) + theme_void(),
            plot_ASP14_d7_FLI1,plot_macs2_ASP14_d7_FLI1,
            plot_ASP14_d10_FLI1,plot_macs2_ASP14_d10_FLI1,
            plot_ASP14_d11_FLI1,plot_macs2_ASP14_d11_FLI1,
            plot_ASP14_d14_FLI1,plot_macs2_ASP14_d14_FLI1,
            plot_ASP14_d17_FLI1,plot_macs2_ASP14_d17_FLI1,
            plot_annot, align = "v", axis="tb", ncol=1, nrow=13, rel_heights = c(0.04,0.03,rep(c(0.20, 0.09),5),0.28)),
        ggdraw(),legend_macs2_FLI1, ncol=3, rel_widths=c(0.85,0.05,0.15)))
dev.off()