
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
  geom_area(aes(y = reads_all), linewidth=0.05,linetype = 0, fill="#0000c0", alpha=0.5, outline.type="full")
  return(plot)
}

make_macs2_plot <- function(macs2_data,start,end,min_gradient,max_gradient){
    #macs2_data <- read.table(macs2_file, h=F, sep="\t")
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
            low = "#8BC0E5",
            high = "#0000c0",
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
            low = "#8BC0E5",
            high = "#0000c0",
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

pdf("KCNA2_transcripts.pdf",width = 5.5,height = 3.15)
ggplot(gtf_KCNA2,aes(xmin=start_transcript, xmax=end_transcript, y=reorder(nickname_2,nickname), forward = orientation)) +
    geom_gene_arrow(fill="white") +
    geom_subgene_arrow(data=subset(gtf_KCNA2, V3=="exon"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=0.6, fill='coral1') +
    geom_subgene_arrow(data=subset(gtf_KCNA2, V3=="CDS"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=1, fill='coral4') +
    theme_genes() +
    theme(legend.position="none") %+replace% theme(
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
    xlim(110583000,110641000) +
    xlab(expression(italic(KCNA2)~" transcripts"))
dev.off()

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

cells <- c('A673','CHLA10','CHLA25','EW1','EW3','EW7','EW22','EW24','MHHES1','MIC','POE','RDES','RH1','SKES1','SKNMC','TC32','TC71','TC106')
#GGAA_coords <- data.frame('start'=c(110611020),'end'=c(110611071))
GGAA_coords <- data.frame('start'=c(110610970),'end'=c(110611121))

for (cell in cells){
    ## Data reading
    H3K4me3_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_H3K4me3_KCNA2_coverage.tsv"),collapse=""), h=T, sep="\t")
    H3K4me3_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_H3K4me3_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr1")
    H3K27ac_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_H3K27ac_KCNA2_coverage.tsv"),collapse=""), h=T, sep="\t")
    H3K27ac_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_H3K27ac_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr1")
    H3K27me3_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_H3K27me3_KCNA2_coverage.tsv"),collapse=""), h=T, sep="\t")
    H3K27me3_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_H3K27me3_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr1")
    transcription_factor_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_transcription_factor_KCNA2_coverage.tsv"),collapse=""), h=T, sep="\t")
    transcription_factor_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_transcription_factor_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr1")

    ## For a plot with all sequencings
    # Coverage plots
    max_cov_plot <- max(c(H3K4me3_cov$reads_all,H3K27ac_cov$reads_all,H3K27me3_cov$reads_all,transcription_factor_cov$reads_all))
    plot_H3K4me3_cov <- make_cov_plot(H3K4me3_cov,110583000,110641000,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K4me3"),"Depth coverage"))) +
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_H3K27ac_cov <- make_cov_plot(H3K27ac_cov,110583000,110641000,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K27ac"),"Depth coverage")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_H3K27me3_cov <- make_cov_plot(H3K27me3_cov,110583000,110641000,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K27me3"),"Depth coverage")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_transcription_factor_cov <- make_cov_plot(transcription_factor_cov,110583000,110641000,max_cov_plot)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("ERG"),"Depth coverage")))}
    else {plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("FLI1"),"Depth coverage")))}

    # Macs2 plots
    all_peaks <- rbind(rbind(rbind(H3K4me3_macs2,transcription_factor_macs2),H3K27ac_macs2),H3K27me3_macs2)
    all_peaks <- subset(all_peaks, V2>=110583000 & V3 <= 110641000)
    min_gradient <- floor(min(all_peaks$V9))
    max_gradient <- ceiling(max(all_peaks$V9))
    legend_macs2 <- make_macs2_legend(H3K4me3_macs2,110583000,110641000,min_gradient,max_gradient)
    plot_macs2_H3K4me3 <- make_macs2_plot(H3K4me3_macs2,110583000,110641000,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K4me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_macs2_H3K27ac <- make_macs2_plot(H3K27ac_macs2,110583000,110641000,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K27ac"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_macs2_H3K27me3 <- make_macs2_plot(H3K27me3_macs2,110583000,110641000,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K27me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_macs2_transcription_factor <- make_macs2_plot(transcription_factor_macs2,110583000,110641000,min_gradient,max_gradient)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("ERG"),"Macs2")))}
    else {plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("FLI1"),"Macs2")))}

    # Merging plots
    pdf(paste(c("plots/KCNA2_",cell,"_all_seq.pdf"),collapse=""),width = 7.87,height = 11.8)
    print(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label=cell) + theme_void(),
                ggplot() + xlim(110583000,110641000) + annotate("text",x=max(GGAA_coords$start), y = 1, size=3, label="(GGAA)[n]",parse=TRUE) + theme_void(),
                plot_H3K4me3_cov,plot_macs2_H3K4me3,
                plot_H3K27ac_cov, plot_macs2_H3K27ac,
                plot_H3K27me3_cov, plot_macs2_H3K27me3,
                plot_transcription_factor_cov,plot_macs2_transcription_factor,
                plot_annot, align = "v", axis="tb", ncol=1, nrow=11, rel_heights = c(0.02,0.02,rep(c(0.12, 0.07),4),0.2)),
            ggdraw(),legend_macs2, ncol=3, rel_widths=c(0.85,0.05,0.15)))
    dev.off()

    ## For a plot with only transcription factor, H3K4me3 and H3K27ac sequencings
    # Coverage plots
    max_cov_plot <- max(c(H3K4me3_cov$reads_all,H3K27ac_cov$reads_all,transcription_factor_cov$reads_all))
    plot_H3K4me3_cov <- make_cov_plot(H3K4me3_cov,110583000,110641000,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K4me3"),"Depth coverage"))) +
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_H3K27ac_cov <- make_cov_plot(H3K27ac_cov,110583000,110641000,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K27ac"),"Depth coverage")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_transcription_factor_cov <- make_cov_plot(transcription_factor_cov,110583000,110641000,max_cov_plot)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("ERG"),"Depth coverage")))}
    if (!(cell %in% c('CHLA25','EW3','TC106'))){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("FLI1"),"Depth coverage")))}

    # Macs2 plots
    all_peaks <- rbind(rbind(H3K4me3_macs2,transcription_factor_macs2),H3K27ac_macs2)
    all_peaks <- subset(all_peaks, V2>=110583000 & V3 <= 110641000)
    min_gradient <- floor(min(all_peaks$V9))
    max_gradient <- ceiling(max(all_peaks$V9))

    legend_macs2 <- make_macs2_legend(H3K4me3_macs2,110583000,110641000,min_gradient,max_gradient)
    plot_macs2_H3K4me3 <- make_macs2_plot(H3K4me3_macs2,110583000,110641000,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K4me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_macs2_H3K27ac <- make_macs2_plot(H3K27ac_macs2,110583000,110641000,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K27ac"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_macs2_transcription_factor <- make_macs2_plot(transcription_factor_macs2,110583000,110641000,min_gradient,max_gradient)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("ERG"),"Macs2")))}
    if (!(cell %in% c('CHLA25','EW3','TC106'))){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("FLI1"),"Macs2")))}

    # Merging plots
    pdf(paste(c("plots/KCNA2_",cell,"_three_seq.pdf"),collapse=""), width = 7.87, height = 8.27)
    print(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label=cell) + theme_void(),
                ggplot() + xlim(110583000,110641000) + annotate("text",x=max(GGAA_coords$start), y = 1, size=3, label="(GGAA)[n]",parse=TRUE) + theme_void(),
                plot_H3K4me3_cov,plot_macs2_H3K4me3,
                plot_H3K27ac_cov, plot_macs2_H3K27ac,
                plot_transcription_factor_cov,plot_macs2_transcription_factor,
                plot_annot, align = "v", axis="tb", ncol=1, nrow=9, rel_heights = c(0.04,0.02,rep(c(0.16, 0.05),3),0.31)),
            ggdraw(),legend_macs2, ncol=3, rel_widths=c(0.85,0.05,0.15)))
    dev.off()

    ## For a plot with only transcription factor and H3K4me3 sequencings
    # Coverage plots
    max_cov_plot <- max(c(H3K4me3_cov$reads_all,H3K27ac_cov$reads_all,transcription_factor_cov$reads_all))
    plot_H3K4me3_cov <- make_cov_plot(H3K4me3_cov,110583000,110641000,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K4me3"),"Depth coverage"))) +
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_transcription_factor_cov <- make_cov_plot(transcription_factor_cov,110583000,110641000,max_cov_plot)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("ERG"),"Depth coverage")))}
    else {plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("FLI1"),"Depth coverage")))}

    # Macs2 plots
    all_peaks <- rbind(H3K4me3_macs2,transcription_factor_macs2)
    all_peaks <- subset(all_peaks, V2>=110583000 & V3 <= 110641000)
    min_gradient <- floor(min(all_peaks$V9))
    max_gradient <- ceiling(max(all_peaks$V9))

    legend_macs2 <- make_macs2_legend(H3K4me3_macs2,110583000,110641000,min_gradient,max_gradient)
    plot_macs2_H3K4me3 <- make_macs2_plot(H3K4me3_macs2,110583000,110641000,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K4me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    plot_macs2_transcription_factor <- make_macs2_plot(transcription_factor_macs2,110583000,110641000,min_gradient,max_gradient)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#dc056b", alpha=1)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("ERG"),"Macs2")))}
    else {plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("FLI1"),"Macs2")))}

    # Merging plots
    pdf(paste(c("plots/KCNA2_",cell,"_two_seq.pdf"),collapse=""), width = 7.87, height = 9.84)
    print(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label=cell) + theme_void(),
                ggplot() + xlim(110583000,110641000) + annotate("text",x=max(GGAA_coords$start), y = 1, size=3, label="(GGAA)[n]",parse=TRUE) + theme_void(),
                plot_H3K4me3_cov,plot_macs2_H3K4me3,
                plot_transcription_factor_cov,plot_macs2_transcription_factor,
                plot_annot, align = "v", axis="tb", ncol=1, nrow=7, rel_heights = c(0.03,0.02,rep(c(0.25, 0.075),2),0.3)),
            ggdraw(),legend_macs2, ncol=3, rel_widths=c(0.85,0.05,0.15)))
    dev.off()
}

