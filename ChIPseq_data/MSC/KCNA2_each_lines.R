
library(ggplot2)
library(gggenes)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)

##################### Functions ###############

get_density_pos <- function(pos,data) {
  density_entry <- data[data$V2 <= pos & data$V3 > pos, "V4"]
  if (length(density_entry) == 1) return(density_entry)
  else return(NA)
}

make_df_cov <- function(file_path,chrom,start,end){
    data_density <- subset(read.table(file_path, sep="\t", h=F),V1 == chrom)
    all_density <- sapply(start:end, get_density_pos, data=data_density)
    return(data.frame("chromosome"=rep(chrom,length(start:end)), position=start:end, reads_all=all_density))
}

make_cov_plot <- function(data_cov,start,end,color_fill,ymax){
  plot <- ggplot(data_cov, aes(x=position)) +
  theme_minimal() +
  theme(legend.position="none") %+replace% theme(
    axis.ticks.y = element_blank(),
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9)) +
  xlim(start,end) +
  ylim(0,ymax) +
  xlab("Position on chr1") +
  geom_line(aes(y = reads_all), linewidth=0.5) +
  geom_area(aes(y = reads_all), linewidth=0.05,linetype = 0, fill=color_fill, alpha=0.5, outline.type="full")
  return(plot)
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
gtf_KCNA2$nickname <- factor(gtf_KCNA2$nickname, levels=c('H','G','F','E','D','C','B','A'))

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

# plasmid control
MSC_control_FLI1_data <- make_df_cov("./Samples/MSC1_control_FLI1/bedgraphs/MSC1_control_FLI1_hg38.bedGraph","chr1",110583000,110641000)
MSC_control_H3K27ac_data <- make_df_cov("./Samples/MSC1_control_H3K27ac/bedgraphs/MSC1_control_H3K27ac_hg38.bedGraph","chr1",110583000,110641000)

# plasmid EWSFLI1
MSC_EWSFLI1_FLI1_data <- make_df_cov("./Samples/MSC1_EWSFLI1_FLI1/bedgraphs/MSC1_EWSFLI1_FLI1_hg38.bedGraph","chr1",110583000,110641000)
MSC_EWSFLI1_H3K27ac_data <- make_df_cov("./Samples/MSC1_EWSFLI1_H3K27ac/bedgraphs/MSC1_EWSFLI1_H3K27ac_hg38.bedGraph","chr1",110583000,110641000)

# plasmid FLI1
MSC_FLI1_FLI1_data <- make_df_cov("./Samples/MSC1_FLI1_FLI1/bedgraphs/MSC1_FLI1_FLI1_hg38.bedGraph","chr1",110583000,110641000)
MSC_FLI1_H3K27ac_data <- make_df_cov("./Samples/MSC1_FLI1_H3K27ac/bedgraphs/MSC1_FLI1_H3K27ac_hg38.bedGraph","chr1",110583000,110641000)

# All plasmids 
max_cov_plot <- max(c(MSC_control_FLI1_data$reads_all,MSC_control_H3K27ac_data$reads_all,MSC_EWSFLI1_FLI1_data$reads_all,MSC_EWSFLI1_H3K27ac_data$reads_all,MSC_EWSFLI1_FLI1_data$reads_all,MSC_EWSFLI1_H3K27ac_data$reads_all),na.rm=T)

MSC_control_FLI1 <- make_cov_plot(MSC_control_FLI1_data,110583000,110641000,"#85d1e4",max_cov_plot) + 
  ylab(bquote(atop(bold("Plasmid control"),"FLI1"))) +
  geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
MSC_control_H3K27ac <- make_cov_plot(MSC_control_H3K27ac_data,110583000,110641000,"#85d1e4",max_cov_plot) + 
  ylab(bquote(atop(bold("Plasmid control"),"H3K27ac"))) +
  geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
MSC_EWSFLI1_FLI1 <- make_cov_plot(MSC_EWSFLI1_FLI1_data,110583000,110641000,"#ec93a6",max_cov_plot) + 
  ylab(bquote(atop(bold("Plasmid EWSFLI1"),"FLI1"))) +
  geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
MSC_EWSFLI1_H3K27ac <- make_cov_plot(MSC_EWSFLI1_H3K27ac_data,110583000,110641000,"#ec93a6",max_cov_plot) + 
  ylab(bquote(atop(bold("Plasmid EWSFLI1"),"H3K27ac"))) +
  geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
MSC_FLI1_FLI1 <- make_cov_plot(MSC_FLI1_FLI1_data,110583000,110641000,"#f49b70",max_cov_plot) + 
  ylab(bquote(atop(bold("Plasmid FLI1"),"FLI1"))) +
  geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)
MSC_FLI1_H3K27ac <- make_cov_plot(MSC_FLI1_H3K27ac_data,110583000,110641000,"#f49b70",max_cov_plot) + 
  ylab(bquote(atop(bold("Plasmid FLI1"),"H3K27ac"))) +
  geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#dc056b", alpha=1)

pdf("./plots/KCNA2_MSC1.pdf", width = 7.87, height = 9.84)
plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label="MSC with plasmids, ChIPseq H3K27ac and FLI1") + theme_void(),
            ggplot() + xlim(110583000,110641000) + annotate("text",x=max(GGAA_coords$start), y = 1, size=3, label="(GGAA)[n]",parse=TRUE) + theme_void(),
            MSC_control_FLI1 ,MSC_control_H3K27ac,
            MSC_FLI1_FLI1 ,MSC_FLI1_H3K27ac,
            MSC_EWSFLI1_FLI1 ,MSC_EWSFLI1_H3K27ac,
            plot_annot, align = "v", axis="tb", ncol=1, nrow=9, rel_heights = c(0.04,0.02,rep(c(0.124, 0.124),3),0.25))
dev.off()