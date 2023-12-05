############ Packages ############

require(DESeq2)
library(biomaRt)
library(genefilter,quietly=TRUE)
library(RColorBrewer)
library(gplots)
library(fdrtool)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(pvclust)
library(circlize)
library(GO.db)
library(GSEABase)
library(clusterProfiler)
library(pathview)
library(DOSE)
library(fgsea)
library(IHW)
library(org.Hs.eg.db)
library(tximport)
library(tidyverse)

############ Arguments & data ############

INFILE <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/RNAseq_ewing/deseq2_kallisto/kallisto_conditions.tab"
OUTDIR <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/RNAseq_ewing/deseq2_kallisto"
BIOMART <- "feb2014.archive.ensembl.org,ENSEMBL_MART_ENSEMBL,hsapiens_gene_ensembl"
ASSEMBLY <- "hg38"
CORRANNOT <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/scripts/corresIDorg.txt"
LOGFCTHRESHOLD <- 1
FDRTHRESHOLD <- 0.05

##################################
# 			with kallisto 		#
##################################

sampleTable <- read.table(INFILE, header=TRUE,sep = '\t')

conditions<-factor( sampleTable[ , 3] )
uniq_conds <- unique(conditions)

corresIDorg <- read.table(CORRANNOT,header=T,row.names=1,sep = '\t')
orgGO <- as.character(corresIDorg[ASSEMBLY,1])
orgKegg <- as.character(corresIDorg[ASSEMBLY,2])

files <- sampleTable[,2]
names(files) <- sampleTable[,1]
txi.kallisto <- tximport(files, type="kallisto", txOut=T)
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
dds <- DESeq(dds)

# Normalized counts matrix
matrix <- counts(dds,normalized=TRUE)
# Transforming raw counts
rld <- rlogTransformation(dds, blind=TRUE)
mrld <- assay(rld)

## Component plot
rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pca = prcomp(t(assay(rld)[select,]))
  
variance = pca$sdev^2 / sum(pca$sdev^2)
variance = round(variance, 3) * 100
data_pca = as.data.frame(pca$x)
data_pca$samplename = rownames(data_pca)
data_pca = merge(data_pca,sampleTable[,c(1,3)])

outPlot = c(paste(OUTDIR,"/PCAplot_without_names.png",sep=""))
png(outPlot,,width=20,height=15,units="cm",res=300)
ggplot(data_pca, aes(x=PC1,y=PC2, label=samplename, color=condition)) +
	geom_point(size=5) +
	theme_light() +
	xlab(paste(c("PC1 (variance = ",variance[1],"%)"),collapse="")) +
	ylab(paste(c("PC2 (variance = ",variance[2],"%)"),collapse="")) +
	guides(color=guide_legend(title="Condition"))
dev.off()

########## Boxplot genes ########

fpkm_exp <- fpkm(dds)
palette_ASP14 <- c("#ece0c2","#090093")
names(palette_ASP14) <- c('ASP14_day0','ASP14_day7')

gene_selected <- "KCNA2"
#transcripts_selected <- correspondance_symbol[which(correspondance_symbol$symbol==gene_selected),1]
transcripts_selected <- c("NM_004974.4","XM_011541396.3","XM_017001213.2","XM_011541400.3","NM_001204269.2","XM_011541398.3")
count_norm <- fpkm_exp
data_boxplot <- t(count_norm[which(rownames(count_norm) %in% transcripts_selected),])
data_boxplot_2 <- data.frame('expression'=c(),'sample'=c(),'condition'=c(),'transcript'=c())

for (transcript in transcripts_selected){
  sub_data <- data.frame('expression'= as.numeric(data_boxplot[,transcript]),'sample'=sampleTable[,1],'condition'=sampleTable[,3],'transcript'=transcript)
  data_boxplot_2 <- rbind(data_boxplot_2,sub_data)
}
data_boxplot <- data_boxplot_2
transcripts_nicknames <- data.frame('transcript'=c("NM_004974.4","XM_011541396.3","XM_017001213.2","XM_011541400.3","NM_001204269.2","XM_011541398.3"),'nickname'=c('F','E','D','C','B','A'))
data_boxplot <- merge(data_boxplot,transcripts_nicknames)

cells <- c('A673','ASP14_day0','ASP14_day7','EW24','EW3','MHHES1','RDES','SKES1','TC71')
for (cell in cells){
  outPlot = c(paste(c(OUTDIR,"/KCNA2_expression_",cell,".png"),collapse=""))
  png(outPlot, width=14, height=12, units="cm", res=300)
  print(ggplot(subset(data_boxplot,condition==cell), aes(x = nickname, y = expression, fill=transcript)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal() +
    guides(fill = "none", size="none") +
    xlab(expression(italic(KCNA2)~" transcripts")) +
    ylab(expression("Normalized expression of "~italic(KCNA2)~" transcripts (FPKM)")) +
    ggtitle(cell) +
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

palette_ASP14 <- c("#ece0c2","#090093")
names(palette_ASP14) <- c('ASP14_day0','ASP14_day7')
outPlot = c(paste(c(OUTDIR,"/KCNA2_expression_ASP14.png"),collapse=""))
png(outPlot, width=18, height=12, units="cm", res=300)
print(ggplot(subset(data_boxplot,condition %in% c('ASP14_day0','ASP14_day7')), aes(x = nickname, y = expression, fill=condition)) +
  geom_boxplot(alpha=0.65) +
  labs(fill="Condition") +
  theme_minimal() +
  guides(size="none") +
  scale_fill_manual(values=palette_ASP14, labels=c('ASP14 day 0','ASP14 day 7')) +
  xlab(expression(italic(KCNA2)~" transcripts")) +
  ylab(expression("Normalized expression of "~italic(KCNA2)~" transcripts (FPKM)")) +
  ggtitle("ASP14") +
  theme(plot.title = element_text(hjust = 0.5)))
dev.off()

######## ASP14 comparison ########
condCol="condition"
logFCthreshold=1
AdjPValthreshold=0.05
GenesInFig=50
bootstrap=FALSE
nboot=30

cond1 <- "ASP14_day0"
cond2 <- "ASP14_day7"

comp <- paste(cond1,cond2,sep="__vs__")

sampleAnnot <- as.data.frame(sampleTable[,3])
rownames(sampleAnnot) <- sampleTable[,2]
colnames(sampleAnnot) <- "condition"
samples <- rownames(sampleAnnot)

res <- results(dds,contrast=c(condCol,cond1,cond2),independentFiltering=T)
res$meanInComp <- rowMeans(mrld[,sampleAnnot[,condCol]%in% c(cond1,cond2)])
res$padj <- adj_pvalues(ihw(pvalues = res$pvalue,covariates = res$meanInComp,alpha=AdjPValthreshold))

DE <- data.frame(res)

DE.sel <- list()
DE.sel$up <- DE[which(DE$padj < AdjPValthreshold & DE$log2FoldChange > logFCthreshold),]
DE.sel$down <- DE[which(DE$padj < AdjPValthreshold & DE$log2FoldChange < -logFCthreshold),]

DE.sel$isDE=rbind(DE.sel$up,DE.sel$down)
DE.sel$notDE=DE[setdiff(rownames(DE),rownames(DE.sel$isDE)),]

DE$DE="NONE"
DE[rownames(DE.sel$up),"DE"]="UP"
DE[rownames(DE.sel$down),"DE"]="DOWN"
DE$DE=factor(DE$DE,levels=c("DOWN","NONE","UP"))
DE$transcript <- rownames(DE)
correspondance_symbol <- read.table("/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/genomes/hg38/table_symbol_refseq.tsv",sep="\t",h=T)
DE <- merge(DE,correspondance_symbol)
write.table(DE, paste(paste(OUTDIR,comp,sep="/"),"DE.tsv",sep="_"),sep="\t",row.names=F,quote=F)


DE$DE_label <- NA
list_genes <- c("KCNA2")
DE$DE_label[which(DE$gene %in% list_genes)] <- DE$transcript[which(DE$gene %in% list_genes)]
DE_palette = c("firebrick3", "gray60", "palegreen3")
names(DE_palette) = c("DOWN","NONE","UP")

png(paste(c(OUTDIR,"/",comp,"_Volcano_with_names.png"),collapse=""),width=25,height=25,units="cm",res=300)
  ggplot(data=DE, aes(x=log2FoldChange, y=-log10(padj), color=DE, label=DE_label)) + 
  geom_point() +
  theme_minimal() +
  geom_text_repel(color="black", min.segment.length=0, nudge_x=2) +
  scale_color_manual(values=DE_palette)
dev.off()
