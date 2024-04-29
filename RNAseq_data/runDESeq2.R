
############ Functions ############

autoGparFontSizeMatrix<-function(n){ #Calcule automatiquement la taille de police selon le nombre de colonnes ou lignes (empirique)
  n=max(n,50)
  n=min(n,1000)
  return(gpar(fontsize=1/n*600))
}

corrDist<-function(x) return(as.dist((1 - cor(Matrix::t(x)))/2))

rowScale<-function(data,center=TRUE,scaled=FALSE){
  data<-t(data)
  data<-t(scale(data,center=center,scale=scaled))
  return(data)
}

unsupervisedClustering<-function(x,transpose=TRUE,method.dist="pearson",method.hclust="ward.D2",bootstrap=FALSE,nboot=10){
  if(transpose) x<-t(x)
  if(bootstrap){
    require(pvclust)
    resClust<-pvclust(t(x),nboot=nboot,method.hclust = method.hclust,parallel = TRUE,method.dist = method.dist)$hclust
  }else{
    if(method.dist=="pearson"){
      resDist<-corrDist(x)
    }else{
      resDist<-dist(x, method = method.dist)
    }
    resClust<-stats::hclust(resDist,method = method.hclust)
  }
  return(resClust)
}

############ Packages ############

suppressMessages(require("DESeq2"))
suppressMessages(library(biomaRt))
suppressMessages(library(genefilter,quietly=TRUE))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
library(tidyverse)
suppressMessages(library(fdrtool))
suppressMessages(library(ggplot2))
#suppressMessages(library(ggsignif))
suppressMessages(library(ggrepel))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(pvclust))
suppressMessages(library(circlize))
suppressMessages(library(GO.db))
suppressMessages(library(GSEABase))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pathview))
suppressMessages(library(DOSE))
suppressMessages(library(fgsea))
suppressMessages(library(IHW))
suppressMessages(library(org.Hs.eg.db))


############ Arguments & data ############


INFILE <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/rnaseq_YAP/DESEQ2/DESEQ2_conditions.tab"
COUNTDIRECTORY <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/rnaseq_YAP/count_dir"
OUTDIR <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/rnaseq_YAP/DESEQ2/results_all"
BIOMART <- "feb2014.archive.ensembl.org,ENSEMBL_MART_ENSEMBL,hsapiens_gene_ensembl"
ASSEMBLY <- "hg38"
CORRANNOT <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/scripts/corresIDorg.txt"
LOGFCTHRESHOLD <- 1
FDRTHRESHOLD <- 0.05

# Data
sampleTable <- read.table(INFILE, header=TRUE,sep = '\t')
directory <- c(COUNTDIRECTORY)

############ Data preparation ############ 

# List of all the combination possible of paired conditions
conditions<-factor( sampleTable[ , 3] )
uniq_conds <- unique(conditions)

palette <- c("black","#0000c0")
names(palette) <- c('A673_siCT','A673_siKCNA2')

############ Data manipulation ############ 

# Counts recovery
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~condition)

# deseq2 analysis
ddsHTSeq <- DESeq(ddsHTSeq)
matrix <- counts(ddsHTSeq,normalized=TRUE)
rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
mrld <- assay(rld)

############ Plots ############

## Component plot
rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pca = prcomp(t(assay(rld)[select,]))
  
variance = pca$sdev^2 / sum(pca$sdev^2)
variance = round(variance, 3) * 100
data_pca = as.data.frame(pca$x)
data_pca$samplename = rownames(data_pca)
data_pca = merge(data_pca,sampleTable[,c(1,3)])

outPlot = c(paste(OUTDIR,"/PCAplot_with_names.pdf",sep=""))
pdf(outPlot,,width=7.87,height=5.9)
ggplot(data_pca, aes(x=PC1,y=PC2, label=samplename, color=condition)) +
	geom_point(size=4) +
	geom_text_repel(color="black") +
	theme_light() +
	xlab(paste(c("PC1 (variance = ",variance[1],"%)"),collapse="")) +
	ylab(paste(c("PC2 (variance = ",variance[2],"%)"),collapse="")) +
	guides(color=guide_legend(title="Condition")) +
  scale_color_manual(values=palette)
dev.off()

#### Comparison

pw_tests <- list()
pw_tests[[1]] <- c('A673_siKCNA2','A673_siCT')

# Variables
condCol="condition"
logFCthreshold=LOGFCTHRESHOLD
AdjPValthreshold=FDRTHRESHOLD
GenesInFig=50
bootstrap=FALSE
nboot=30

cond1 <- pw_tests[[1]][1]
cond2 <- pw_tests[[1]][2]
comp <- paste(cond1,cond2,sep="__vs__")

sampleAnnot <- as.data.frame(sampleTable[,3])
rownames(sampleAnnot) <- sampleTable[,2]
colnames(sampleAnnot) <- "Condition"
samples <- rownames(sampleAnnot)

res <- results(ddsHTSeq,contrast=c(condCol,cond1,cond2),independentFiltering=T)

## Expression matrices : normalised and transformed
exprDat <- matrix
exprDatT <- mrld

res$meanInComp <- rowMeans(exprDat[,sampleAnnot[,condCol]%in% c(cond1,cond2)])
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
DE$symbol <- rownames(DE)

write.table(DE, file=paste(c(OUTDIR,"/",comp,"_","all_genes.tsv"),collapse=""), sep="\t", row.names=F, quote=F)

DE$DE_label <- NA
#DE$DE_label[which(DE$symbol != "NONE")] <- rownames(DE[which(DE$DE != "NONE"),])
list_genes <- c("KCNA2","CCNE2","CCNDBP1","GLI1","CCN1","CCN2","ITGB2","ITGB3")
DE$DE_label[which(rownames(DE) %in% list_genes)] <- rownames(DE[which(rownames(DE) %in% list_genes),])
DE_palette = c("firebrick3", "gray60", "palegreen3")
names(DE_palette) = c("DOWN","NONE","UP")

#Volcano-plot
pdf(paste(c(OUTDIR,"/",comp,"_Volcano_with_names.pdf"),collapse=""),width=5.9,height=5.9)
  ggplot(data=DE, aes(x=log2FoldChange, y=-log10(padj), col=DE, label=DE_label)) + 
  geom_point(size=0.6) +
  theme_minimal() +
  geom_text_repel(color="black", min.segment.length=0, nudge_x=0.5) +
  scale_color_manual(values=DE_palette)
dev.off()

####### Heatmap cell cycle #########

kegg_gmt <- read.table("KEGG_hsa_pathways_2023.gmt",sep=";",h=F)
get_gene_and_ontology <- function(line){
    line <- unlist(strsplit(line,"\t"))
    this_data <- data.frame("gene"=line[3:length(line)],"vector"=rep(line[1],length(line)-2))
    return(this_data)
}
data_ontologies <- do.call("rbind",apply(kegg_gmt, 1, function(x) get_gene_and_ontology(x)))


these_genes <- unique(data_ontologies[which(data_ontologies[,2] == "Cell cycle"),1])
exprDat <- matrix
exprDatT <- mrld
annot=sampleAnnot["Condition"]
annot[,condCol] <- as.factor(annot[,condCol])
colTopAnnot <- vector("list", ncol(annot))
names(colTopAnnot) <- "Condition"
i<-1
for(col in colnames(annot)){
  colTopAnnot[[col]]<-palette
  i<-i+1
}
ha<-HeatmapAnnotation(df=annot,col = colTopAnnot,name = "Condition")
DEgenes.names=these_genes
sampleHt<-colnames(exprDatT)
haByComp<-ha
exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
exprDE.scaled=na.omit(rowScale(exprDE,center=T,scaled=T))
bootTemp=bootstrap
if(nrow(exprDE.scaled>10)){bootTemp=FALSE}
hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
#hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
kmeans_samples <- kmeans(t(exprDE.scaled), 2)
split <-factor(kmeans_samples$cluster, levels=c(kmeans_samples$cluster[1],kmeans_samples$cluster[4]))
rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01),na.rm = T)
colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
row_font_size <- autoGparFontSizeMatrix(nrow(exprDE.scaled))
pdf(paste(OUTDIR,"Heatmap_cell_cycle.pdf",sep="/"),width=6.7,height=8.66)
draw(Heatmap(exprDE.scaled,top_annotation = haByComp, row_names_gp = row_font_size,
              cluster_rows = hclustGeneDE,col = colHA,name="Z-score",
              column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled)),
              column_split=split,cluster_column_slices = F,column_title = NULL),
              merge_legend = TRUE) 
dev.off()

########## TGF beta signaling pathway ###########

these_genes <- unique(data_ontologies[which(data_ontologies[,2] == "TGF-beta signaling pathway"),1])
exprDat <- matrix
exprDatT <- mrld
annot=sampleAnnot[condCol]
annot[,condCol] <- as.factor(annot[,condCol])
colTopAnnot <- vector("list", ncol(annot))
names(colTopAnnot) <- "Condition"
i<-1
for(col in colnames(annot)){
  colTopAnnot[[col]]<-palette
  i<-i+1
}
ha<-HeatmapAnnotation(df=annot,col = colTopAnnot,name = "Condition")
DEgenes.names=these_genes
sampleHt<-colnames(exprDatT)
haByComp<-ha
exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
exprDE.scaled=na.omit(rowScale(exprDE,center=T,scaled=T))
bootTemp=bootstrap
if(nrow(exprDE.scaled>10)){bootTemp=FALSE}
hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
#hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
kmeans_samples <- kmeans(t(exprDE.scaled), 2)
split <-factor(kmeans_samples$cluster, levels=c(kmeans_samples$cluster[1],kmeans_samples$cluster[4]))
rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01),na.rm = T)
colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
row_font_size <- autoGparFontSizeMatrix(nrow(exprDE.scaled))
pdf(paste(OUTDIR,"Heatmap_tgf_beta.pdf",sep="/"),width=6.7,height=8.66)
draw(Heatmap(exprDE.scaled,top_annotation = haByComp, row_names_gp = row_font_size,
              cluster_rows = hclustGeneDE,col = colHA,name="Z-score",
              column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled)),
              column_split=split,cluster_column_slices = F,column_title = NULL),
              merge_legend = TRUE) 
dev.off()

########## Wnt signaling pathway ###########

these_genes <- unique(data_ontologies[which(data_ontologies[,2] == "Wnt signaling pathway"),1])
exprDat <- matrix
exprDatT <- mrld
annot=sampleAnnot["Condition"]
annot[,condCol] <- as.factor(annot[,condCol])
colTopAnnot <- vector("list", ncol(annot))
names(colTopAnnot) <- "Condition"
i<-1
for(col in colnames(annot)){
  colTopAnnot[[col]]<-palette
  i<-i+1
}
ha<-HeatmapAnnotation(df=annot,col = colTopAnnot,name = "Condition")
DEgenes.names=these_genes
sampleHt<-colnames(exprDatT)
haByComp<-ha
exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
exprDE.scaled=na.omit(rowScale(exprDE,center=T,scaled=T))
bootTemp=bootstrap
if(nrow(exprDE.scaled>10)){bootTemp=FALSE}
hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
#hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
kmeans_samples <- kmeans(t(exprDE.scaled), 2)
split <-factor(kmeans_samples$cluster, levels=c(kmeans_samples$cluster[1],kmeans_samples$cluster[4]))
rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01),na.rm = T)
colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
row_font_size <- autoGparFontSizeMatrix(nrow(exprDE.scaled))
pdf(paste(OUTDIR,"Heatmap_Wnt.pdf",sep="/"),width=6.7,height=8.66)
draw(Heatmap(exprDE.scaled,top_annotation = haByComp, row_names_gp = row_font_size,
              cluster_rows = hclustGeneDE,col = colHA,name="Z-score",
              column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled)),
              column_split=split,cluster_column_slices = F,column_title = NULL),
              merge_legend = TRUE) 
dev.off()

####### Heatmap hippo pathway ########

these_genes <- unique(data_ontologies[which(data_ontologies[,2] == "Hippo signaling pathway"),1])
exprDat <- matrix
exprDatT <- mrld
annot=sampleAnnot["Condition"]
annot[,condCol] <- as.factor(annot[,condCol])
colTopAnnot <- vector("list", ncol(annot))
names(colTopAnnot) <- "Condition"
i<-1
for(col in colnames(annot)){
  colTopAnnot[[col]]<-palette
  i<-i+1
}
ha<-HeatmapAnnotation(df=annot,col = colTopAnnot,name = "Condition")
DEgenes.names=these_genes
sampleHt<-colnames(exprDatT)
haByComp<-ha
exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
exprDE.scaled=na.omit(rowScale(exprDE,center=T,scaled=T))
bootTemp=bootstrap
if(nrow(exprDE.scaled>10)){bootTemp=FALSE}
hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
#hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
kmeans_samples <- kmeans(t(exprDE.scaled), 2)
split <-factor(kmeans_samples$cluster, levels=c(kmeans_samples$cluster[1],kmeans_samples$cluster[4]))
rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01),na.rm = T)
colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
row_font_size <- autoGparFontSizeMatrix(nrow(exprDE.scaled))
pdf(paste(OUTDIR,"Heatmap_hippo.pdf",sep="/"),width=6.7,height=8.66)
draw(Heatmap(exprDE.scaled,top_annotation = haByComp, row_names_gp = row_font_size,
              cluster_rows = hclustGeneDE,col = colHA,name="Z-score",
              column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled)),
              column_split=split,cluster_column_slices = F,column_title = NULL),
              merge_legend = TRUE) 
dev.off()
