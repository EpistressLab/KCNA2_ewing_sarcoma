######### Packages ########

library(ggplot2)
library(tidyverse)

######### Params ########

palette <- c("#68BD8E","red")
names(palette) <- c("FALSE","TRUE")

# genes <- c('KCNA2','KCNA1','KCNA3','KCNA4','KCNA5','KCNA6','KCNA7','KCNA10')
# KCNA6 is not available
 genes <- c('KCNA2','KCNA1','KCNA3','KCNA4','KCNA5','KCNA7','KCNA10')

for (gene in genes){
    ######### Get data #########

    exp <- read.table("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", h=T)
    samples <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t", h=T, quote="")[,c(1,7)]

    exp <- as.data.frame(t(exp[which(exp$Description == gene),]))
    exp$SAMPID <- gsub("\\.","-",rownames(exp))
    exp <- exp[-c(1,2),]
    exp$log2TPM <- log(as.numeric(exp[,1])+1, 2)
    exp <- merge(exp,samples)


    ###### Adding ewing sarcoma data #####

    exp_CCLE <- read.table("../CCLE_Depmap_22Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv", sep=",", h=T)
    samples_CCLE <- read.table("../CCLE_Depmap_22Q4/Model.csv", sep=",", h=T)[,c(1,23)]

    colnames(exp_CCLE) <- sapply(strsplit(colnames(exp_CCLE), "..", fixed=T),'[[',1)
    colnames(exp_CCLE)[1] <- "ModelID"

    exp_CCLE <- exp_CCLE[,colnames(exp_CCLE) %in% c("ModelID",gene)]
    data_exp_CCLE <- merge(exp_CCLE,samples_CCLE)
    colnames(data_exp_CCLE) <- c("SAMPID","log2TPM","SMTSD")
    data_exp <- rbind(exp[,c(1,3,4)],data_exp_CCLE[which(data_exp_CCLE$SMTSD == "Ewing Sarcoma"),])
    data_exp$is_ewing <- as.character(data_exp$SMTSD == "Ewing Sarcoma")

    png(paste(c("GTEx_boxplot_",gene,"_with_ewing.png"),collapse=""),height=21,width=13,units="cm",res=300)
    print(ggplot(data_exp, aes(x=log2TPM, y=reorder(SMTSD,log2TPM,decreasing=F), fill=is_ewing)) +
        geom_boxplot() +
        scale_fill_manual(values=palette) +
        xlab(bquote(log[2]("TPM+1"))) +
        ylab("Tissue Site Detail field") +
        theme_minimal() +
        ggtitle("GTEx version 8 and\nEwing Sarcoma DepMap 22Q4") +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill="none") +
        labs(subtitle=bquote(italic(.(gene))~"expression")))
    dev.off()
}