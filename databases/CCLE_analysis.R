######### Packages ########

library(ggplot2)
library(tidyverse)

######### Params ########

palette <- c("#0000c0","#dc056b")
names(palette) <- c("FALSE","TRUE")

genes <- c('KCNA2','KCNA1','KCNA3','KCNA4','KCNA5','KCNA6','KCNA7','KCNA10')

for (gene in genes){
    ######### Get data #########

    exp <- read.table("OmicsExpressionProteinCodingGenesTPMLogp1.csv", sep=",", h=T)
    samples <- read.table("Model.csv", sep=",", h=T)[,c(1,23)]

    colnames(exp) <- sapply(strsplit(colnames(exp), "..", fixed=T),'[[',1)
    colnames(exp)[1] <- "ModelID"

    exp <- exp[,colnames(exp) %in% c("ModelID",gene)]
    data_exp <- merge(exp,samples)
    colnames(data_exp)[2] <- "expression"
    table_disease <- table(samples$OncotreePrimaryDisease)
    diseases_to_keep <- names(table_disease[which(table_disease!=1)])
    data_exp <- data_exp[which(data_exp$OncotreePrimaryDisease %in% diseases_to_keep),]


    ######## Plot #########

    data_exp$is_ewing <- as.character(data_exp$OncotreePrimaryDisease == "Ewing Sarcoma")

    pdf(paste(c("CCLE_boxplot_",gene,"_v2.pdf"),collapse=""),height=8.25,width=5.9)
    print(ggplot(data_exp, aes(x=expression, y=reorder(OncotreePrimaryDisease,expression,decreasing=F), fill=is_ewing)) +
        geom_boxplot() + 
        scale_fill_manual(values=palette) +
        xlab(bquote(log[2]("TPM+1"))) +
        ylab("Primary disease") +
        theme_minimal() +
        guides(fill="none") +
        ggtitle("DepMap 22Q4") +
        labs(subtitle=bquote(italic(.(gene))~"expression")))
    dev.off()
}

# To have KCNA2 plot with x and y axis inverted
gene <- "KCNA2"
exp <- read.table("OmicsExpressionProteinCodingGenesTPMLogp1.csv", sep=",", h=T)
samples <- read.table("Model.csv", sep=",", h=T)[,c(1,23)]

colnames(exp) <- sapply(strsplit(colnames(exp), "..", fixed=T),'[[',1)
colnames(exp)[1] <- "ModelID"

exp <- exp[,colnames(exp) %in% c("ModelID",gene)]
data_exp <- merge(exp,samples)
colnames(data_exp)[2] <- "expression"

######## Plot #########

data_exp$is_ewing <- as.character(data_exp$OncotreePrimaryDisease == "Ewing Sarcoma")

pdf(paste(c("CCLE_boxplot_",gene,"_v2.pdf"),collapse=""),width=8.25,height=5.9)
    ggplot(data_exp, aes(x=reorder(OncotreePrimaryDisease,expression,decreasing=T), y=expression, fill=is_ewing)) +
    geom_boxplot() + 
    scale_fill_manual(values=palette) +
    ylab(bquote(log[2]("TPM+1"))) +
    xlab("Primary disease") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill="none") +
    ggtitle("DepMap 22Q4") +
    labs(subtitle=bquote(italic(.(gene))~"expression"))
dev.off()