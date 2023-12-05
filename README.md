# KCNA2_ewing_sarcoma
## Chipseq data

ChIPseq data are from GSE176400, GSE129155, and GSE94278 from the Gene Expression Omnibus database. GSE176400 contains ChIPseq from 18 Ewing's sarcoma cell lines targeting FLI1 or ERG based on cell line fusion; histone marks H3K4me3, H3K27me3, and H3K27ac; and an input serving as a control. GSE129155 contains ChIPseq targeting the H3K27ac mark of the ASP14 cell line without sh induction targeting EWS-FLI1 (day 0) and after 7 days of sh induction, and ChIPseq targeting FLI1 on the same cell line starting at the end of sh induction (day 7), and then on days 9, 10, 11, 14, and 17; in order to observe the effects of gradual reappearance of EWS-FLI1. These data are supplemented with an input from day 0 that served as a control. GSE94278 contains Bigwig files of normalized ChIPseq values of Mesenchymal Stem Cells (MSCs) transfected with a control plasmid, FLI1 plasmid, and EWS-FLI1 plasmid.

A few corrupted reads (nucleotide sequence length different from the phred score) were removed using Awk. ChIPseq data were all filtered with fastp by removing reads with more than 40% of bases with a phred score below 30 and removing PCR adapters if detected. Reads were aligned with BWA-MEM against the GRCh38.p14 reference genome "analysis set" version from NCBI. After conversion to BAM format, peaks were detected with Macs2 using the controls associated with each lineage. Sequencing depths of the alignment files to construct the graphs were retrieved through homemade python script. Macs2 sequencing depths and peaks coordinates were plotted, with transcript representations from the RefSeq GTF file associated with the reference genome used
As for the number of microsatellite repeats at the EWS-FLI1 or EWS- ERG binding site, the consensus sequence of the region was obtained with bcftools mpileup by keeping only mutations present at least 5 times, then the counting was done by TandemRepeatFinder by strongly penalizing the mismatches and indels.

|  Tool  |   Version  |  Publication  |
|:------:|:----------:|:-------------:|
| Snakemake-minimal | 7.19.1 | https://doi.org/10.12688/f1000research.29032.2 |
| Awk | 4.0.2 | https://doi.org/10.1002/spe.4380090403 |
| fastp | 0.23.2 | https://doi.org/10.1093/bioinformatics/bty560 |
| BWA-MEM | 0.7.17 | https://doi.org/10.1093/bioinformatics/btp324 |
| samtools | 1.16.1 | https://pubmed.ncbi.nlm.nih.gov/33590861/ |
| Macs | 2.2.7.1 | https://doi.org/10.1186/gb-2008-9-9-r137 |
| Pysam | 0.20.0 | https://github.com/pysam-developers/pysam |
| ggplot2 | 3.4.0 | https://doi.org/10.1007/978-0-387-98141-3 |
| gggenes | 0.4.1 | https://github.com/wilkox/gggenes |
| Cowplot | 1.1.1 | https://github.com/wilkelab/cowplot |
| bcftools mpileup | 1.16 | https://pubmed.ncbi.nlm.nih.gov/33590861/ |
| TandemRepeatFinder | 4.09.1 | https://academic.oup.com/nar/article/27/2/573/1061099 |


## RNAseq data

RNAseq data from the A-673, EW-3, EW-24, MHH-ES-1, RD-ES, SK-ES-1, and TC-71 lines were obtained from the ActivMotif Inc. platform. The library preparation and sequencing were performed in this platform. Each of the lines was prepared as biological triplicates. Libraries were prepared using the Illumina TruSeq Stranded mRNA Sample Preparation kit. Then sequencing was performed on Illumina NextSeq 500 in a 2 x 42 bp format (PE42). ASP14 RNAseq data without induction (day 0) and with induction of sh targeting EWS-FLI1 (day 7) came from the Gene Expression Omnibus GSE164373, also in triplicates for each condition. The reads were filtered with fastp 0.20.1 on the same criteria as the ChIPseq data. They were then aligned to the transcriptome associated with RefSeq GRCh38.p14 via kallisto. Counts were normalized with DESeq2. Transcripts were considered differentially expressed if their expression levels between conditions were affected by at least 2-fold, and had an adjusted p-value less than 0.05 (|log2 FC| > 1 & adj.p.value < 0.05), according to a Wald test.


|  Tool  |   Version  |  Publication  |
|:------:|:----------:|:-------------:|
| Snakemake-minimal | 7.19.1 | https://doi.org/10.12688/f1000research.29032.2 |
| kallisto | 0.46.2 | https://doi.org/10.1038/nbt.3519 |
| DESeq2 | 1.30.1 | https://doi.org/10.1186/s13059-014-0550-8 |
