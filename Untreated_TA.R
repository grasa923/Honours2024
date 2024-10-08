##Date: 2/09/2024
##Author; Sarah Grant
##Description: Compare untreated MCF-7 cell lines with RNA treated MCF-7 samples

library(DESeq2)
#library(pathfindR)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(paletteer)
library(cowplot)
library(EnhancedVolcano)
library(tidyverse)
library(biomaRt)


#General sheet cleanup
sampleSheet <- read.delim("V350134869_sampleSheet.txt")
raw_counts <- read.delim("V350134869_raw_count_matrix.txt", row.names = 1)

raw_counts <- raw_counts[, match(sampleSheet$filename, colnames(raw_counts))]

sampleSheet$Treatment <- gsub("-| ", "_", sampleSheet$Treatment)
sampleSheet$Cell.line <- gsub("-| ", "_", sampleSheet$Cell.line)
sampleSheet$condition <- paste0(sampleSheet$Cell.line,sampleSheet$Treatment)
sampleSheet$Batch <- gsub("#","",stringr::str_extract(sampleSheet$Sample, "#[0-9]"))
sampleSheet <- sampleSheet %>% mutate(condition = gsub("Untreated", "", condition),
                                      condition = gsub("MCF_", "MCF", condition),
                                      condition = gsub("_"," ", condition),
                                      condition = gsub("RNA","", condition),
                                      condition = gsub("7C","7 C", condition),
                                      condition = gsub("7U","7 U", condition))


raw_counts <- raw_counts[,sampleSheet$Cell.line != "unknown"]
sampleSheet <- sampleSheet[sampleSheet$Cell.line != "unknown",]

prac_sampleSheet <- sampleSheet[!grepl("MCF_10",sampleSheet$Cell.line) & !grepl("MCF_7_Het5",sampleSheet$Cell.line),]

#Transcriptional Adaptation
# Design with MCF7 cell lines
dds_MB <- DESeqDataSetFromMatrix(raw_counts[rownames(raw_counts)!="ENSG00000012048.25",!grepl("MCF_10",sampleSheet$Cell.line) & !grepl("MCF_7_Het5",sampleSheet$Cell.line)], 
                                   sampleSheet[!grepl("MCF_10",sampleSheet$Cell.line) & !grepl("MCF_7_Het5",sampleSheet$Cell.line),],
                                   design = ~Treatment) 

## Remove lowly expressed genes - 10 reads in at least 4 samples
keep <- rowSums(counts(dds_MB) >=10) >=4
dds_MB <- dds_MB[keep,]

dds_MB <- DESeq(dds_MB)

## do pairwise DE analysis
untreated_capped <- results(dds_MB, contrast = c("Treatment", "Capped_BRCA1_RNA","Untreated"))
untreated_uncapped <- results(dds_MB, contrast = c("Treatment", "Uncapped_BRCA1_RNA","Untreated"))

## shrink log foldchange estimates (recommended) - more accurate for lowly expressed genes
untreated_capped_shrink <- lfcShrink(dds_MB, contrast = c("Treatment", "Untreated","Capped_BRCA1_RNA"), type='ashr')
untreated_uncapped_shrink <- lfcShrink(dds_MB, contrast = c("Treatment", "Untreated","Uncapped_BRCA1_RNA"), type='ashr')

#PCA plot
vsd <- vst(dds_MB)

pca <- prcomp(t(assay(vsd)))
plotPCA(vsd, intgroup=c("Treatment","Cell.line"))

#Volcano  plot
p.df <- as.data.frame(untreated_capped) %>% filter(!is.na(padj))
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -2 | p.df$log2FoldChange > 2), T,F)
p7<-ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color))+geom_point() + geom_vline(xintercept=c(-2,2), linetype='dotted') +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("Capped_BRCA1_RNA vs Untreated")

p.df <- as.data.frame(untreated_uncapped) %>% filter(!is.na(padj))
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -2 | p.df$log2FoldChange > 2), T,F)
p8<-ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color))+geom_point() + geom_vline(xintercept=c(-2,2), linetype='dotted') +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("Uncapped_BRCA1_RNA vs Untreated")

#Individual Gene Plots
plotCounts(dds_MB, gene = 'ENSG00000151150.22', intgroup = "Treatment")
plotCounts(dds_TA, gene = 'ENSG00000126562.17', intgroup = "Treatment")

#Write to file
write.csv(as.data.frame(untreated_capped)%>% filter(!is.na(padj)), 
          file="untreated_capped.csv")

write.csv(as.data.frame(untreated_uncapped)%>% filter(!is.na(padj)), 
          file="untreated_uncapped.csv")

#converting EnsembleID to Gene name
untreated_uncapped_file <- read_csv("untreated_uncapped.csv")

untreated_uncapped_file$GeneID <- as.character(untreated_uncapped_file$GeneID)
untreated_uncapped_file$GeneID <- sub("[.][0-9]*","",untreated_uncapped_file$GeneID)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gencode_ids <- untreated_uncapped_file$GeneID

conversion <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = gencode_ids,
  mart = mart
)

untreated_uncapped_symbol <- left_join(untreated_uncapped_file, conversion, by = c("GeneID" = "ensembl_gene_id"))

write.csv(untreated_uncapped_symbol,
          file = "untreated_uncapped_symbol.csv")


