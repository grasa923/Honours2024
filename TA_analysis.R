##Date: 28/08/2024
##Author; Sarah Grant
##Description: Sub-setting RNA-seq data isolated from MCF-7 cells treated with BRCA1 RNA to analyse 
#differentially expressed genes for transcriptional adaptation response

library(DESeq2)
#library(pathfindR)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(paletteer)
library(cowplot)
library(EnhancedVolcano)


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

#Transcriptional Adaptation
# Transcriptional adaptation design, sub-setting for treated MCF-7 cells and including batch as a co variate
## dds_TA is the DESeq object for the counts (minus BRCA1) for the cells transfected
dds_TA <- DESeqDataSetFromMatrix(raw_counts[rownames(raw_counts)!="ENSG00000012048.25",sampleSheet$Treatment != "Untreated"], 
                                 sampleSheet[sampleSheet$Treatment != "Untreated",],
                                 design = ~Treatment+Batch) 

## Remove lowly expressed genes - 10 reads in at least 6 samples
keep <- rowSums(counts(dds_TA) >=10) >=6
dds_TA <- dds_TA[keep,]

dds_TA <- DESeq(dds_TA)

## do pairwise DE analysis
uncapped.B1_capped <- results(dds_TA, contrast = c("Treatment", "Uncapped_BRCA1_RNA","Capped_BRCA1_RNA"))
uncapped.anti_capped <- results(dds_TA, contrast = c("Treatment","Capped_BRCA1_RNA", "Uncapped_anti_sense_RNA"))
uncapped <- results(dds_TA, contrast = c("Treatment", "Uncapped_BRCA1_RNA", "Uncapped_anti_sense_RNA"))

## shrink log foldchange estimates (recommended) - more accurate for lowly expressed genes
uncapped.B1_capped_shrink <- lfcShrink(dds_TA, contrast = c("Treatment", "Uncapped_BRCA1_RNA","Capped_BRCA1_RNA"), type='ashr')
uncapped.anti_capped_shrink <- lfcShrink(dds_TA, contrast = c("Treatment","Capped_BRCA1_RNA", "Uncapped_anti_sense_RNA"), type='ashr')
uncapped_shrink <- lfcShrink(dds_TA, contrast = c("Treatment","Uncapped_BRCA1_RNA", "Uncapped_anti_sense_RNA"), type='ashr')


#Volcano plots
p.df <- as.data.frame(uncapped.B1_capped)
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -2 | p.df$log2FoldChange > 2), T,F)
p1<-ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color))+geom_point() +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("BRCA1 capped vs BRCA1 uncapped")

p.df <- as.data.frame(uncapped.anti_capped) %>% filter(!is.na(padj))
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -1 | p.df$log2FoldChange > 1), T,F)
p2<-ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color))+geom_point() + geom_vline(xintercept=c(-1,1), linetype='dotted') +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("Capped BRCA1  vs Uncapped Anti-sense BRCA1")

p.df <- as.data.frame(uncapped) %>% filter(!is.na(padj))
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -1 | p.df$log2FoldChange > 1), T,F)
p3<-ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color)) +geom_point() + geom_vline(xintercept=c(-1,1), linetype='dotted') +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("Uncapped BRCA1 vs Uncapped Anti-sense BRCA1 ")

#Writing to files
write.csv(as.data.frame(uncapped)%>% filter(!is.na(padj)), 
          file="uncappedvsantisense.csv")

write.csv(as.data.frame(uncapped.anti_capped)%>% filter(!is.na(padj)), 
          file="Cappedvsantisense.csv")
