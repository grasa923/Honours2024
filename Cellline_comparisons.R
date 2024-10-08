##Date: 28/08/2024
##Author; Sarah Grant
##Description: Comparing MCF7 and MCF7 het5 cells, and MCF10 and MCF10del to identify up-regulated genes from PTC causing variants.

## Cell comparison

#library(pathfindR)
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(paletteer)
library(cowplot)

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

##Sub-setting for untreated MCF-7 and Het-5 cell line
dds_mcf7 <- DESeqDataSetFromMatrix(raw_counts[,sampleSheet$Treatment == "Untreated" & grepl("MCF_7",sampleSheet$Cell.line)], 
                                   sampleSheet[sampleSheet$Treatment == "Untreated" & grepl("MCF_7",sampleSheet$Cell.line),],
                                   design = ~Cell.line) 
#PCA plot for MCF-7
vsd <- vst(dds_mcf7)
pca <- prcomp(t(assay(vsd)))
p.df <- as.data.frame(pca$x[,1:2])
p.df$Cell.line <- sampleSheet$Cell.line[sampleSheet$Treatment == "Untreated" & grepl("MCF_7",sampleSheet$Cell.line)]

p2 <- ggplot(p.df, aes(PC1, PC2, color=Cell.line)) + geom_point() + theme_bw()

#filtering lowly expressed genes
keep <- rowSums(counts(dds_mcf7) >=10) >=4
dds_mcf7 <- dds_mcf7[keep,]

dds_mcf7 <- DESeq(dds_mcf7)


## do pairwise DE analysis
mcf7 <- results(dds_mcf7)

## shrink log fold change estimates (recommended) - more accurate for lowly expressed genes
mcf7_shrink <- lfcShrink(dds_mcf7, type='ashr', coef="Cell.line_MCF_7_Het5_vs_MCF_7")

#Volcano plot for MCF-7 data
p.df <- as.data.frame(mcf7) %>% filter(!is.na(padj))
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -2 | p.df$log2FoldChange > 2), T,F)
ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color))+geom_point() + geom_vline(xintercept=c(-2,2), linetype='dotted') +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("MCF7 vs MCF7 BRCA1-/+ untreated")


## Sub-setting for untreated MCF10A and MCF10 del185AG
dds_mcf10 <- DESeqDataSetFromMatrix(raw_counts[,grepl("MCF_10",sampleSheet$Cell.line)], 
                                    sampleSheet[grepl("MCF_10",sampleSheet$Cell.line),],
                                    design = ~Cell.line) 
#PCA plot for MCF10 
vsd <- vst(dds_mcf10)
pca <- prcomp(t(assay(vsd)))
p.df <- as.data.frame(pca$x[,1:2])
p.df$Cell.line <- sampleSheet$Cell.line[grepl("MCF_10",sampleSheet$Cell.line)]


p3 <- ggplot(p.df, aes(PC1, PC2, color=Cell.line)) + geom_point() + theme_bw() 


## Remove lowly expressed genes - 10 reads in at least 6 samples
keep <- rowSums(counts(dds_mcf10) >=10) >=6
dds_mcf10 <- dds_mcf10[keep,]

dds_mcf10 <- DESeq(dds_mcf10)

## do pairwise DE analysis
mcf10 <- results(dds_mcf10)


## shrink log foldchange estimates (recommended) - more accurate for lowly expressed genes
mcf10_shrink <- lfcShrink(dds_mcf10, type='ashr', coef="Cell.line_MCF_10A_del185AG_vs_MCF_10A")

p.df <- as.data.frame(mcf10) %>% filter(!is.na(padj))
p.df$color <- ifelse(p.df$padj <0.05 & (p.df$log2FoldChange < -2 | p.df$log2FoldChange > 2), T,F)
ggplot(p.df, aes(x=log2FoldChange, y=-log10(pvalue), color = color))+geom_point() + geom_vline(xintercept=c(-2,2), linetype='dotted') +
  geom_hline(yintercept=-log10(max(p.df$pvalue[p.df$padj <0.05], na.rm=T)), linetype='dotted')  + theme_bw() + scale_colour_manual(values=c('grey','red'))+
  theme(legend.position='none') + ggtitle("MCF10A vs MCF10A BRCA1-/+ untreated")



### Overlapping DE genes
merged <- merge(as.data.frame(mcf7_shrink),as.data.frame(mcf10_shrink), by="row.names", all=T, suffixes = c(".MCF7", ".MCF10"))

symbols <- mapIds(org.Hs.eg.db, keys = gsub("\\.[0-9]+","", merged$Row.names), keytype = "ENSEMBL", column="SYMBOL")

merged$symbol <- symbols


merged %>% filter(padj.MCF7 < 0.05, padj.MCF10 < 0.05) %>% dplyr::select(symbol, log2FoldChange.MCF7, log2FoldChange.MCF10, padj.MCF7, padj.MCF10) %>% mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>% knitr::kable(col.names = c("Gene", "MCF-7", "MCF-10A", "MCF-7", "MCF-10A"), digits=2, align= c('l','r','r','r','r'))  %>%kableExtra::kable_paper(full_width=T) %>%
  kableExtra::add_header_above(c(" " = 1, "LogFoldChange" = 2, "pvalue (adjusted)" = 2))%>%
  kableExtra::scroll_box(width = "100%", height = "200px")

#write file
write.csv(as.data.frame(merged)%>% filter(!is.na(padj.MCF7))%>% filter(!is.na(padj.MCF10)), 
  file="merged.csv")
