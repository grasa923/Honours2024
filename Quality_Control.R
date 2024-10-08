##Date: 28/08/2024
##Author; Sarah Grant
##Description: Quality control of data isolated from MCF-7 and MCF10A cells 

## General QC
# PCA, MAplot, BRCA1 exprs, 

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

#BRCA expression
b1_exprs <- data.frame(BRCA1=as.numeric(raw_counts[rownames(raw_counts) =="ENSG00000012048.25",]), group = sampleSheet$condition, Treated = ifelse(sampleSheet$Treatment == "Untreated", "Untreated", "Treated"))

b1_exprs %>% mutate(group = gsub("Untreated", "", group),
                    group = gsub("MCF_", "MCF", group),
                    group = gsub("_"," ", group),
                    group = gsub("RNA","", group),
                    group = gsub("7C","T C", group),
                    group = gsub("7U","T U", group))%>%
  ggplot(aes(x=group, y=log10(BRCA1), group=group, fill=Treated)) + geom_boxplot()+theme_bw()+ylab("BRCA1 log10(Raw counts)") + xlab(NULL) + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1) )

# PCA plot to look at clustering based on all genes and removing BRCA1 counts
col.pal <- c()

dds <- DESeqDataSetFromMatrix(raw_counts, 
                              sampleSheet,
                              design = ~condition) 
keep <- rowSums(counts(dds) >=10) >=6
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)


plotMA(resLFC, ylim =c(-2,2))

vsd <- vst(dds)

pca <- prcomp(t(assay(vsd)))
p.df <- as.data.frame(pca$x[,1:2])
p.df$condition <- sampleSheet$condition


p1 <- ggplot(p.df, aes(PC1, PC2, color=condition)) + geom_point() + theme_bw() + scale_color_paletteer_d("IslamicArt::shiraz") + ggtitle("All data")

dds <- DESeqDataSetFromMatrix(raw_counts[rownames(raw_counts)!="ENSG00000012048.25",], 
                              sampleSheet,
                              design = ~condition) 

vsd <- vst(dds)

pca <- prcomp(t(assay(vsd)))
p.df <- as.data.frame(pca$x[,1:2])
p.df$condition <- sampleSheet$condition

p2 <- ggplot(p.df, aes(PC1, PC2, color=condition)) + geom_point() + theme_bw() + scale_color_paletteer_d("IslamicArt::shiraz")+ ggtitle("BRCA1 counts removed")

leg <- get_legend(p1)

plot_grid(p1 + theme(legend.position = 'none'),p2+ theme(legend.position = 'none'),leg, nrow=1, rel_widths = c(1,1,.8))



#Splitting samples into sub experiments
dds_TA <- DESeqDataSetFromMatrix(raw_counts[rownames(raw_counts)!="ENSG00000012048.25",sampleSheet$Treatment != "Untreated"], 
                                 sampleSheet[sampleSheet$Treatment != "Untreated",],
                                 design = ~Treatment) 
keep <- rowSums(counts(dds_TA) >=10) >=6
dds_TA <- dds_TA[keep,]

dds_TA <- DESeq(dds_TA)
res_TA <-results(dds_TA)
plotMA(res_TA, ylim =c(-2,2))

vsd <- vst(dds_TA)

pca <- prcomp(t(assay(vsd)))
p.df <- as.data.frame(pca$x[,1:2])
p.df$Treatment <- sampleSheet$Treatment[sampleSheet$Treatment != "Untreated"]
p.df$Batch <- sampleSheet$Batch[sampleSheet$Treatment != "Untreated"]

p1 <- ggplot(p.df, aes(PC1, PC2, color=Treatment, shape=Batch)) + geom_point() + theme_bw()