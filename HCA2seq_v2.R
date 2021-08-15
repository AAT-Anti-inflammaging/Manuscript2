setwd("C:/Users/yuany/OneDrive - University of Florida/Song's Lab/Projects/Senescence/RNASeq/DESeq2")
#setwd("C:/Users/yuanyemd/OD/Courses/12. 2019Spring/PHC6937 Genetic Data Analysis/Project/RNASeq/HCA2")
Ctr1 <- read.delim("R053_Ctr1_count.txt", header = F)
Ctr2 <- read.delim("R054_Ctr3_count.txt", header = F)
AAT1 <- read.delim("R055_AAT1_count.txt", header = F)
AAT2 <- read.delim("R056_AAT2_count.txt", header = F)
colnames(Ctr1) <- c("Gene", "Ctr1")
colnames(Ctr2) <- c("Gene", "Ctr2")
colnames(AAT1) <- c("Gene", "AAT1")
colnames(AAT2) <- c("Gene", "AAT2")
library(dplyr)
HCA2seqcount <- full_join(Ctr1,Ctr2)
HCA2seqcount <- full_join(HCA2seqcount,AAT1)
HCA2seqcount <- full_join(HCA2seqcount,AAT2)
#write.csv(HCA2seqcount, file = "HCA2count.csv", row.names = F) 

countdata <- HCA2seqcount
rownames(countdata) <- HCA2seqcount[,1] # gene ID need to be the row name
countdata <- countdata[,-1]
countdata <- countdata[-c(23711:23715),]
colSums(countdata) # the total reads for each sample is different, indicating this is not normalized counts, which can be used for DESeq2
which(is.na(countdata)==TRUE) # no missing value
# this countdata will be the countData argument in DESeqDataSetFromMatrix function
coldata <- data.frame('condition'=c("Control","Control","hAAT","hAAT"))
rownames(coldata) <- c("Ctr1","Ctr2","AAT1","AAT2")
# coldata is the colData argument in DESeqDataSetFromMatrix function

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.8")
# The above 3 lines of code is the right way to install DESeq2 package if your R version is after 3.5.0. Just delete the # and run all 3 lines
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# pre-filtering the row sum smaller than 10

dds$condition <- factor(dds$condition, levels = c("Control","hAAT"))
#  set the reference condition as Control

dds <- DESeq(dds)
# running analysis in single function DESeq
res <- results(dds)
res
# result from the analysis is stored in res

res <- results(dds, name = "condition_hAAT_vs_Control")

resultsNames(dds)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("apeglm", version = "3.8")
# shrinkage of log2 fold change requires installing package apeglm
resLFC <- lfcShrink(dds, coef="condition_hAAT_vs_Control", type = "apeglm")
# shrinked log2 fold change is in resLFC. The apeglm method is used for effect size shrinkage. The shrinkage can help to looks at the largest fold changes that are not due to low counts. The large fold changes from genes with lots of counts would not shrunk, while the imprecise fold changes are shrunk. For publication, need to cite(Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions forsequence count data: removing the noise and preserving large differences.Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895)

resOrdered <- res[order(res$pvalue),]
summary(res)
# order the results table by p-values

res01 <- results(dds, alpha = 0.01)
# change the FDR alpha from 0.1 to 0.05
sum(res01$padj < 0.01, na.rm = T)
# 44 adjusted p-values were less than 0.01

mcols(res01)
# wald test was used to perform inference on coef in negative binomial regression, BH method was used to adjust p-value

plotMA(resLFC, alpha = 0.01)
# plot the shrunken log2 fold changes. The red dots are the ones with p-adj < 0.05

#idx <- identify(res05$baseMean, res05$log2FoldChange)
#significantgenenames <- rownames(res05)[idx]
# After calling plotMA, we can use the function identify to itnteractively detect the row number of individual genes by clicking on the plot. Then save the index.
countdata[which(rownames(countdata)=="NFKBIA"),]
d <- plotCounts(dds, gene=which(rownames(res)=="NFKBIA"), intgroup = "condition", returnData = T)
library(ggplot2)
ggplot(d) + aes(x=condition,y=count, color=condition) + geom_point(size=4) + 
  ggtitle("IkB-alpha") + xlab("") + ylab("Count") + 
  scale_color_manual(name="Condition", values=c("orange","blue")) + theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text =  element_text(size = 16),
        title = element_text(size = 20),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5))
# plot the gene with smallest p-adj. Can edit plot using ggplot2

#write.csv(as.data.frame(resOrdered), file = "HCA2Seq_DESeq2.csv")
resSig <- subset(resOrdered, padj < 0.01)
#write.csv(as.data.frame(resSig), file = "HCA2Seq_DESeq2_significant.csv")
library(ggrepel)
ggplot(as.data.frame(resSig)) + aes(x=-log(padj),y=log2FoldChange) + geom_jitter(color = "orange", size = 2) + 
  xlab("-log(Adjusted P-value)") + ylab("Log2(Fold Change)") +
  geom_hline(yintercept=0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_text_repel(aes(label=rownames(resSig)), color = "blue") + ggtitle("DEGs (p-adj < 0.01)") + 
  theme_minimal() + theme(axis.title = element_text(size = 18),
                          axis.text =  element_text(size = 16),
                          title = element_text(size = 20),
                          plot.title = element_text(hjust = 0.5))

resSigup <- subset(resSig, log2FoldChange > 0)
#write.csv(as.data.frame(resSigup), file = "HCA2Seq_DESeq2_significant_upregulation.csv")
resSigdown <- subset(resSig, log2FoldChange < 0)
#write.csv(as.data.frame(resSigdown), file = "HCA2Seq_DESeq2_significant_downregulation.csv")
# output the result into csv file


vsd <- vst(dds, blind = F)
head(assay(vsd), 3)
# variance stabilizing transformation to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low. Blind = F so that the estimated dispersion will be used and no need to re-estimate the dispersion.
# y = log2(n + n0), n: count values, n0: positive constant


library("pheatmap")
siglabel <- which(match(rownames(res), rownames(resSig), nomatch = 1)!=1) # get the row number of the significant genes
sigorder <- order(rowMeans(counts(dds, normalized=T)[siglabel,]), decreasing=T) # get the order of the mean count of significant genes
siglabel <- siglabel[sigorder] # reorder the significant gene row number according to row count
annotation_col <- data.frame(Condition=c("Ctr1", "Ctr2", "AAT1", "AAT2"))
pheatmap(assay(vsd)[siglabel,], cluster_rows = F, show_rownames = F, cluster_cols = F, 
         main = "Normalized Count of DEGs", fontsize = 16) # get heatmap


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db", version = "3.8")
# install the package for enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
# Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers. Carlson M (2018). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.7.0.


DEG <- as.data.frame(resSig)
DEG$gene <- rownames(DEG)
DEG <- DEG %>% arrange(padj, log2FoldChange)
DEG_up <- as.data.frame(resSigup)
DEG_down <- as.data.frame(resSigdown)

hs <- org.Hs.eg.db
symbols <- as.character(rownames(DEG_down))
GOlist <- select(hs, keys = symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
GOlist[c(34, 25),2] <- c("3576", "645638")
GOlist <- GOlist[-27,]
# find the Entrez ID of two NAs and delete the one without Entrez ID

library(DOSE)
ggo_mf <- groupGO(gene = GOlist$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "MF",
               level = 3,
               readable = TRUE)
ggo_bp <- groupGO(gene = GOlist$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  level = 3,
                  readable = TRUE)
barplot(ggo_mf, drop=TRUE) + ggtitle("DEGs Classification Based on Molecular Function")
barplot(ggo_bp, drop=TRUE) + ggtitle("DEGs Classification Based on Biological Process")
# classification of the DEGs into molecular activity (ont=MF) and biological process (ont=BP)


ego_mf <- enrichGO(gene = GOlist$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
ego_bp <- enrichGO(gene = GOlist$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 1e-8)
# Over-representation test(Boyle et al. 2004)
dotplot(ego_mf) + ggtitle("GO Molecular Function")
dotplot(ego_bp) + ggtitle("GO Biological Process")
emapplot(ego_mf) + ggtitle("GO Molecular Function Enrichment Map") + theme(plot.title = element_text(hjust = 0.5))
emapplot(ego_bp) + ggtitle("GO Biological Process Enrichment Map") + theme(plot.title = element_text(hjust = 0.5))
# enrichment map
ego_mf2 <- setReadable(ego_mf, OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
ego_bp2 <- setReadable(ego_bp, OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
foldchange <- DEG_down$log2FoldChange[-27]
# remove the one without ID
names(foldchange) <- GOlist$ENTREZID
cnetplot(ego_mf2, categorySize="pvalue", foldChange=foldchange, readable = TRUE) + ggtitle("DEGs and Molecular Function") + theme(plot.title = element_text(hjust = 0.5))
cnetplot(ego_bp2, categorySize="pvalue", foldChange=foldchange, readable = TRUE) + ggtitle("DEGs and Biological Process") + theme(plot.title = element_text(hjust = 0.5))


#setwd("C:/Users/yuanyemd/OD/Courses/12. 2019Spring/PHC6937 Genetic Data Analysis/Project/RNASeq/HCA2")
cuffdifflist <- read.csv("CuffdiffList.csv", header = T)
cuffdifflist <- c(as.character(cuffdifflist[1:9,1]), as.character(cuffdifflist[,3]))
commongenes <- intersect(DEG$gene, cuffdifflist)
commongenes
length(commongenes)
# 23 genes were found by both Cuffdiff 2 and DESeq2

setwd("C:/Users/yuany/OneDrive - University of Florida/Song's Lab/Projects/Senescence/RNASeq/DESeq2/down-regulation/TranscriptionFactor")
#setwd("C:/Users/yuanyemd/OD/Courses/12. 2019Spring/PHC6937 Genetic Data Analysis/Project/RNASeq/HCA2/down-regulation/TranscriptionFactor")
TF2019 <- read.delim("TRRUST_Transcription_Factors_2019_table.txt", header = T)
humanlabel <- grep("human", TF2019$Term)
TF2019 <- TF2019[humanlabel,]
#TF2019$Overlap <- as.character(TF2019$Overlap)
#TF2019$Overlap <- gsub("/", replacement = "|", TF2019$Overlap)
#write.csv(TF2019, file = "Down_Trust_Transcription_Factor_2019.csv", row.names = F)

TFJASPAR <- read.delim("TRANSFAC_and_JASPAR_PWMs_table.txt", header = T)
humanlabel <- grep("human", TFJASPAR$Term)
TFJASPAR <- TFJASPAR[humanlabel,]
#TFJASPAR$Overlap <- as.character(TFJASPAR$Overlap)
#TFJASPAR$Overlap <- gsub("/", replacement = "|", TFJASPAR$Overlap)
#write.csv(TFJASPAR, file = "Down_TRANSFAC_and_JASPAR_PWMs.csv", row.names = F)

TFPPI <- read.delim("Transcription_Factor_PPIs_table.txt", header = T)
#TFPPI$Overlap <- as.character(TFPPI$Overlap)
#TFPPI$Overlap <- gsub("/", replacement = "|", TFPPI$Overlap)
#write.csv(TFPPI, file = "Down_Transcription_Factor_PPIs.csv", row.names = F)

setwd("C:/Users/yuany/OneDrive - University of Florida/Song's Lab/Projects/Senescence/RNASeq/DESeq2/down-regulation/Pathway")
#setwd("C:/Users/yuanyemd/OD/Courses/12. 2019Spring/PHC6937 Genetic Data Analysis/Project/RNASeq/HCA2/down-regulation/Pathway")
KEGG2019 <- read.delim("KEGG_2019_Human_table.txt", header = T)
#KEGG2019$Overlap <- as.character(KEGG2019$Overlap)
#KEGG2019$Overlap <- gsub("/", replacement = "|", KEGG2019$Overlap)
#write.csv(KEGG2019, file = "Down_KEGG2019.csv", row.names = F)

WikiPathways2019 <- read.delim("WikiPathways_2019_Human_table.txt", header = T)
#WikiPathways2019$Overlap <- as.character(WikiPathways2019$Overlap)
#WikiPathways2019$Overlap <- gsub("/", replacement = "|", WikiPathways2019$Overlap)
#write.csv(KEGG2019, file = "Down_WikiPathways2019.csv", row.names = F)

setwd("C:/Users/yuany/OneDrive - University of Florida/Song's Lab/Projects/Senescence/RNASeq/DESeq2/up-regulation")
KEGG2019up <- read.delim("KEGG_2019_Human_table.txt", header = T)
#KEGG2019up$Overlap <- as.character(KEGG2019up$Overlap)
#KEGG2019up$Overlap <- gsub("/", replacement = "|", KEGG2019up$Overlap)
#write.csv(KEGG2019up, file = "Up_KEGG2019.csv", row.names = F)



tfplot <- TF2019[1:20,]
tftable <- as.table(tfplot$Combined.Score)
names(tftable) <- tfplot$Term
tflevel <- names(tftable)[order(tftable, decreasing = T)]
ggplot(tfplot, aes(x=factor(tflevel, levels = tflevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("Trust Transcription Factor 2019 (Down-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()

jplot <- TFJASPAR[1:20,]
jtable <- as.table(jplot$Combined.Score)
names(jtable) <- jplot$Term
jlevel <- names(jtable)[order(jtable, decreasing = T)]
ggplot(jplot, aes(x=factor(jlevel, levels = jlevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("TFJASPAR (Down-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()


ppiplot <- TFPPI[1:20,]
ppitable <- as.table(ppiplot$Combined.Score)
names(ppitable) <- ppiplot$Term
ppilevel <- names(ppitable)[order(ppitable, decreasing = T)]
ggplot(ppiplot, aes(x=factor(ppilevel, levels = ppilevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("Transcription Factor PPIs (Down-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()


kplot <- KEGG2019[1:20,]
ktable <- as.table(kplot$Combined.Score)
names(ktable) <- kplot$Term
klevel <- names(ktable)[order(ktable, decreasing = T)]
ggplot(kplot, aes(x=factor(klevel, levels = klevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("KEGG Pathway 2019 (Down-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()

wikiplot <- WikiPathways2019[1:20,]
wikitable <- as.table(wikiplot$Combined.Score)
names(wikitable) <- wikiplot$Term
wikilevel <- names(wikitable)[order(wikitable, decreasing = T)]
ggplot(wikiplot, aes(x=factor(wikilevel, levels = wikilevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("Wiki Pathway 2019 (Down-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()


kupplot <- KEGG2019up
kuptable <- as.table(kupplot$Combined.Score)
names(kuptable) <- kupplot$Term
kuplevel <- names(kuptable)[order(kuptable, decreasing = T)]
ggplot(kupplot, aes(x=factor(kuplevel, levels = kuplevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("KEGG Pathway 2019 (Up-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()



setwd("C:/Users/yuany/OneDrive - University of Florida/Song's Lab/Projects/Senescence/RNASeq/DESeq2/down-regulation/GO")
onto <- read.delim("GO_Biological_Process_2018_table.txt", header = T)
#onto$Overlap <- as.character(onto$Overlap)
#onto$Overlap <- gsub("/", replacement = "|", onto$Overlap)
#write.csv(onto, file = "Down_GO_Ontology_Biological_Process_2018.csv", row.names = F)


setwd("C:/Users/yuany/OneDrive - University of Florida/Song's Lab/Projects/Senescence/RNASeq/DESeq2/up-regulation")
ontoup <- read.delim("GO_Biological_Process_2018_table.txt", header = T)
#ontoup$Overlap <- as.character(ontoup$Overlap)
#ontoup$Overlap <- gsub("/", replacement = "|", ontoup$Overlap)
#write.csv(ontoup, file = "Down_GO_Ontology_Biological_Process_2018.csv", row.names = F)


ontoplot <- onto[1:20,]
ontotable <- as.table(ontoplot$Combined.Score)
names(ontotable) <- ontoplot$Term
ontolevel <- names(ontotable)[order(ontotable, decreasing = T)]
ggplot(ontoplot, aes(x=factor(ontolevel, levels = ontolevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("GO Biological Process 2018 (Down-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()


ontoupplot <- ontoup
ontouptable <- as.table(ontoupplot$Combined.Score)
names(ontouptable) <- ontoupplot$Term
ontouplevel <- names(ontouptable)[order(ontouptable, decreasing = T)]
ggplot(ontoupplot, aes(x=factor(ontouplevel, levels = ontouplevel), y=Combined.Score, color=-log(Adjusted.P.value)))  + 
  geom_point(size = 3) + coord_flip() + ylab("Combined Score") + xlab("") + ggtitle("GO Biological Process 2018 (Up-regulated)") + 
  scale_colour_gradient(name = "Negative\nLogarithm\nof Adjusted\nP-value", 
                        low = "orange", high = "blue") + theme_bw()
