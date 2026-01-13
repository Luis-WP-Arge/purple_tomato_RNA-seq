setwd("/media/willian/HDCaddy/Work/Chaiane/New-pipeline-Star/Annot-4.1/")

#BiocManager::install("DESeq2")
library("DESeq2")

#Gerar a seguinte tabela DEGs-Flesh1-vs-FleshNC.xlsx
count.gene<-read.csv("Flesh-Peel-4.1.featureCounts", header = TRUE, sep = "\t")[,c(1,6:30)]
count.gene.specific<-read.csv("Interest-genes/Anthocyanin-pathway/Solyc11g066580-featureCounts.txt", header = TRUE, sep = "\t")[,c(1,6:30)]

count.gene.1<-rbind.data.frame(count.gene, count.gene.specific)
count.gene<-count.gene.1

colnames(count.gene)

count.gene[grep("Solyc05g008250", count.gene$Geneid),]


colnames(count.gene)<-c("Geneid", "Length",
                        "Flesh1_1", "Flesh1_2", "Flesh1_3",
                        "Flesh2_1", "Flesh2_2", "Flesh2_3",
                        "Flesh5_1", "Flesh5_2", "Flesh5_3",
                        "Flesh_NC1", "Flesh_NC2", "Flesh_NC3",
                        "Peel_1_1", "Peel_1_2", "Peel_1_3",
                        "Peel_2_1", "Peel_2_2", "Peel_2_3",
                        "Peel_5_1", "Peel_5_2", "Peel_5_3",
                        "Peel_NC1", "Peel_NC2", "Peel_NC3")

count.gene$Geneid<-gsub("gene\\:", "", count.gene$Geneid)
#count.gene$Geneid<-gsub("mRNA\\:", "", count.gene$Geneid)
rownames(count.gene)<-count.gene$Geneid
gene.length<-count.gene[,c(1,2)]
rownames(gene.length)<-gene.length$Geneid
gene.length[,1]<-NULL
count.gene<-count.gene[,3:26]
count.gene<-as.matrix(count.gene)

condition <- factor(c(rep("Flesh1",3), rep("Flesh2",3), rep("Flesh5",3), rep("Flesh_NC",3),
                      rep("Peel_1",3), rep("Peel_2",3), rep("Peel_5",3),
                      rep("Peel_NC",3)))
coldata <- data.frame(row.names=colnames(count.gene), condition)

dds <- DESeqDataSetFromMatrix(countData=count.gene, colData=coldata, design=~condition)
dds

dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 4 #normalized counts >=5, and filter out genes where there are less than 3 samples


###################### normalização por TPM
source("~/Dropbox/UFRJ/Scripts/tpm.R")
tpm.gene<-tpm(count.gene, data.frame(gene.length))

tpm.gene.2<-tpm.gene
tpm.gene.2[tpm.gene.2<1]<-0
tpm.gene.2$sum<-rowSums(tpm.gene.2[,1:24]!=0)
tpm.gene.3<-tpm.gene.2[which(tpm.gene.2$sum>=3),1:24]
idx<-rownames(tpm.gene.3)

###################### 

dds <- dds[idx,]
dds <- DESeq(dds)

#png("qc-dispersions-1.png", 1000, 1000, pointsize=20)
pdf("qc-dispersions.pdf", width = 8, height = 8)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
sampleDists <- as.matrix(dist(t(assay(rld))))

#install.packages("gplots")
library(gplots)

#png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
pdf("qc-heatmap-samples.pdf", width = 8, height = 8)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

#png("qc-pca.png", 1000, 1000, pointsize=20)
pdf("qc-pca.pdf", width = 8, height = 8)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

#install.packages("remotes")
#remotes::install_github("twbattaglia/btools")
#BiocManager::install("phyloseq")
#BiocManager::install("vctrs", force = TRUE)
library("vctrs")
library(btools)

options(rgl.useNULL = TRUE)

pdf("3d-plot.pdf")
plotPCA3D(rld, intgroup = "condition", ntop = 500,
          returnData = FALSE)
dev.off()

or
#install.packages("webshot2")
#library("webshot2")
rgl::snapshot3d(filename = "3d.png", width = 900, height = 900)

resultsNames(dds)
Flesh1.vs.FleshNC<-results(dds, contrast = c("condition", "Flesh1", "Flesh_NC"))
Flesh2.vs.FleshNC<-results(dds, contrast = c("condition", "Flesh2", "Flesh_NC"))
Flesh5.vs.FleshNC<-results(dds, contrast = c("condition", "Flesh5", "Flesh_NC"))
Flesh2.vs.Flesh1<-results(dds, contrast = c("condition", "Flesh2", "Flesh1"))
Flesh5.vs.Flesh2<-results(dds, contrast = c("condition", "Flesh5", "Flesh2"))
PeelNC.vs.FleshNC<-results(dds, contrast = c("condition", "Peel_NC", "Flesh_NC"))
Peel1.vs.PeelNC<-results(dds, contrast = c("condition", "Peel_1", "Peel_NC"))
Peel2.vs.PeelNC<-results(dds, contrast = c("condition", "Peel_2", "Peel_NC"))
Peel5.vs.PeelNC<-results(dds, contrast = c("condition", "Peel_5", "Peel_NC"))
Peel2.vs.Peel1<-results(dds, contrast = c("condition", "Peel_2", "Peel_1"))
Peel5.vs.Peel2<-results(dds, contrast = c("condition", "Peel_5", "Peel_2"))

#Flesh1.vs.FleshNC <- results(dds, name="condition_Flesh1_vs_Fles_NC")
#Flesh2.vs.FleshNC <- results(dds, name="condition_Flesh2_vs_Fles_NC")
#Flesh5.vs.FleshNC <- results(dds, name="condition_Flesh5_vs_Fles_NC")
#PeelNC.vs.FleshNC <- results(dds, name="condition_Peel_NC_vs_Fles_NC")

table(Flesh1.vs.FleshNC$padj<=0.05)
table(Flesh2.vs.FleshNC$padj<=0.05)
table(Flesh5.vs.FleshNC$padj<=0.05)
table(Flesh2.vs.Flesh1$padj<=0.05)
table(Flesh5.vs.Flesh2$padj<=0.05)

table(PeelNC.vs.FleshNC$padj<=0.05)
table(Peel1.vs.PeelNC$padj<=0.05)
table(Peel2.vs.PeelNC$padj<=0.05)
table(Peel5.vs.PeelNC$padj<=0.05)
table(Peel2.vs.Peel1$padj<=0.05)
table(Peel5.vs.Peel2$padj<=0.05)

Gene.function<-read.csv("ITAG4.1_descriptions-tab-2-with-dropped-genes.txt", header= FALSE, sep = "\t")
colnames(Gene.function)<-c("Gene.ID", "Function")

count.gene.2<-data.frame(count.gene)
#count.gene.2$Gene<-gsub("\\.\\d", "", row.names(count.gene.2))
count.gene.2$Gene.ID<-row.names(count.gene.2)
count.gene.2$Gene.function<-Gene.function[match(count.gene.2$Gene.ID, Gene.function$Gene.ID), 2]
count.gene.2<-count.gene.2[,c(25,26,1:24)]

count.gene.2$TPM.Flesh1.1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh1_1"]
count.gene.2$TPM.Flesh1.2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh1_2"]
count.gene.2$TPM.Flesh1.3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh1_3"]
count.gene.2$TPM.Flesh2.1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
count.gene.2$TPM.Flesh2.2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
count.gene.2$TPM.Flesh2.3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
count.gene.2$TPM.Flesh5.1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh5_1"]
count.gene.2$TPM.Flesh5.2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh5_2"]
count.gene.2$TPM.Flesh5.3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh5_3"]
count.gene.2$TPM.Flesh.NC1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
count.gene.2$TPM.Flesh.NC2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
count.gene.2$TPM.Flesh.NC3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]

count.gene.2$TPM.Peel_1_1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_1_1"]
count.gene.2$TPM.Peel_1_2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_1_2"]
count.gene.2$TPM.Peel_1_3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_1_3"]
count.gene.2$TPM.Peel_2_1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
count.gene.2$TPM.Peel_2_2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
count.gene.2$TPM.Peel_2_3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
count.gene.2$TPM.Peel_5_1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_5_1"]
count.gene.2$TPM.Peel_5_2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_5_2"]
count.gene.2$TPM.Peel_5_3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_5_3"]
count.gene.2$TPM.Peel.NC1<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
count.gene.2$TPM.Peel.NC2<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
count.gene.2$TPM.Peel.NC3<-tpm.gene[match(count.gene.2$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]

count.gene.2$TPM.Flesh1.mean<-rowMeans(count.gene.2[,c(27:29)])
count.gene.2$TPM.Flesh1.sd<-rowSds(as.matrix(count.gene.2[,c(27:29)]))
count.gene.2$TPM.Flesh2.mean<-rowMeans(count.gene.2[,c(30:32)])
count.gene.2$TPM.Flesh2.sd<-rowSds(as.matrix(count.gene.2[,c(30:32)]))
count.gene.2$TPM.Flesh5.mean<-rowMeans(count.gene.2[,c(33:35)])
count.gene.2$TPM.Flesh5.sd<-rowSds(as.matrix(count.gene.2[,c(33:35)]))
count.gene.2$TPM.Flesh.NC.mean<-rowMeans(count.gene.2[,c(36:38)])
count.gene.2$TPM.Flesh.NC.sd<-rowSds(as.matrix(count.gene.2[,c(36:38)]))

count.gene.2$TPM.Peel1.mean<-rowMeans(count.gene.2[,c(39:41)])
count.gene.2$TPM.Peel1.sd<-rowSds(as.matrix(count.gene.2[,c(39:41)]))
count.gene.2$TPM.Peel2.mean<-rowMeans(count.gene.2[,c(42:44)])
count.gene.2$TPM.Peel2.sd<-rowSds(as.matrix(count.gene.2[,c(42:44)]))
count.gene.2$TPM.Peel5.mean<-rowMeans(count.gene.2[,c(45:47)])
count.gene.2$TPM.Peel5.sd<-rowSds(as.matrix(count.gene.2[,c(45:47)]))
count.gene.2$TPM.Peel.NC.mean<-rowMeans(count.gene.2[,c(48:50)])
count.gene.2$TPM.Peel.NC.sd<-rowSds(as.matrix(count.gene.2[,c(48:50)]))

#install.packages("xlsx")
# library("xlsx")
# write.xlsx(count.gene.2[,c(1:2, 3:26)], "Gene-expression.xlsx", sheetName = "Count - raw", row.names = FALSE)
# write.xlsx(count.gene.2[,c(1:2, 27:50)], "Gene-expression.xlsx", sheetName = "TPM", append = TRUE, row.names = FALSE)
# write.xlsx(count.gene.2[,c(1:2, 51:66)], "Gene-expression.xlsx", sheetName = "Mean + SD", append = TRUE, row.names = FALSE)

library("writexl")
write_xlsx(count.gene.2[,c(1:2, 3:26)], "Gene-expression-raw-counts.xlsx")
write_xlsx(count.gene.2[,c(1:2, 27:50)], "Gene-expression-TPM.xlsx")
write_xlsx(count.gene.2[,c(1:2, 51:66)], "Gene-expression-TPM-mean-sd.xlsx")

write.table(count.gene.2[,c(1:2, 3:26)], "Gene-expression-raw-counts.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(count.gene.2[,c(1:2, 27:50)], "Gene-expression-TPM.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(count.gene.2[,c(1:2, 51:66)], "Gene-expression-TPM-mean-sd.txt", quote = FALSE, row.names = FALSE, sep = "\t")

##################################################################################

Flesh1.vs.FleshNC.1<-data.frame(Flesh1.vs.FleshNC[which(Flesh1.vs.FleshNC$padj<=0.05),])
Flesh1.vs.FleshNC.1$Gene.function<-Gene.function[match(row.names(Flesh1.vs.FleshNC.1), Gene.function$Gene.ID), "Function"]
Flesh1.vs.FleshNC.1$Gene.ID<-rownames(Flesh1.vs.FleshNC.1)
Flesh1.vs.FleshNC.1$TPM.Flesh1.1<-tpm.gene[match(Flesh1.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh1_1"]
Flesh1.vs.FleshNC.1$TPM.Flesh1.2<-tpm.gene[match(Flesh1.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh1_2"]
Flesh1.vs.FleshNC.1$TPM.Flesh1.3<-tpm.gene[match(Flesh1.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh1_3"]
Flesh1.vs.FleshNC.1$TPM.Flesh1.mean<-rowMeans(Flesh1.vs.FleshNC.1[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")])
Flesh1.vs.FleshNC.1$TPM.Flesh1.sd<-rowSds(as.matrix(Flesh1.vs.FleshNC.1[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")]))
Flesh1.vs.FleshNC.1$TPM.Flesh.NC1<-tpm.gene[match(Flesh1.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
Flesh1.vs.FleshNC.1$TPM.Flesh.NC2<-tpm.gene[match(Flesh1.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
Flesh1.vs.FleshNC.1$TPM.Flesh.NC3<-tpm.gene[match(Flesh1.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
Flesh1.vs.FleshNC.1$TPM.Flesh.NC.mean<-rowMeans(Flesh1.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
Flesh1.vs.FleshNC.1$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh1.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))
Flesh1.vs.FleshNC.1.logC1<-Flesh1.vs.FleshNC.1[which(Flesh1.vs.FleshNC.1$log2FoldChange>=1 | Flesh1.vs.FleshNC.1$log2FoldChange<=-1),]

Flesh2.vs.FleshNC.1<-data.frame(Flesh2.vs.FleshNC[which(Flesh2.vs.FleshNC$padj<=0.05),])
Flesh2.vs.FleshNC.1$Gene.function<-Gene.function[match(row.names(Flesh2.vs.FleshNC.1), Gene.function$Gene.ID), "Function"]
Flesh2.vs.FleshNC.1$Gene.ID<-rownames(Flesh2.vs.FleshNC.1)
Flesh2.vs.FleshNC.1$TPM.Flesh2.1<-tpm.gene[match(Flesh2.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
Flesh2.vs.FleshNC.1$TPM.Flesh2.2<-tpm.gene[match(Flesh2.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
Flesh2.vs.FleshNC.1$TPM.Flesh2.3<-tpm.gene[match(Flesh2.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
Flesh2.vs.FleshNC.1$TPM.Flesh2.mean<-rowMeans(Flesh2.vs.FleshNC.1[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")])
Flesh2.vs.FleshNC.1$TPM.Flesh2.sd<-rowSds(as.matrix(Flesh2.vs.FleshNC.1[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")]))
Flesh2.vs.FleshNC.1$TPM.Flesh.NC1<-tpm.gene[match(Flesh2.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
Flesh2.vs.FleshNC.1$TPM.Flesh.NC2<-tpm.gene[match(Flesh2.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
Flesh2.vs.FleshNC.1$TPM.Flesh.NC3<-tpm.gene[match(Flesh2.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
Flesh2.vs.FleshNC.1$TPM.Flesh.NC.mean<-rowMeans(Flesh2.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
Flesh2.vs.FleshNC.1$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh2.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))
Flesh2.vs.FleshNC.1.logC1<-Flesh2.vs.FleshNC.1[which(Flesh2.vs.FleshNC.1$log2FoldChange>=1 | Flesh2.vs.FleshNC.1$log2FoldChange<=-1),]

Flesh5.vs.FleshNC.1<-data.frame(Flesh5.vs.FleshNC[which(Flesh5.vs.FleshNC$padj<=0.05),])
Flesh5.vs.FleshNC.1$Gene.function<-Gene.function[match(row.names(Flesh5.vs.FleshNC.1), Gene.function$Gene.ID), "Function"]
Flesh5.vs.FleshNC.1$Gene.ID<-rownames(Flesh5.vs.FleshNC.1)
Flesh5.vs.FleshNC.1$TPM.Flesh5.1<-tpm.gene[match(Flesh5.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh5_1"]
Flesh5.vs.FleshNC.1$TPM.Flesh5.2<-tpm.gene[match(Flesh5.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh5_2"]
Flesh5.vs.FleshNC.1$TPM.Flesh5.3<-tpm.gene[match(Flesh5.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh5_3"]
Flesh5.vs.FleshNC.1$TPM.Flesh5.mean<-rowMeans(Flesh5.vs.FleshNC.1[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")])
Flesh5.vs.FleshNC.1$TPM.Flesh5.sd<-rowSds(as.matrix(Flesh5.vs.FleshNC.1[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")]))
Flesh5.vs.FleshNC.1$TPM.Flesh.NC1<-tpm.gene[match(Flesh5.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
Flesh5.vs.FleshNC.1$TPM.Flesh.NC2<-tpm.gene[match(Flesh5.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
Flesh5.vs.FleshNC.1$TPM.Flesh.NC3<-tpm.gene[match(Flesh5.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
Flesh5.vs.FleshNC.1$TPM.Flesh.NC.mean<-rowMeans(Flesh5.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
Flesh5.vs.FleshNC.1$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh5.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))
Flesh5.vs.FleshNC.1.logC1<-Flesh5.vs.FleshNC.1[which(Flesh5.vs.FleshNC.1$log2FoldChange>=1 | Flesh5.vs.FleshNC.1$log2FoldChange<=-1),]

Flesh2.vs.Flesh1<-data.frame(Flesh2.vs.Flesh1[which(Flesh2.vs.Flesh1$padj<=0.05),])
Flesh2.vs.Flesh1$Gene.function<-Gene.function[match(row.names(Flesh2.vs.Flesh1), Gene.function$Gene.ID), "Function"]
Flesh2.vs.Flesh1$Gene.ID<-rownames(Flesh2.vs.Flesh1)
Flesh2.vs.Flesh1$TPM.Flesh2.1<-tpm.gene[match(Flesh2.vs.Flesh1$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
Flesh2.vs.Flesh1$TPM.Flesh2.2<-tpm.gene[match(Flesh2.vs.Flesh1$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
Flesh2.vs.Flesh1$TPM.Flesh2.3<-tpm.gene[match(Flesh2.vs.Flesh1$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
Flesh2.vs.Flesh1$TPM.Flesh2.mean<-rowMeans(Flesh2.vs.Flesh1[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")])
Flesh2.vs.Flesh1$TPM.Flesh2.sd<-rowSds(as.matrix(Flesh2.vs.Flesh1[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")]))
Flesh2.vs.Flesh1$TPM.Flesh1.1<-tpm.gene[match(Flesh2.vs.Flesh1$Gene.ID, rownames(tpm.gene)), "Flesh1_1"]
Flesh2.vs.Flesh1$TPM.Flesh1.2<-tpm.gene[match(Flesh2.vs.Flesh1$Gene.ID, rownames(tpm.gene)), "Flesh1_2"]
Flesh2.vs.Flesh1$TPM.Flesh1.3<-tpm.gene[match(Flesh2.vs.Flesh1$Gene.ID, rownames(tpm.gene)), "Flesh1_3"]
Flesh2.vs.Flesh1$TPM.Flesh.NC.mean<-rowMeans(Flesh2.vs.Flesh1[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")])
Flesh2.vs.Flesh1$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh2.vs.Flesh1[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")]))
Flesh2.vs.Flesh1.logC1<-Flesh2.vs.Flesh1[which(Flesh2.vs.Flesh1$log2FoldChange>=1 | Flesh2.vs.Flesh1$log2FoldChange<=-1),]

Flesh5.vs.Flesh2<-data.frame(Flesh5.vs.Flesh2[which(Flesh5.vs.Flesh2$padj<=0.05),])
Flesh5.vs.Flesh2$Gene.function<-Gene.function[match(row.names(Flesh5.vs.Flesh2), Gene.function$Gene.ID), "Function"]
Flesh5.vs.Flesh2$Gene.ID<-rownames(Flesh5.vs.Flesh2)
Flesh5.vs.Flesh2$TPM.Flesh5.1<-tpm.gene[match(Flesh5.vs.Flesh2$Gene.ID, rownames(tpm.gene)), "Flesh5_1"]
Flesh5.vs.Flesh2$TPM.Flesh5.2<-tpm.gene[match(Flesh5.vs.Flesh2$Gene.ID, rownames(tpm.gene)), "Flesh5_2"]
Flesh5.vs.Flesh2$TPM.Flesh5.3<-tpm.gene[match(Flesh5.vs.Flesh2$Gene.ID, rownames(tpm.gene)), "Flesh5_3"]
Flesh5.vs.Flesh2$TPM.Flesh5.mean<-rowMeans(Flesh5.vs.Flesh2[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")])
Flesh5.vs.Flesh2$TPM.Flesh5.sd<-rowSds(as.matrix(Flesh5.vs.Flesh2[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")]))
Flesh5.vs.Flesh2$TPM.Flesh2.1<-tpm.gene[match(Flesh5.vs.Flesh2$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
Flesh5.vs.Flesh2$TPM.Flesh2.2<-tpm.gene[match(Flesh5.vs.Flesh2$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
Flesh5.vs.Flesh2$TPM.Flesh2.3<-tpm.gene[match(Flesh5.vs.Flesh2$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
Flesh5.vs.Flesh2$TPM.Flesh.NC.mean<-rowMeans(Flesh5.vs.Flesh2[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")])
Flesh5.vs.Flesh2$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh5.vs.Flesh2[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")]))
Flesh5.vs.Flesh2.logC1<-Flesh5.vs.Flesh2[which(Flesh5.vs.Flesh2$log2FoldChange>=1 | Flesh5.vs.Flesh2$log2FoldChange<=-1),]


PeelNC.vs.FleshNC.1<-data.frame(PeelNC.vs.FleshNC[which(PeelNC.vs.FleshNC$padj<=0.05),])
PeelNC.vs.FleshNC.1$Gene.function<-Gene.function[match(row.names(PeelNC.vs.FleshNC.1), Gene.function$Gene.ID), "Function"]
PeelNC.vs.FleshNC.1$Gene.ID<-rownames(PeelNC.vs.FleshNC.1)
PeelNC.vs.FleshNC.1$TPM.Flesh.NC1<-tpm.gene[match(PeelNC.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
PeelNC.vs.FleshNC.1$TPM.Flesh.NC2<-tpm.gene[match(PeelNC.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
PeelNC.vs.FleshNC.1$TPM.Flesh.NC3<-tpm.gene[match(PeelNC.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
PeelNC.vs.FleshNC.1$TPM.Flesh.NC.mean<-rowMeans(PeelNC.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
PeelNC.vs.FleshNC.1$TPM.Flesh.NC.sd<-rowSds(as.matrix(PeelNC.vs.FleshNC.1[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))
PeelNC.vs.FleshNC.1$TPM.Peel.NC1<-tpm.gene[match(PeelNC.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
PeelNC.vs.FleshNC.1$TPM.Peel.NC2<-tpm.gene[match(PeelNC.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
PeelNC.vs.FleshNC.1$TPM.Peel.NC3<-tpm.gene[match(PeelNC.vs.FleshNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
PeelNC.vs.FleshNC.1$TPM.Peel.NC.mean<-rowMeans(PeelNC.vs.FleshNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
PeelNC.vs.FleshNC.1$TPM.Peel.NC.sd<-rowSds(as.matrix(PeelNC.vs.FleshNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))
PeelNC.vs.FleshNC.1.logC1<-PeelNC.vs.FleshNC.1[which(PeelNC.vs.FleshNC.1$log2FoldChange>=1 | PeelNC.vs.FleshNC.1$log2FoldChange<=-1),]

Peel1.vs.PeelNC.1<-data.frame(Peel1.vs.PeelNC[which(Peel1.vs.PeelNC$padj<=0.05),])
Peel1.vs.PeelNC.1$Gene.function<-Gene.function[match(row.names(Peel1.vs.PeelNC.1), Gene.function$Gene.ID), "Function"]
Peel1.vs.PeelNC.1$Gene.ID<-rownames(Peel1.vs.PeelNC.1)
Peel1.vs.PeelNC.1$TPM.Peel_1_1<-tpm.gene[match(Peel1.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_1_1"]
Peel1.vs.PeelNC.1$TPM.Peel_1_2<-tpm.gene[match(Peel1.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_1_2"]
Peel1.vs.PeelNC.1$TPM.Peel_1_3<-tpm.gene[match(Peel1.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_1_3"]
Peel1.vs.PeelNC.1$TPM.Peel_1.mean<-rowMeans(Peel1.vs.PeelNC.1[,c("TPM.Peel_1_1", "TPM.Peel_1_2", "TPM.Peel_1_3")])
Peel1.vs.PeelNC.1$TPM.Peel_1.sd<-rowSds(as.matrix(Peel1.vs.PeelNC.1[,c("TPM.Peel_1_1", "TPM.Peel_1_2", "TPM.Peel_1_3")]))
Peel1.vs.PeelNC.1$TPM.Peel.NC1<-tpm.gene[match(Peel1.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
Peel1.vs.PeelNC.1$TPM.Peel.NC2<-tpm.gene[match(Peel1.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
Peel1.vs.PeelNC.1$TPM.Peel.NC3<-tpm.gene[match(Peel1.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
Peel1.vs.PeelNC.1$TPM.Peel.NC.mean<-rowMeans(Peel1.vs.PeelNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
Peel1.vs.PeelNC.1$TPM.Peel.NC.sd<-rowSds(as.matrix(Peel1.vs.PeelNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))
Peel1.vs.PeelNC.1.logC1<-Peel1.vs.PeelNC.1[which(Peel1.vs.PeelNC.1$log2FoldChange>=1 | Peel1.vs.PeelNC.1$log2FoldChange<=-1),]

Peel2.vs.PeelNC.1<-data.frame(Peel2.vs.PeelNC[which(Peel2.vs.PeelNC$padj<=0.05),])
Peel2.vs.PeelNC.1$Gene.function<-Gene.function[match(row.names(Peel2.vs.PeelNC.1), Gene.function$Gene.ID), "Function"]
Peel2.vs.PeelNC.1$Gene.ID<-rownames(Peel2.vs.PeelNC.1)
Peel2.vs.PeelNC.1$TPM.Peel_2_1<-tpm.gene[match(Peel2.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
Peel2.vs.PeelNC.1$TPM.Peel_2_2<-tpm.gene[match(Peel2.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
Peel2.vs.PeelNC.1$TPM.Peel_2_3<-tpm.gene[match(Peel2.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
Peel2.vs.PeelNC.1$TPM.Peel_2.mean<-rowMeans(Peel2.vs.PeelNC.1[,c("TPM.Peel_2_1", "TPM.Peel_2_2", "TPM.Peel_2_3")])
Peel2.vs.PeelNC.1$TPM.Peel_2.sd<-rowSds(as.matrix(Peel2.vs.PeelNC.1[,c("TPM.Peel_2_1", "TPM.Peel_2_2", "TPM.Peel_2_3")]))
Peel2.vs.PeelNC.1$TPM.Peel.NC1<-tpm.gene[match(Peel2.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
Peel2.vs.PeelNC.1$TPM.Peel.NC2<-tpm.gene[match(Peel2.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
Peel2.vs.PeelNC.1$TPM.Peel.NC3<-tpm.gene[match(Peel2.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
Peel2.vs.PeelNC.1$TPM.Peel.NC.mean<-rowMeans(Peel2.vs.PeelNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
Peel2.vs.PeelNC.1$TPM.Peel.NC.sd<-rowSds(as.matrix(Peel2.vs.PeelNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))
Peel2.vs.PeelNC.1.logC1<-Peel2.vs.PeelNC.1[which(Peel2.vs.PeelNC.1$log2FoldChange>=1 | Peel2.vs.PeelNC.1$log2FoldChange<=-1),]

Peel5.vs.PeelNC.1<-data.frame(Peel5.vs.PeelNC[which(Peel5.vs.PeelNC$padj<=0.05),])
Peel5.vs.PeelNC.1$Gene.function<-Gene.function[match(row.names(Peel5.vs.PeelNC.1), Gene.function$Gene.ID), "Function"]
Peel5.vs.PeelNC.1$Gene.ID<-rownames(Peel5.vs.PeelNC.1)
Peel5.vs.PeelNC.1$TPM.Peel_5_1<-tpm.gene[match(Peel5.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_5_1"]
Peel5.vs.PeelNC.1$TPM.Peel_5_2<-tpm.gene[match(Peel5.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_5_2"]
Peel5.vs.PeelNC.1$TPM.Peel_5_3<-tpm.gene[match(Peel5.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_5_3"]
Peel5.vs.PeelNC.1$TPM.Peel_5.mean<-rowMeans(Peel5.vs.PeelNC.1[,c("TPM.Peel_5_1", "TPM.Peel_5_2", "TPM.Peel_5_3")])
Peel5.vs.PeelNC.1$TPM.Peel_5.sd<-rowSds(as.matrix(Peel5.vs.PeelNC.1[,c("TPM.Peel_5_1", "TPM.Peel_5_2", "TPM.Peel_5_3")]))
Peel5.vs.PeelNC.1$TPM.Peel.NC1<-tpm.gene[match(Peel5.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
Peel5.vs.PeelNC.1$TPM.Peel.NC2<-tpm.gene[match(Peel5.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
Peel5.vs.PeelNC.1$TPM.Peel.NC3<-tpm.gene[match(Peel5.vs.PeelNC.1$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
Peel5.vs.PeelNC.1$TPM.Peel.NC.mean<-rowMeans(Peel5.vs.PeelNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
Peel5.vs.PeelNC.1$TPM.Peel.NC.sd<-rowSds(as.matrix(Peel5.vs.PeelNC.1[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))
Peel5.vs.PeelNC.1.logC1<-Peel5.vs.PeelNC.1[which(Peel5.vs.PeelNC.1$log2FoldChange>=1 | Peel5.vs.PeelNC.1$log2FoldChange<=-1),]

Peel2.vs.Peel1<-data.frame(Peel2.vs.Peel1[which(Peel2.vs.Peel1$padj<=0.05),])
Peel2.vs.Peel1$Gene.function<-Gene.function[match(row.names(Peel2.vs.Peel1), Gene.function$Gene.ID), "Function"]
Peel2.vs.Peel1$Gene.ID<-rownames(Peel2.vs.Peel1)
Peel2.vs.Peel1$TPM.Peel2.1<-tpm.gene[match(Peel2.vs.Peel1$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
Peel2.vs.Peel1$TPM.Peel2.2<-tpm.gene[match(Peel2.vs.Peel1$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
Peel2.vs.Peel1$TPM.Peel2.3<-tpm.gene[match(Peel2.vs.Peel1$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
Peel2.vs.Peel1$TPM.Peel2.mean<-rowMeans(Peel2.vs.Peel1[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")])
Peel2.vs.Peel1$TPM.Peel2.sd<-rowSds(as.matrix(Peel2.vs.Peel1[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")]))
Peel2.vs.Peel1$TPM.Peel1.1<-tpm.gene[match(Peel2.vs.Peel1$Gene.ID, rownames(tpm.gene)), "Peel_1_1"]
Peel2.vs.Peel1$TPM.Peel1.2<-tpm.gene[match(Peel2.vs.Peel1$Gene.ID, rownames(tpm.gene)), "Peel_1_2"]
Peel2.vs.Peel1$TPM.Peel1.3<-tpm.gene[match(Peel2.vs.Peel1$Gene.ID, rownames(tpm.gene)), "Peel_1_3"]
Peel2.vs.Peel1$TPM.Peel.NC.mean<-rowMeans(Peel2.vs.Peel1[,c("TPM.Peel1.1", "TPM.Peel1.2", "TPM.Peel1.3")])
Peel2.vs.Peel1$TPM.Peel1.sd<-rowSds(as.matrix(Peel2.vs.Peel1[,c("TPM.Peel1.1", "TPM.Peel1.2", "TPM.Peel1.3")]))
Peel2.vs.Peel1.logC1<-Peel2.vs.Peel1[which(Peel2.vs.Peel1$log2FoldChange>=1 | Peel2.vs.Peel1$log2FoldChange<=-1),]

Peel5.vs.Peel2<-data.frame(Peel5.vs.Peel2[which(Peel5.vs.Peel2$padj<=0.05),])
Peel5.vs.Peel2$Gene.function<-Gene.function[match(row.names(Peel5.vs.Peel2), Gene.function$Gene.ID), "Function"]
Peel5.vs.Peel2$Gene.ID<-rownames(Peel5.vs.Peel2)
Peel5.vs.Peel2$TPM.Peel5.1<-tpm.gene[match(Peel5.vs.Peel2$Gene.ID, rownames(tpm.gene)), "Peel_5_1"]
Peel5.vs.Peel2$TPM.Peel5.2<-tpm.gene[match(Peel5.vs.Peel2$Gene.ID, rownames(tpm.gene)), "Peel_5_2"]
Peel5.vs.Peel2$TPM.Peel5.3<-tpm.gene[match(Peel5.vs.Peel2$Gene.ID, rownames(tpm.gene)), "Peel_5_3"]
Peel5.vs.Peel2$TPM.Peel5.mean<-rowMeans(Peel5.vs.Peel2[,c("TPM.Peel5.1", "TPM.Peel5.2", "TPM.Peel5.3")])
Peel5.vs.Peel2$TPM.Peel5.sd<-rowSds(as.matrix(Peel5.vs.Peel2[,c("TPM.Peel5.1", "TPM.Peel5.2", "TPM.Peel5.3")]))
Peel5.vs.Peel2$TPM.Peel2.1<-tpm.gene[match(Peel5.vs.Peel2$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
Peel5.vs.Peel2$TPM.Peel2.2<-tpm.gene[match(Peel5.vs.Peel2$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
Peel5.vs.Peel2$TPM.Peel2.3<-tpm.gene[match(Peel5.vs.Peel2$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
Peel5.vs.Peel2$TPM.Peel2.mean<-rowMeans(Peel5.vs.Peel2[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")])
Peel5.vs.Peel2$TPM.Peel2.sd<-rowSds(as.matrix(Peel5.vs.Peel2[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")]))
Peel5.vs.Peel2.logC1<-Peel5.vs.Peel2[which(Peel5.vs.Peel2$log2FoldChange>=1 | Peel5.vs.Peel2$log2FoldChange<=-1),]

#install.packages("writexl")
library("writexl")
write_xlsx(data.frame(Flesh1.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh1-vs-FleshNC.xlsx")
write_xlsx(data.frame(Flesh2.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh2-vs-FleshNC.xlsx")
write_xlsx(data.frame(Flesh5.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh5-vs-FleshNC.xlsx")
write_xlsx(data.frame(Flesh2.vs.Flesh1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh2-vs-Flesh1.xlsx")
write_xlsx(data.frame(Flesh5.vs.Flesh2[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh5-vs-Flesh2.xlsx")
write_xlsx(data.frame(PeelNC.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-PeelNC-vs-FleshNC.xlsx")
write_xlsx(data.frame(Peel1.vs.PeelNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel1-vs-PeelNC.xlsx")
write_xlsx(data.frame(Peel2.vs.PeelNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel2-vs-PeelNC.xlsx")
write_xlsx(data.frame(Peel5.vs.PeelNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel5-vs-PeelNC.xlsx") 
write_xlsx(data.frame(Peel2.vs.Peel1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel2-vs-Peel1.xlsx")
write_xlsx(data.frame(Peel5.vs.Peel2[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel5-vs-Peel2.xlsx")

write_xlsx(data.frame(Flesh1.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-Flesh1-vs-FleshNC-logFC1.xlsx")
write_xlsx(data.frame(Flesh2.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-Flesh2-vs-FleshNC-logFC1.xlsx")
write_xlsx(data.frame(Flesh5.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-Flesh5-vs-FleshNC-logFC1.xlsx")
write_xlsx(data.frame(Flesh2.vs.Flesh1.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh2-vs-Flesh1-logFC1.xlsx")
write_xlsx(data.frame(Flesh5.vs.Flesh2.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Flesh5-vs-Flesh2-logFC1.xlsx")
write_xlsx(data.frame(PeelNC.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-PeelNC-vs-FleshNC-logFC1.xlsx")
write_xlsx(data.frame(Peel1.vs.PeelNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-Peel1-vs-PeelNC-logFC1.xlsx")
write_xlsx(data.frame(Peel2.vs.PeelNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-Peel2-vs-PeelNC-logFC1.xlsx")
write_xlsx(data.frame(Peel5.vs.PeelNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-excel/DEGs-Peel5-vs-PeelNC-logFC1.xlsx") 
write_xlsx(data.frame(Peel2.vs.Peel1.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel2-vs-Peel1-logFC1.xlsx")
write_xlsx(data.frame(Peel5.vs.Peel2.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-excel/DEGs-Peel5-vs-Peel2-logFC1.xlsx")

####################################################
write.table(data.frame(Flesh1.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh1-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh2.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh2-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh5.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh5-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh2.vs.Flesh1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh2-vs-Flesh1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh5.vs.Flesh2[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh5-vs-Flesh2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(PeelNC.vs.FleshNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-PeelNC-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel1.vs.PeelNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel1-vs-PeelNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel2.vs.PeelNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel2-vs-PeelNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel5.vs.PeelNC.1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel5-vs-PeelNC.txt", quote = FALSE, row.names = FALSE, sep = "\t") 
write.table(data.frame(Peel2.vs.Peel1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel2-vs-Peel1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel5.vs.Peel2[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel5-vs-Peel2.txt", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(data.frame(Flesh1.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-Flesh1-vs-FleshNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh2.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-Flesh2-vs-FleshNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh5.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-Flesh5-vs-FleshNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh2.vs.Flesh1.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh2-vs-Flesh1-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh5.vs.Flesh2.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Flesh5-vs-Flesh2-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(PeelNC.vs.FleshNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-PeelNC-vs-FleshNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel1.vs.PeelNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-Peel1-vs-PeelNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel2.vs.PeelNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-Peel2-vs-PeelNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel5.vs.PeelNC.1.logC1[,c(8,7,1:6,9:18)]), "DEGs-scale-1/DEGs-txt/DEGs-Peel5-vs-PeelNC-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t") 
write.table(data.frame(Peel2.vs.Peel1.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel2-vs-Peel1-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel5.vs.Peel2.logC1[,c(8,7,1:6,9:18)]), "DEGs-raw/DEGs-txt/DEGs-Peel5-vs-Peel2-logFC1.txt", quote = FALSE, row.names = FALSE, sep = "\t")

########################################################################################################


Flesh1.vs.FleshNC.1.no.filter<-data.frame(Flesh1.vs.FleshNC)
Flesh1.vs.FleshNC.1.no.filter$Gene.function<-Gene.function[match(row.names(Flesh1.vs.FleshNC.1.no.filter), Gene.function$Gene.ID), "Function"]
Flesh1.vs.FleshNC.1.no.filter$Gene.ID<-rownames(Flesh1.vs.FleshNC.1.no.filter)
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh1.1<-tpm.gene[match(Flesh1.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh1_1"]
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh1.2<-tpm.gene[match(Flesh1.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh1_2"]
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh1.3<-tpm.gene[match(Flesh1.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh1_3"]
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh1.mean<-rowMeans(Flesh1.vs.FleshNC.1.no.filter[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")])
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh1.sd<-rowSds(as.matrix(Flesh1.vs.FleshNC.1.no.filter[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")]))
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh.NC1<-tpm.gene[match(Flesh1.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh.NC2<-tpm.gene[match(Flesh1.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh.NC3<-tpm.gene[match(Flesh1.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh.NC.mean<-rowMeans(Flesh1.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
Flesh1.vs.FleshNC.1.no.filter$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh1.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))


Flesh2.vs.FleshNC.1.no.filter<-data.frame(Flesh2.vs.FleshNC)
Flesh2.vs.FleshNC.1.no.filter$Gene.function<-Gene.function[match(row.names(Flesh2.vs.FleshNC.1.no.filter), Gene.function$Gene.ID), "Function"]
Flesh2.vs.FleshNC.1.no.filter$Gene.ID<-rownames(Flesh2.vs.FleshNC.1.no.filter)
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh2.1<-tpm.gene[match(Flesh2.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh2.2<-tpm.gene[match(Flesh2.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh2.3<-tpm.gene[match(Flesh2.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh2.mean<-rowMeans(Flesh2.vs.FleshNC.1.no.filter[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")])
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh2.sd<-rowSds(as.matrix(Flesh2.vs.FleshNC.1.no.filter[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")]))
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh.NC1<-tpm.gene[match(Flesh2.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh.NC2<-tpm.gene[match(Flesh2.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh.NC3<-tpm.gene[match(Flesh2.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh.NC.mean<-rowMeans(Flesh2.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
Flesh2.vs.FleshNC.1.no.filter$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh2.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))

Flesh5.vs.FleshNC.1.no.filter<-data.frame(Flesh5.vs.FleshNC)
Flesh5.vs.FleshNC.1.no.filter$Gene.function<-Gene.function[match(row.names(Flesh5.vs.FleshNC.1.no.filter), Gene.function$Gene.ID), "Function"]
Flesh5.vs.FleshNC.1.no.filter$Gene.ID<-rownames(Flesh5.vs.FleshNC.1.no.filter)
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh5.1<-tpm.gene[match(Flesh5.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh5_1"]
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh5.2<-tpm.gene[match(Flesh5.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh5_2"]
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh5.3<-tpm.gene[match(Flesh5.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh5_3"]
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh5.mean<-rowMeans(Flesh5.vs.FleshNC.1.no.filter[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")])
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh5.sd<-rowSds(as.matrix(Flesh5.vs.FleshNC.1.no.filter[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")]))
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh.NC1<-tpm.gene[match(Flesh5.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh.NC2<-tpm.gene[match(Flesh5.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh.NC3<-tpm.gene[match(Flesh5.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh.NC.mean<-rowMeans(Flesh5.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
Flesh5.vs.FleshNC.1.no.filter$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh5.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))

Flesh2.vs.Flesh1.no.filter<-data.frame(Flesh2.vs.Flesh1)
Flesh2.vs.Flesh1.no.filter$Gene.function<-Gene.function[match(row.names(Flesh2.vs.Flesh1.no.filter), Gene.function$Gene.ID), "Function"]
Flesh2.vs.Flesh1.no.filter$Gene.ID<-rownames(Flesh2.vs.Flesh1.no.filter)
Flesh2.vs.Flesh1.no.filter$TPM.Flesh2.1<-tpm.gene[match(Flesh2.vs.Flesh1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
Flesh2.vs.Flesh1.no.filter$TPM.Flesh2.2<-tpm.gene[match(Flesh2.vs.Flesh1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
Flesh2.vs.Flesh1.no.filter$TPM.Flesh2.3<-tpm.gene[match(Flesh2.vs.Flesh1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
Flesh2.vs.Flesh1.no.filter$TPM.Flesh2.mean<-rowMeans(Flesh2.vs.Flesh1.no.filter[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")])
Flesh2.vs.Flesh1.no.filter$TPM.Flesh2.sd<-rowSds(as.matrix(Flesh2.vs.Flesh1.no.filter[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")]))
Flesh2.vs.Flesh1.no.filter$TPM.Flesh1.1<-tpm.gene[match(Flesh2.vs.Flesh1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh1_1"]
Flesh2.vs.Flesh1.no.filter$TPM.Flesh1.2<-tpm.gene[match(Flesh2.vs.Flesh1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh1_2"]
Flesh2.vs.Flesh1.no.filter$TPM.Flesh1.3<-tpm.gene[match(Flesh2.vs.Flesh1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh1_3"]
Flesh2.vs.Flesh1.no.filter$TPM.Flesh.NC.mean<-rowMeans(Flesh2.vs.Flesh1.no.filter[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")])
Flesh2.vs.Flesh1.no.filter$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh2.vs.Flesh1.no.filter[,c("TPM.Flesh1.1", "TPM.Flesh1.2", "TPM.Flesh1.3")]))

Flesh5.vs.Flesh2.no.filter<-data.frame(Flesh5.vs.Flesh2)
Flesh5.vs.Flesh2.no.filter$Gene.function<-Gene.function[match(row.names(Flesh5.vs.Flesh2.no.filter), Gene.function$Gene.ID), "Function"]
Flesh5.vs.Flesh2.no.filter$Gene.ID<-rownames(Flesh5.vs.Flesh2.no.filter)
Flesh5.vs.Flesh2.no.filter$TPM.Flesh5.1<-tpm.gene[match(Flesh5.vs.Flesh2.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh5_1"]
Flesh5.vs.Flesh2.no.filter$TPM.Flesh5.2<-tpm.gene[match(Flesh5.vs.Flesh2.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh5_2"]
Flesh5.vs.Flesh2.no.filter$TPM.Flesh5.3<-tpm.gene[match(Flesh5.vs.Flesh2.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh5_3"]
Flesh5.vs.Flesh2.no.filter$TPM.Flesh5.mean<-rowMeans(Flesh5.vs.Flesh2.no.filter[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")])
Flesh5.vs.Flesh2.no.filter$TPM.Flesh5.sd<-rowSds(as.matrix(Flesh5.vs.Flesh2.no.filter[,c("TPM.Flesh5.1", "TPM.Flesh5.2", "TPM.Flesh5.3")]))
Flesh5.vs.Flesh2.no.filter$TPM.Flesh2.1<-tpm.gene[match(Flesh5.vs.Flesh2.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_1"]
Flesh5.vs.Flesh2.no.filter$TPM.Flesh2.2<-tpm.gene[match(Flesh5.vs.Flesh2.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_2"]
Flesh5.vs.Flesh2.no.filter$TPM.Flesh2.3<-tpm.gene[match(Flesh5.vs.Flesh2.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh2_3"]
Flesh5.vs.Flesh2.no.filter$TPM.Flesh.NC.mean<-rowMeans(Flesh5.vs.Flesh2.no.filter[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")])
Flesh5.vs.Flesh2.no.filter$TPM.Flesh.NC.sd<-rowSds(as.matrix(Flesh5.vs.Flesh2.no.filter[,c("TPM.Flesh2.1", "TPM.Flesh2.2", "TPM.Flesh2.3")]))


PeelNC.vs.FleshNC.1.no.filter<-data.frame(PeelNC.vs.FleshNC)
PeelNC.vs.FleshNC.1.no.filter$Gene.function<-Gene.function[match(row.names(PeelNC.vs.FleshNC.1.no.filter), Gene.function$Gene.ID), "Function"]
PeelNC.vs.FleshNC.1.no.filter$Gene.ID<-rownames(PeelNC.vs.FleshNC.1.no.filter)
PeelNC.vs.FleshNC.1.no.filter$TPM.Flesh.NC1<-tpm.gene[match(PeelNC.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC1"]
PeelNC.vs.FleshNC.1.no.filter$TPM.Flesh.NC2<-tpm.gene[match(PeelNC.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC2"]
PeelNC.vs.FleshNC.1.no.filter$TPM.Flesh.NC3<-tpm.gene[match(PeelNC.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Flesh_NC3"]
PeelNC.vs.FleshNC.1.no.filter$TPM.Flesh.NC.mean<-rowMeans(PeelNC.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")])
PeelNC.vs.FleshNC.1.no.filter$TPM.Flesh.NC.sd<-rowSds(as.matrix(PeelNC.vs.FleshNC.1.no.filter[,c("TPM.Flesh.NC1", "TPM.Flesh.NC2", "TPM.Flesh.NC3")]))
PeelNC.vs.FleshNC.1.no.filter$TPM.Peel.NC1<-tpm.gene[match(PeelNC.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
PeelNC.vs.FleshNC.1.no.filter$TPM.Peel.NC2<-tpm.gene[match(PeelNC.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
PeelNC.vs.FleshNC.1.no.filter$TPM.Peel.NC3<-tpm.gene[match(PeelNC.vs.FleshNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
PeelNC.vs.FleshNC.1.no.filter$TPM.Peel.NC.mean<-rowMeans(PeelNC.vs.FleshNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
PeelNC.vs.FleshNC.1.no.filter$TPM.Peel.NC.sd<-rowSds(as.matrix(PeelNC.vs.FleshNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))

Peel1.vs.PeelNC.1.no.filter<-data.frame(Peel1.vs.PeelNC)
Peel1.vs.PeelNC.1.no.filter$Gene.function<-Gene.function[match(row.names(Peel1.vs.PeelNC.1.no.filter), Gene.function$Gene.ID), "Function"]
Peel1.vs.PeelNC.1.no.filter$Gene.ID<-rownames(Peel1.vs.PeelNC.1.no.filter)
Peel1.vs.PeelNC.1.no.filter$TPM.Peel_1_1<-tpm.gene[match(Peel1.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_1_1"]
Peel1.vs.PeelNC.1.no.filter$TPM.Peel_1_2<-tpm.gene[match(Peel1.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_1_2"]
Peel1.vs.PeelNC.1.no.filter$TPM.Peel_1_3<-tpm.gene[match(Peel1.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_1_3"]
Peel1.vs.PeelNC.1.no.filter$TPM.Peel_1.mean<-rowMeans(Peel1.vs.PeelNC.1.no.filter[,c("TPM.Peel_1_1", "TPM.Peel_1_2", "TPM.Peel_1_3")])
Peel1.vs.PeelNC.1.no.filter$TPM.Peel_1.sd<-rowSds(as.matrix(Peel1.vs.PeelNC.1.no.filter[,c("TPM.Peel_1_1", "TPM.Peel_1_2", "TPM.Peel_1_3")]))
Peel1.vs.PeelNC.1.no.filter$TPM.Peel.NC1<-tpm.gene[match(Peel1.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
Peel1.vs.PeelNC.1.no.filter$TPM.Peel.NC2<-tpm.gene[match(Peel1.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
Peel1.vs.PeelNC.1.no.filter$TPM.Peel.NC3<-tpm.gene[match(Peel1.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
Peel1.vs.PeelNC.1.no.filter$TPM.Peel.NC.mean<-rowMeans(Peel1.vs.PeelNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
Peel1.vs.PeelNC.1.no.filter$TPM.Peel.NC.sd<-rowSds(as.matrix(Peel1.vs.PeelNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))

Peel2.vs.PeelNC.1.no.filter<-data.frame(Peel2.vs.PeelNC)
Peel2.vs.PeelNC.1.no.filter$Gene.function<-Gene.function[match(row.names(Peel2.vs.PeelNC.1.no.filter), Gene.function$Gene.ID), "Function"]
Peel2.vs.PeelNC.1.no.filter$Gene.ID<-rownames(Peel2.vs.PeelNC.1.no.filter)
Peel2.vs.PeelNC.1.no.filter$TPM.Peel_2_1<-tpm.gene[match(Peel2.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
Peel2.vs.PeelNC.1.no.filter$TPM.Peel_2_2<-tpm.gene[match(Peel2.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
Peel2.vs.PeelNC.1.no.filter$TPM.Peel_2_3<-tpm.gene[match(Peel2.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
Peel2.vs.PeelNC.1.no.filter$TPM.Peel_2.mean<-rowMeans(Peel2.vs.PeelNC.1.no.filter[,c("TPM.Peel_2_1", "TPM.Peel_2_2", "TPM.Peel_2_3")])
Peel2.vs.PeelNC.1.no.filter$TPM.Peel_2.sd<-rowSds(as.matrix(Peel2.vs.PeelNC.1.no.filter[,c("TPM.Peel_2_1", "TPM.Peel_2_2", "TPM.Peel_2_3")]))
Peel2.vs.PeelNC.1.no.filter$TPM.Peel.NC1<-tpm.gene[match(Peel2.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
Peel2.vs.PeelNC.1.no.filter$TPM.Peel.NC2<-tpm.gene[match(Peel2.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
Peel2.vs.PeelNC.1.no.filter$TPM.Peel.NC3<-tpm.gene[match(Peel2.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
Peel2.vs.PeelNC.1.no.filter$TPM.Peel.NC.mean<-rowMeans(Peel2.vs.PeelNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
Peel2.vs.PeelNC.1.no.filter$TPM.Peel.NC.sd<-rowSds(as.matrix(Peel2.vs.PeelNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))

Peel5.vs.PeelNC.1.no.filter<-data.frame(Peel5.vs.PeelNC)
Peel5.vs.PeelNC.1.no.filter$Gene.function<-Gene.function[match(row.names(Peel5.vs.PeelNC.1.no.filter), Gene.function$Gene.ID), "Function"]
Peel5.vs.PeelNC.1.no.filter$Gene.ID<-rownames(Peel5.vs.PeelNC.1.no.filter)
Peel5.vs.PeelNC.1.no.filter$TPM.Peel_5_1<-tpm.gene[match(Peel5.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_5_1"]
Peel5.vs.PeelNC.1.no.filter$TPM.Peel_5_2<-tpm.gene[match(Peel5.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_5_2"]
Peel5.vs.PeelNC.1.no.filter$TPM.Peel_5_3<-tpm.gene[match(Peel5.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_5_3"]
Peel5.vs.PeelNC.1.no.filter$TPM.Peel_5.mean<-rowMeans(Peel5.vs.PeelNC.1.no.filter[,c("TPM.Peel_5_1", "TPM.Peel_5_2", "TPM.Peel_5_3")])
Peel5.vs.PeelNC.1.no.filter$TPM.Peel_5.sd<-rowSds(as.matrix(Peel5.vs.PeelNC.1.no.filter[,c("TPM.Peel_5_1", "TPM.Peel_5_2", "TPM.Peel_5_3")]))
Peel5.vs.PeelNC.1.no.filter$TPM.Peel.NC1<-tpm.gene[match(Peel5.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC1"]
Peel5.vs.PeelNC.1.no.filter$TPM.Peel.NC2<-tpm.gene[match(Peel5.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC2"]
Peel5.vs.PeelNC.1.no.filter$TPM.Peel.NC3<-tpm.gene[match(Peel5.vs.PeelNC.1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_NC3"]
Peel5.vs.PeelNC.1.no.filter$TPM.Peel.NC.mean<-rowMeans(Peel5.vs.PeelNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")])
Peel5.vs.PeelNC.1.no.filter$TPM.Peel.NC.sd<-rowSds(as.matrix(Peel5.vs.PeelNC.1.no.filter[,c("TPM.Peel.NC1", "TPM.Peel.NC2", "TPM.Peel.NC3")]))

Peel2.vs.Peel1.no.filter<-data.frame(Peel2.vs.Peel1)
Peel2.vs.Peel1.no.filter$Gene.function<-Gene.function[match(row.names(Peel2.vs.Peel1.no.filter), Gene.function$Gene.ID), "Function"]
Peel2.vs.Peel1.no.filter$Gene.ID<-rownames(Peel2.vs.Peel1.no.filter)
Peel2.vs.Peel1.no.filter$TPM.Peel2.1<-tpm.gene[match(Peel2.vs.Peel1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
Peel2.vs.Peel1.no.filter$TPM.Peel2.2<-tpm.gene[match(Peel2.vs.Peel1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
Peel2.vs.Peel1.no.filter$TPM.Peel2.3<-tpm.gene[match(Peel2.vs.Peel1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
Peel2.vs.Peel1.no.filter$TPM.Peel2.mean<-rowMeans(Peel2.vs.Peel1.no.filter[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")])
Peel2.vs.Peel1.no.filter$TPM.Peel2.sd<-rowSds(as.matrix(Peel2.vs.Peel1.no.filter[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")]))
Peel2.vs.Peel1.no.filter$TPM.Peel1.1<-tpm.gene[match(Peel2.vs.Peel1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_1_1"]
Peel2.vs.Peel1.no.filter$TPM.Peel1.2<-tpm.gene[match(Peel2.vs.Peel1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_1_2"]
Peel2.vs.Peel1.no.filter$TPM.Peel1.3<-tpm.gene[match(Peel2.vs.Peel1.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_1_3"]
Peel2.vs.Peel1.no.filter$TPM.Peel.NC.mean<-rowMeans(Peel2.vs.Peel1.no.filter[,c("TPM.Peel1.1", "TPM.Peel1.2", "TPM.Peel1.3")])
Peel2.vs.Peel1.no.filter$TPM.Peel1.sd<-rowSds(as.matrix(Peel2.vs.Peel1.no.filter[,c("TPM.Peel1.1", "TPM.Peel1.2", "TPM.Peel1.3")]))

Peel5.vs.Peel2.no.filter<-data.frame(Peel5.vs.Peel2)
Peel5.vs.Peel2.no.filter$Gene.function<-Gene.function[match(row.names(Peel5.vs.Peel2.no.filter), Gene.function$Gene.ID), "Function"]
Peel5.vs.Peel2.no.filter$Gene.ID<-rownames(Peel5.vs.Peel2.no.filter)
Peel5.vs.Peel2.no.filter$TPM.Peel5.1<-tpm.gene[match(Peel5.vs.Peel2.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_5_1"]
Peel5.vs.Peel2.no.filter$TPM.Peel5.2<-tpm.gene[match(Peel5.vs.Peel2.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_5_2"]
Peel5.vs.Peel2.no.filter$TPM.Peel5.3<-tpm.gene[match(Peel5.vs.Peel2.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_5_3"]
Peel5.vs.Peel2.no.filter$TPM.Peel5.mean<-rowMeans(Peel5.vs.Peel2.no.filter[,c("TPM.Peel5.1", "TPM.Peel5.2", "TPM.Peel5.3")])
Peel5.vs.Peel2.no.filter$TPM.Peel5.sd<-rowSds(as.matrix(Peel5.vs.Peel2.no.filter[,c("TPM.Peel5.1", "TPM.Peel5.2", "TPM.Peel5.3")]))
Peel5.vs.Peel2.no.filter$TPM.Peel2.1<-tpm.gene[match(Peel5.vs.Peel2.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_1"]
Peel5.vs.Peel2.no.filter$TPM.Peel2.2<-tpm.gene[match(Peel5.vs.Peel2.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_2"]
Peel5.vs.Peel2.no.filter$TPM.Peel2.3<-tpm.gene[match(Peel5.vs.Peel2.no.filter$Gene.ID, rownames(tpm.gene)), "Peel_2_3"]
Peel5.vs.Peel2.no.filter$TPM.Peel2.mean<-rowMeans(Peel5.vs.Peel2.no.filter[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")])
Peel5.vs.Peel2.no.filter$TPM.Peel2.sd<-rowSds(as.matrix(Peel5.vs.Peel2.no.filter[,c("TPM.Peel2.1", "TPM.Peel2.2", "TPM.Peel2.3")]))

write_xlsx(data.frame(Flesh1.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Flesh1-vs-FleshNC.xlsx")
write_xlsx(data.frame(Flesh2.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Flesh2-vs-FleshNC.xlsx")
write_xlsx(data.frame(Flesh5.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Flesh5-vs-FleshNC.xlsx")
write_xlsx(data.frame(Flesh2.vs.Flesh1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Flesh2-vs-Flesh1.xlsx")
write_xlsx(data.frame(Flesh5.vs.Flesh2.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Flesh5-vs-Flesh2.xlsx")
write_xlsx(data.frame(PeelNC.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/PeelNC-vs-FleshNC.xlsx")
write_xlsx(data.frame(Peel1.vs.PeelNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Peel1-vs-PeelNC.xlsx")
write_xlsx(data.frame(Peel2.vs.PeelNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Peel2-vs-PeelNC.xlsx")
write_xlsx(data.frame(Peel5.vs.PeelNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Peel5-vs-PeelNC.xlsx") 
write_xlsx(data.frame(Peel2.vs.Peel1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Peel2-vs-Peel1.xlsx")
write_xlsx(data.frame(Peel5.vs.Peel2.no.filter[,c(8,7,1:6,9:18)]), "No-filter/excel/Peel5-vs-Peel2.xlsx")

####################################################
write.table(data.frame(Flesh1.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Flesh1-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh2.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Flesh2-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh5.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Flesh5-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh2.vs.Flesh1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Flesh2-vs-Flesh1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Flesh5.vs.Flesh2.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Flesh5-vs-Flesh2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(PeelNC.vs.FleshNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/PeelNC-vs-FleshNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel1.vs.PeelNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Peel1-vs-PeelNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel2.vs.PeelNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Peel2-vs-PeelNC.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel5.vs.PeelNC.1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Peel5-vs-PeelNC.txt", quote = FALSE, row.names = FALSE, sep = "\t") 
write.table(data.frame(Peel2.vs.Peel1.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Peel2-vs-Peel1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(data.frame(Peel5.vs.Peel2.no.filter[,c(8,7,1:6,9:18)]), "No-filter/txt/Peel5-vs-Peel2.txt", quote = FALSE, row.names = FALSE, sep = "\t")

########################################################################################################
