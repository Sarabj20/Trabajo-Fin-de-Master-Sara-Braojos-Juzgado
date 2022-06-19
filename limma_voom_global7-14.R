####################################################################
## Análisis de los resultados de STAR y featureCounts usando      ## 
## limma + voom                                                   ##
##                                                                ##
##               Aligner: STAR                                    ##
##               Quantification: featureCounts                    ##
##               Differential expression: limma-voom              ##
##               Some graphical representations: ggplot2          ##
##                                                                ##
##                                                                ##
## Autores: Pedro de los Reyes Rodriguez pedro.reyes@ibvf.csic.es ##
##          Francisco J. Romero Campero                           ##
####################################################################

## Required packages
library(edgeR) #also load limma as a dependency
library(ggplot2)
library(ggrepel)
library(sva) #batch correction
library(dplyr)
library(pheatmap)
library(readr)
library(viridis) #Default Color Maps
library(MetBrewer) #cool palette
library(DESeq2) #just to calculate fpkms
library(reshape2) #reorganize data frames
library(tidyr) #tidy messy data


# Leer las dos  matrices de samples.counts: 7 y 14 dias
sample.counts.7 <- read.table(file = '../7dias/sample_counts_7.txt')
head(sample.counts.7)
nrow(sample.counts.7)
sample.counts.14 <- read.table(file = '../14dias/sample_counts_14_corrected.txt')
head(sample.counts.14)
nrow(sample.counts.14)

#Unir matrices (cbind())
sample.counts.7.14 <- cbind(sample.counts.7,sample.counts.14)
colnames(sample.counts.7.14)
colnames(sample.counts.7.14) <- c("co10_1_7d", "co10_2_7d", "co10_3_7d",
                                  "Col0_1_7d", "Col0_2_7d", "Col0_3_7d",
                                  "x35SCO_1_7d", "x35SCO_2_7d", "x35SCO_3_7d",
                                  "co10_1_14d", "co10_2_14d", "co10_3_14d",
                                  "Col0_1_14d", "Col0_2_14d", "Col0_3_14d",
                                  "x35SCO_1_14d", "x35SCO_2_14d", "x35SCO_3_14d")

## Sample info. Same sample order in coldata than in count matrix
condition <- c(rep("co10_7d", 3),rep("Col0_7d",3), rep("x35SCO_7d",3), 
               rep("co10_14d", 3),rep("Col0_14d",3), rep("x35SCO_14d",3)) #Modificar
type <- c(rep("paired-end",9), rep("single-end", 9))
coldata <- data.frame(condition, type)
rownames(coldata) <- colnames(sample.counts.7.14)
coldata

##################################################################
########### Differential expression analysis ####################
################################################################

## Create DGEList object
d0 <- DGEList(sample.counts.7.14)
dim(d0)
## Calculate normalization factors. It doesn't normalize! 
##Just calculates normalization factors for use downstream
d0.norm <- calcNormFactors(d0) # tmm as default
# d0.norm <- calcNormFactors(d0, method = "TMM") 
# d0.norm <- calcNormFactors(d0, method = "upperquartile")

## Filter low expressed genes
cutoff <- 1
max.value <- apply(cpm(d0.norm), 1, max) #Maximo valor para cada gen = 1
drop <- which(max.value < cutoff) #Son los genes que no tienen al menos un valor de 1 CPM en alguna muestra
d <- d0.norm[-drop,] 
dim(d) #number of genes left

colnames(sample.counts.7.14)
genotype <- c(rep("co10_7d", 3),rep("Col0_7d",3), rep("x35SCO_7d",3), 
              rep("co10_14d", 3),rep("Col0_14d",3), rep("x35SCO_14d",3)) #MODIFICAR
group <- as.factor(genotype) #Create a new variable “group” as factor

plotMDS(d, col = as.numeric(group), labels = genotype) #Multidimensional scaling (MDS) plot


###### Voom transformation and calculation of variance weights
## Specify the model to be fitted. We do this before using voom since voom 
# uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + group)

# When operating on a DGEList-object, voom converts raw counts 
# to log-CPM values by automatically extracting library sizes 
# and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified 
# within voom using the normalize.method argument.

## ¿What is voom doing?
# Counts are transformed to log2 counts per million reads (CPM),
# where “per million reads” is defined based on the normalization 
# factors we calculated earlier
# A linear model is fitted to the log2 CPM for each gene, and the 
# residuals are calculated. A smoothed curve is fitted to the 
# sqrt(residual standard deviation) by average expression 
# (see red line in the plot). The smoothed curve is used to 
# obtain weights for each gene and sample that are passed 
# into limma along with the log2 CPMs.


# The read counts are processed by the voom function in limma to 
# convert them into log2 counts per million (logCPM) with associated 
# precision weights. If the data are very noisy, The logCPM values 
# can be normalized between samples by the voom function or can be 
# pre-normalized by adding normalization factors within edgeR (calcNormFactors).
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)


#########################################
#######        HEATMAP      #############
#########################################
# Draw heatmap for DEGS
voom.counts <- y$E
head(voom.counts) #counts transformed by voom 

#Leer todas las tablas de los DEG, almacenarlas en variables (total 8)
activated.genes.35SCO.7d <- read.table(file = "../7dias/tablas/activated_genes_35SCO_wt.txt")
activated.genes.35SCO.7d <- activated.genes.35SCO.7d$V1
repressed.genes.35SCO.7d <- read.table(file = "../7dias/tablas/repressed_genes_35SCO_wt.txt")
repressed.genes.35SCO.7d <- repressed.genes.35SCO.7d$V1
activated.genes.co10.7d <- read.table(file = "../7dias/tablas/activated_genes_co10_wt.txt")
activated.genes.co10.7d <- activated.genes.co10.7d$V1
repressed.genes.co10.7d <- read.table(file = "../7dias/tablas/repressed_genes_co10_wt.txt")
repressed.genes.co10.7d <- repressed.genes.co10.7d$V1

activated.genes.35SCO.14d <- read.table(file = "../14dias/tablas/activated_genes_35SCO_wt.txt")
activated.genes.35SCO.14d <- activated.genes.35SCO.14d$V1
repressed.genes.35SCO.14d <- read.table(file = "../14dias/tablas/repressed_genes_35SCO_wt.txt")
repressed.genes.35SCO.14d <- repressed.genes.35SCO.14d$V1
activated.genes.co10.14d <- read.table(file = "../14dias/tablas/activated_genes_co10_wt.txt")
activated.genes.co10.14d <- activated.genes.co10.14d$V1
repressed.genes.co10.14d <- read.table(file = "../14dias/tablas/repressed_genes_co10_wt.txt")
repressed.genes.co10.14d <- repressed.genes.co10.14d$V1

degs <- c(activated.genes.35SCO.7d, repressed.genes.35SCO.7d,
          activated.genes.co10.7d, repressed.genes.co10.7d,
          activated.genes.35SCO.14d, repressed.genes.35SCO.14d,
          activated.genes.co10.14d, repressed.genes.co10.14d)
length(degs)
degs <- unique(degs)
length(degs)

int.degs.tabla <- intersect(degs, rownames(voom.counts))
length(int.degs.tabla)

degs.table <- voom.counts[int.degs.tabla,]
#colnames(degs.table) <- colnames(voom.counts)
colnames(degs.table) <- c("co10_1_7d",  "co10_2_7d",  "co10_3_7d",
                          "Col-0_1_7d", "Col-0_2_7d", "Col-0_3_7d",
                          "35SCO_1_7d", "35SCO_2_7d", "35SCO_3_7d",
                          "co10_1_14d", "co10_2_14d", "co10_3_14d",
                          "Col-0_1_14d", "Col-0_2_14d", "Col-0_3_14d",
                          "35SCO_1_14d", "35SCO_2_14d", "35SCO_3_14d")
coldata

tiff("images/heatmap_degs.tiff", height = 6, width = 8, # Warning message: In grid.newpage() : unable to open TIFF file 'images/heatmap_degs.tiff'
     units = 'in', res=400, compression="lzw")
pheatmap(as.matrix(degs.table), cluster_rows = T, cluster_cols = T,
         scale = "row", clustering_distance_rows = "correlation", 
         clustering_method = "complete", #annotation_col = coldata, 
         main="DEGs",fontsize_col=14, fontsize_row = 4, color = magma(28),
         show_rownames = F)
dev.off()


## Heatmap for selected genes

# selected.genes <- read.table(file="selected_genes/list.txt", sep = "\t", header = T)
# head(selected.genes)
# 
# up.selected <- selected.genes$up
# length(up.selected)
# up.selected <- up.selected[1:22]
# up.table <- filtered.data[as.vector(up.selected),]
# tiff("selected_genes/genes_up_cluster.tiff", height = 4, width = 6, 
#      units = 'in', res=300, compression="lzw")
# pheatmap(as.matrix(up.table), cluster_rows = F, cluster_cols = T,
#          scale = "row", clustering_distance_rows = "correlation", 
#          clustering_method = "complete", #annotation_col = pheno.data, 
#          main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
# dev.off()
# 
# down.selected <- selected.genes$down
# length(down.selected)
# down.table <- filtered.data[as.vector(down.selected),]
# tiff("selected_genes/genes_down_cluster.tiff", height = 4, width = 6, 
#      units = 'in', res=300, compression="lzw")
# pheatmap(as.matrix(down.table), cluster_rows = F, cluster_cols = T,
#          scale = "row", clustering_distance_rows = "correlation", 
#          clustering_method = "complete", #annotation_col = pheno.data, 
#          main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
# dev.off()



##########################################
#######         PCA plot         #########
##########################################

## Draw PCA plot
# transpose the data and compute principal components
# pca.data <- prcomp(t(degs.table)) #Only for DEGs

# Log transform the data before performing PCA
# the log-transform is used to reduce the influence 
# of extreme values or outliers. 
pca.data <- prcomp(t(voom.counts))
nrow(pca.data$rotation)

# Calculate PCA component percentages
pca.data.perc <- round(100*pca.data$sdev^2/sum(pca.data$sdev^2),1)

# Extract 1 and 2 principle components and create a data frame with sample names, first and second principal components and group information
sample.names <- c("co10_1_7d",  "co10_2_7d",  "co10_3_7d",
                  "Col-0_1_7d", "Col-0_2_7d", "Col-0_3_7d",
                  "35SCO_1_7d", "35SCO_2_7d", "35SCO_3_7d",
                  "co10_1_14d", "co10_2_14d", "co10_3_14d",
                  "Col-0_1_14d", "Col-0_2_14d", "Col-0_3_14d",
                  "35SCO_1_14d", "35SCO_2_14d", "35SCO_3_14d")
df.pca.data <- data.frame(PC1 = pca.data$x[,1], PC2 = pca.data$x[,2], sample = sample.names, 
                          condition = c(rep("co10_7d", 3),rep("Col-0_7d",3), rep("35SCO_7d",3), 
                                        rep("co10_14d", 3),rep("Col-0_14d",3), rep("35SCO_14d",3)))
head(df.pca.data)
# tissue <- c("root", "air", "air")
# tissue <- rep(tissue, 3)
# df.pca.data$tissue <- tissue

# color by sample
ggplot(df.pca.data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+ 
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")")) 


# color by condition/group
tiff("images/PCAcondition.tiff", height = 6, width = 8,
     units = 'in', res=600, compression="lzw")
ggplot(df.pca.data, aes(PC1,PC2, color = condition))+
  geom_point(size=4)+
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")"))+
  geom_text_repel(aes(label=sample),point.padding = 0.5) 
# stat_ellipse(geom = "polygon", aes(fill=condition), alpha=0.2, show.legend = FALSE, level=0.95)
 dev.off()

