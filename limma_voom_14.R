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
##          Sara Braojos Juzgado                                  ##
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
library(reshape2) #reorganize data frames
library(tidyr) #tidy messy data

## Import the read count matrix data into R.
counts <- read.table(file="counts14.txt", sep = "\t", header = T)
rownames(counts) <- counts$Geneid
head(counts)
counts["AT5G15840",]#CO
counts["AT1G65480",]#FT


## Get only the counts
samples <- c(paste("co10", 1:3, sep = "_"), paste("Col-0", 1:3, sep = "_"), 
             paste("35SCO", 1:3, sep = "_"))
sample.counts <- counts[,7:15]
rownames(sample.counts) <- counts$Geneid
head(sample.counts,2)
colnames(sample.counts) <- samples
#write.table(sample.counts, file = "sample_counts_14.txt", sep = "\t", quote = FALSE)


## Library sizes
barplot(colSums(sample.counts)*1e-6, names = colnames(sample.counts), 
        ylab="Library size (millions)",las = 2, cex.names = 0.7)


## Sample info. Same sample order in coldata than in count matrix
condition <- c(rep("co10", 3),rep("col0",3), rep("x35SCO",3))
type <- rep("single-end",9)
coldata <- data.frame(condition, type)
rownames(coldata) <- colnames(sample.counts)
coldata


## Visualize similarity between replicates 
plot(log2(sample.counts[,1]+1),log2(sample.counts[,2]+1),pch=19,cex=0.7,xlab="co10_1",ylab="co10_2",cex.lab=1.25)
plot(log2(sample.counts[,2]+1),log2(sample.counts[,3]+1),pch=19,cex=0.7,xlab="co10_2",ylab="co10_3",cex.lab=1.25)
plot(log2(sample.counts[,1]+1),log2(sample.counts[,3]+1),pch=19,cex=0.7,xlab="co10_1",ylab="co10_3",cex.lab=1.25)

plot(log2(sample.counts[,4]+1),log2(sample.counts[,5]+1),pch=19,cex=0.7,xlab="wt_1",ylab="wt_2",cex.lab=1.25)
plot(log2(sample.counts[,5]+1),log2(sample.counts[,6]+1),pch=19,cex=0.7,xlab="wt_2",ylab="wt_3",cex.lab=1.25)
plot(log2(sample.counts[,4]+1),log2(sample.counts[,6]+1),pch=19,cex=0.7,xlab="wt_1",ylab="wt_3",cex.lab=1.25)

plot(log2(sample.counts[,7]+1),log2(sample.counts[,8]+1),pch=19,cex=0.7,xlab="35SCO_1",ylab="35SCO_2",cex.lab=1.25)
plot(log2(sample.counts[,8]+1),log2(sample.counts[,9]+1),pch=19,cex=0.7,xlab="35SCO_2",ylab="35SCO_3",cex.lab=1.25)
plot(log2(sample.counts[,7]+1),log2(sample.counts[,9]+1),pch=19,cex=0.7,xlab="35SCO_1",ylab="35SCO_3",cex.lab=1.25)

#Vemos que las counts de la réplica 1 en las tres condiciones se desvían de la 
#línea de tendencia, estamos ante un efecto lote debido a que en el caso de la 
#réplica 1 se obtuvo RNA de la planta completa (también raíz), mientras que de
#las otras dos réplicas se obtivo sólo de la parte aérea (hojas).


##########################################################
######### Batch correction using Combat-seq #############
########################################################

# devtools::install_github("zhangyuqing/sva-devel")
# https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/library(sva)
tissue <- rep(c("root", "noroot", "noroot"),3)
combat1 <- ComBat_seq(counts = as.matrix(sample.counts), batch = tissue, group = NULL) 
# write.table(combat1, file = "combat1.txt", sep = "\t", quote = FALSE)


## Use batch corrected counts from now on
sample.counts <- as.data.frame(combat1) 

## Important: check again the scatterplots between replicates

## Boxplot to check global sample distributions are similar and comparable.
boxplot(log2(sample.counts+1),col=met.brewer(n=9, name="Cassatt2"),ylab="log2(counts + 1)",cex.lab=1.5)

##################################################################
########### Differential expression analysis ####################
################################################################

## Create DGEList object
d0 <- DGEList(sample.counts)
dim(d0)
## Calculate normalization factors (it doesn't normalize) 
d0.norm <- calcNormFactors(d0) # TMM as default


## Filter low expressed genes
cutoff <- 1
max.value <- apply(cpm(d0.norm), 1, max) #Maximo valor para cada gen = 1
drop <- which(max.value < cutoff) #Son los genes que no tienen al menos un valor de 1 CPM en alguna muestra
d <- d0.norm[-drop,] 
dim(d) #number of genes left

colnames(sample.counts)
genotype <- c(rep("co10", 3),rep("col0",3), rep("x35SCO",3))
group <- as.factor(genotype) #Create a new variable “group” as factor

plotMDS(d, col = as.numeric(group), labels = genotype) #Multidimensional scaling (MDS) plot


########## Voom transformation and calculation of variance weights ##########

## Specify the model to be fitted. 
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
nofilter.voom <- voom(d0, mm, plot = T) #Sin filtrar

#png(filename = 'voom_mean_variance_trend_14d.png', width =10, height = 7,units = 'in',res = 300)
#y <- voom(d, mm, plot = T)
#dev.off()


## Boxplot after voom WITHOUT normalization
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
y.no.norm <- voom(d, mm, plot = T)

#png(filename = 'boxplot_no_normalized_data_14days.png', width =10, height = 7,units = 'in',res = 300)
boxplot(y.no.norm$E, col=met.brewer(n=9, name="Hiroshige"), ylab="log2CPM voom",cex.lab=1.2, main= "No normalized data")
#dev.off()


## Boxplot after voom WITH normalization.

#png(filename = 'boxplot_normalized_data_14days.png', width =10, height = 7,units = 'in',res = 300)
boxplot(y$E, col=met.brewer(n=9, name="Hiroshige"),ylab="log2CPM voom",cex.lab=1.2, main= "Normalized data")
#dev.off()


## Plot density
reshaped.y <- melt(y$E) #just reshaping data frame to input to ggplot
p <- ggplot(aes(x=value, colour=Var2), data=reshaped.y)
p + geom_density() + xlab("Log2-cpm") + ylab("Density")


## Ajuste lineal
fit <- lmFit(y, mm)
head(coef(fit))

contrast.matrix <- makeContrasts(groupx35SCO-groupcol0, groupco10-groupcol0,
                                 levels=colnames(coef(fit))) #Hace contrastes
head(contrast.matrix)

contrast.linear.fit <- contrasts.fit(fit, contrast.matrix) #Ajusta contrastes
contrast.results <- eBayes(contrast.linear.fit) #Calcula la bondad del ajuste

###########################################
#######Contraste 35SCO vs wt (Col0)########
###########################################
x35SCO.wt <- topTable(contrast.results, number=nrow(contrast.results),coef=1,sort.by="logFC")
head(x35SCO.wt)
#write.table(x35SCO.wt, file = "tablas/contraste_35SCO_wt_14.txt", sep = "\t", quote = FALSE) #Para gráfico expresión FC de los genes

x35SCO.wt["AT1G65480",] #FT
x35SCO.wt["AT5G15840",] #CO
nrow(x35SCO.wt)

fold.change.35SCO.wt <- x35SCO.wt$logFC
genes.ids.35SCO.wt <- rownames(x35SCO.wt)
p.value <- x35SCO.wt$adj.P.Val

activated.genes.35SCO <- genes.ids.35SCO.wt[fold.change.35SCO.wt > log2(2) & p.value < 0.05] 
repressed.genes.35SCO <- genes.ids.35SCO.wt[fold.change.35SCO.wt < -log2(2) & p.value < 0.05]
length(activated.genes.35SCO)
length(repressed.genes.35SCO)
# write.table(activated.genes.35SCO, file="tablas/activated_genes_35SCO_wt_14d.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.35SCO, file="tablas/repressed_genes_35SCO_wt_14d.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

x35SCO.wt.activated.info <- x35SCO.wt[activated.genes.35SCO,]
x35SCO.wt.repressed.info <- x35SCO.wt[repressed.genes.35SCO,]

# write.table(x35SCO.wt.activated.info, file="tablas/x35SCO_wt_activated_info_14d.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# write.table(x35SCO.wt.repressed.info, file="tablas/x35SCO_wt_repressed_info_14d.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


############################ Scatterplot 35SCO vs wt ############################
voom.counts <- y$E
co10.counts<- (voom.counts[,1] + voom.counts[,2] + voom.counts[,3])/3
col0.counts<- (voom.counts[,4] + voom.counts[,5] + voom.counts[,6])/3
x35SCO.counts<- (voom.counts[,7] + voom.counts[,8] + voom.counts[,9])/3

names(x35SCO.counts) <- rownames(voom.counts)
names(col0.counts) <- rownames(voom.counts)
names(co10.counts) <- rownames(voom.counts)
mean.counts <- matrix(c(co10.counts, col0.counts, x35SCO.counts), ncol=3)
colnames(mean.counts) <- c("co10", "Col0", "x35SCO")
rownames(mean.counts) <- rownames(voom.counts)
head(mean.counts)
mean.counts["AT5G15840",]

mean.counts.1.2 <- as.data.frame(mean.counts[,1:2]) # co10 vs. col0
mean.counts.2.3 <- as.data.frame(mean.counts[,2:3]) # 35SCO vs. col0


ggplot(data = as.data.frame(mean.counts.2.3), aes(x=Col0, y=x35SCO)) + geom_point() + theme_minimal()

# add a column of NAs
mean.counts.2.3$diffexpressed <- "NO"
x35SCO.wt$diffexpressed <- "NO"
# if log2Foldchange > 1  set as "UP" 
for (i in 1:nrow(x35SCO.wt))
{
  if (x35SCO.wt[i,"logFC"]>1)
  {
    up.gene <- rownames(x35SCO.wt[i,])
    mean.counts.2.3[up.gene,"diffexpressed"] <- "UP"
    x35SCO.wt[up.gene,"diffexpressed"] <- "UP"
  }
}
# if log2Foldchange < -1  set as "DOWN"
for(i in 1:nrow(x35SCO.wt))
{
  if (x35SCO.wt[i,"logFC"] < -1)
  {
    down.gene <- rownames(x35SCO.wt[i,])
    mean.counts.2.3[down.gene,"diffexpressed"] <- "DOWN"
    x35SCO.wt[down.gene,"diffexpressed"] <- "DOWN"
  }
}
# Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

p <- ggplot(data = as.data.frame(mean.counts.2.3), aes(x=Col0, y=x35SCO, col=diffexpressed)) + geom_point() + theme_minimal() +scale_color_manual(values=mycolors)
p + ggtitle("Scatterplot") + xlab("log2(CPM) Col0") + ylab("log2(CPM) 35SCO")


## Quick Volcano plot
gene.names <- rownames(contrast.results$coefficients)
volcanoplot(contrast.results, coef=1, names = gene.names, highlight = 5)

## MA plot: To check whether data are comparable among groups
limma::plotMA(contrast.results, coef=1, main="35SCO vs WT")
plotWithHighlights(x35SCO.wt$AveExpr, x35SCO.wt$logFC,
                   status=x35SCO.wt$diffexpressed,
                   hl.col = c("red", "blue"), hl.cex = 0.5,
                   bg.col = "grey", xlab="Average log2-expression",
                   ylab= "log2 FC", main = "35SCO vs WT")
abline(h=0, col="darkred", lwd=2)



#################################################
###############Contraste co10 vs wt#############
################################################
co10.wt <- topTable(contrast.results, number=nrow(contrast.results),coef=2,sort.by="logFC")
#write.table(co10.wt, file = "tablas/contraste_co10_wt_14.txt", sep = "\t", quote = FALSE) #Para gráfico expresión FC de los genes
# co10.wt <- topTable(contrast.results, number=33602,coef=2,sort.by="logFC", p.value = 0.05)
head(co10.wt)

co10.wt["AT5G15840",]#CO
co10.wt["AT1G65480",]#FT


fold.change.co10.wt <- co10.wt$logFC
p.value <- co10.wt$adj.P.Val
genes.ids.co10.wt <- rownames(co10.wt)

activated.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt > log2(1.5)]
repressed.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt < -log2(1.5)]
# activated.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt > log2(1.5) & p.value < 0.05 ]
# repressed.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt < -log2(1.5) & p.value < 0.05]

#Como las diferencias entre co10 y col0 son mucho más sutiles que entre 35sco y col0, 
#en este caso sólo se establece el umbral con el foldchange.
length(activated.genes.co10)
length(repressed.genes.co10)

# write.table(activated.genes.co10, file="tablas/activated_genes_co10_wt_14d.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.co10, file="tablas/repressed_genes_co10_wt_14d.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

co10.wt.activated.info <- co10.wt[activated.genes.co10,]
co10.wt.repressed.info <- co10.wt[repressed.genes.co10,]

# write.table(co10.wt.activated.info, file="tablas/co10_wt_activated_info_14d.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# write.table(co10.wt.repressed.info, file="tablas/co10_wt_repressed_info_14d.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


###################### Scatterplot co10 vs WT #################################
ggplot(data = as.data.frame(mean.counts.1.2), aes(x=Col0, y=co10)) + geom_point() + theme_minimal()

# add a column of NAs
mean.counts.1.2$diffexpressed <- "NO"
co10.wt$diffexpressed <- "NO"
# if log2Foldchange > 1  set as "UP" 
for (i in 1:nrow(co10.wt))
{
  if (co10.wt[i,"logFC"]>1)
  {
    up.gene <- rownames(co10.wt[i,])
    mean.counts.1.2[up.gene,"diffexpressed"] <- "UP"
    co10.wt[up.gene,"diffexpressed"] <- "UP"
  }
}
# if log2Foldchange < -1  set as "DOWN"
for(i in 1:nrow(co10.wt))
{
  if (co10.wt[i,"logFC"] < -1)
  {
    down.gene <- rownames(co10.wt[i,])
    mean.counts.1.2[down.gene,"diffexpressed"] <- "DOWN"
    co10.wt[down.gene,"diffexpressed"] <- "DOWN"
  }
}
# Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

p <- ggplot(data = as.data.frame(mean.counts.1.2), aes(x=Col0, y=co10, col=diffexpressed)) + geom_point() + theme_minimal() +scale_color_manual(values=mycolors)
p + ggtitle("Scatterplot") + xlab("log2(CPM) Col-0") + ylab("log2(CPM) co-10")


## Quick Volcano plot 
gene.names <- rownames(contrast.results$coefficients)
volcanoplot(contrast.results, coef=2, names = gene.names, highlight = 5)

## MA plot : To check whether data are comparable among groups
limma::plotMA(contrast.results, coef=2, main="co10 vs WT")
plotWithHighlights(co10.wt$AveExpr, co10.wt$logFC,
                   status=co10.wt$diffexpressed,
                   hl.col = c("red", "blue"), hl.cex = 0.5,
                   bg.col = "grey", xlab="Average log2-expression",
                   ylab= "log2 FC", main = "co-10 vs WT")
abline(h=0, col="darkred", lwd=2)

##############################################
##########      Volcano plots       ###########
##############################################

#### 35SCO vs WT
volcano <- ggplot(data=x35SCO.wt, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
x35SCO.wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
x35SCO.wt$diffexpressed[x35SCO.wt$logFC > log2(2) & x35SCO.wt$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
x35SCO.wt$diffexpressed[x35SCO.wt$logFC < -log2(2) & x35SCO.wt$P.Value< 0.05] <- "DOWN"

#Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add a color in aes section
volcano <- ggplot(data=x35SCO.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  scale_color_manual(values=mycolors) 

x35SCO.wt$agi <- NA
x35SCO.wt$agi[x35SCO.wt$diffexpressed != "NO"] <- rownames(x35SCO.wt)[x35SCO.wt$diffexpressed != "NO"]

png(filename = 'volcanoplot_35sco_wt_14d.png', width =7, height = 5,units = 'in',res = 300)
ggplot(data=x35SCO.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=agi)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 0) +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()

#### co10 vs WT
volcano <- ggplot(data=co10.wt, aes(x=logFC, y=-log10(P.Value))) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
co10.wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
co10.wt$diffexpressed[co10.wt$logFC > log2(2) & co10.wt$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
co10.wt$diffexpressed[co10.wt$logFC < -log2(2) & co10.wt$P.Value < 0.05] <- "DOWN"

#Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add a color in aes section
volcano <- ggplot(data=co10.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  scale_color_manual(values=mycolors)

co10.wt$agi <- NA
co10.wt$agi[co10.wt$diffexpressed != "NO"] <- rownames(co10.wt)[co10.wt$diffexpressed != "NO"]

png(filename = 'volcanoplot_co10_wt_14d.png', width =7, height = 5,units = 'in',res = 300)
ggplot(data=co10.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=agi)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 0) +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-log2(1.2), log2(1.2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()


#########################################
#######        HEATMAP      #############
#########################################

# Draw heatmap for DEGS
head(voom.counts) #counts transformed by voom 

degs <- c(activated.genes.35SCO, repressed.genes.35SCO,
          activated.genes.co10, repressed.genes.co10)
length(degs)
degs <- unique(degs)
length(degs)

degs.table <- voom.counts[degs,]
colnames(degs.table) <- colnames(voom.counts)
coldata

#tiff("images/heatmap_degs.tiff", height = 4, width = 8, units = 'in', res=400, compression="lzw")
pheatmap(as.matrix(degs.table), cluster_rows = T, cluster_cols = T,
         scale = "row", clustering_distance_rows = "correlation", 
         clustering_method = "complete",  
         main="DEGs",fontsize_col=14, fontsize_row = 4, color = magma(28),
         show_rownames = F)
#dev.off()


##########################################
#######         PCA plot         #########
##########################################

# Draw PCA plot
pca.data <- prcomp(t(voom.counts))
nrow(pca.data$rotation)

# Calculate PCA component percentages
pca.data.perc <- round(100*pca.data$sdev^2/sum(pca.data$sdev^2),1)

# Extract 1 and 2 principle components and create a data frame with sample names, first and second principal components and group information
df.pca.data <- data.frame(PC1 = pca.data$x[,1], PC2 = pca.data$x[,2], sample = colnames(sample.counts), 
                          condition = c(rep("35SCO",3), rep("Col0",3), rep("co10",3)))
head(df.pca.data)
# tissue <- c("root", "air", "air")
# tissue <- rep(tissue, 3)
# df.pca.data$tissue <- tissue

# color by sample
ggplot(df.pca.data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+ 
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")")) 


# color by condition/group
# tiff("images/PCAcondition.tiff", height = 4, width = 6, units = 'in', res=600, compression="lzw")
ggplot(df.pca.data, aes(PC1,PC2, color = condition))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")"))+
  geom_text_repel(aes(label=sample),point.padding = 0.75) 
# dev.off()
