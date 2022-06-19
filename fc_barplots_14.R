##This script allows to make barplots with FOLD CHANGE values for specific genes
## Author: Pedro de los Reyes and Sara Braojos Juzgado
## Date: june 2022
## mail: pedro.reyes@ibvf.csic.es

library(ggplot2)

## Reading contrasts tables
x35SCO.wt <- read.table(file = "tablas/contraste_35SCO_wt_14.txt", sep = "\t", header = TRUE)
head(x35SCO.wt)
co10.wt <- read.table(file = "tablas/contraste_co10_wt_14.txt", sep = "\t", header = TRUE)
head(co10.wt)

##############################################################################
############# Flower development and meristem determinacy 7 DAG ##############
##############################################################################

## Getting the fc values for each gene for the two comparisons
co <- c(x35SCO.wt["AT5G15840",]$logFC, co10.wt["AT5G15840",]$logFC)
gi <- c(x35SCO.wt["AT1G22770",]$logFC, co10.wt["AT1G22770",]$logFC)
ap1 <- c(x35SCO.wt["AT1G69120",]$logFC, co10.wt["AT1G69120",]$logFC)
ap3 <- c(x35SCO.wt["AT3G54340",]$logFC, co10.wt["AT3G54340",]$logFC)
ft <- c(x35SCO.wt["AT1G65480",]$logFC, co10.wt["AT1G65480",]$logFC)
sep1 <- c(x35SCO.wt["AT5G15800",]$logFC, co10.wt["AT5G15800",]$logFC)
sep2 <- c(x35SCO.wt["AT3G02310",]$logFC, co10.wt["AT3G02310",]$logFC)
sep3 <- c(x35SCO.wt["AT1G24260",]$logFC, co10.wt["AT1G24260",]$logFC)
sep4 <- c(x35SCO.wt["AT2G03710",]$logFC, co10.wt["AT2G03710",]$logFC)
ful <- c(x35SCO.wt["AT5G60910",]$logFC, co10.wt["AT5G60910",]$logFC)
flc <- c(x35SCO.wt["AT5G10140",]$logFC, co10.wt["AT5G10140",]$logFC)
stm <- c(x35SCO.wt["AT1G62360",]$logFC, co10.wt["AT1G62360",]$logFC) # este es el único sólo de meristem determinacy

## merge fc values for the table
fc.values <- c(co, gi, ap1, ap3, ft, sep1, sep2, sep3, sep4, ful, flc, stm)

## set comparison column for the table
number.of.genes.to.plot <- 12
comparison <- rep(c("35SCO", "co10"), number.of.genes.to.plot)

## set gene column for the table
gene <- c(rep("CO",2), rep("GI", 2), rep("AP1", 2), rep("AP3", 2), rep("FT", 2),
          rep("SEP1", 2), rep("SEP2", 2), rep("SEP3", 2), rep("SEP4", 2), 
          rep("FUL", 2), rep("FLC", 2), rep("STM", 2))

# Build the table
fc.barplot.data <- data.frame(comparison, gene, fc.values )
head(fc.barplot.data)

## Plot the table
png(filename = "barplot_meristem_determinacy_and_flower_development_14d.png", width = 10, height = 4, units = "in", res = 300)
ggplot(fc.barplot.data, aes(x = gene, y=fc.values, fill=comparison)) +
  scale_x_discrete(limits=gene)+
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_manual(values=c("darkolivegreen3","darkorange2")) +
  labs(title = "Meristem determinacy and flower development 14 DAG") +
  xlab("") + 
  ylab("Log2(fold-change)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5)
dev.off()

###############################################################################
########################## Photoprotection  14 DAG ############################
###############################################################################

lhcb2.1 <- c(x35SCO.wt["AT2G05100",]$logFC, co10.wt["AT2G05100",]$logFC)
deg13 <- c(x35SCO.wt["AT3G27690",]$logFC, co10.wt["AT3G27690",]$logFC)
lhcb4.2 <- c(x35SCO.wt["AT3G08940",]$logFC, co10.wt["AT3G08940",]$logFC)
lhcb1.4 <- c(x35SCO.wt["AT2G34430",]$logFC, co10.wt["AT2G34430",]$logFC)
lhcb2.2 <- c(x35SCO.wt["AT2G05070",]$logFC, co10.wt["AT2G05070",]$logFC)
cab3 <- c(x35SCO.wt["AT1G29910",]$logFC, co10.wt["AT1G29910",]$logFC)
lhcb1.1 <- c(x35SCO.wt["AT1G29920",]$logFC, co10.wt["AT1G29920",]$logFC)
lhcb6 <- c(x35SCO.wt["AT1G15820",]$logFC, co10.wt["AT1G15820",]$logFC)
lhcb3 <- c(x35SCO.wt["AT5G54270",]$logFC, co10.wt["AT5G54270",]$logFC)
lhca4 <- c(x35SCO.wt["AT3G47470",]$logFC, co10.wt["AT3G47470",]$logFC)

## merge fc values for the table
fc.values <- c(lhcb2.1, deg13, lhcb4.2, lhcb1.4, lhcb2.2, cab3, lhcb1.1, lhcb6, lhcb3, lhca4)

## set comparison column for the table
number.of.genes.to.plot <- 10
comparison <- rep(c("35SCO", "co10"), number.of.genes.to.plot)

## set gene column for the table
gene <- c(rep("LHCB2.1",2), rep("DEG13", 2), rep("LHCB4.2", 2), rep("LHCB1.4", 2), rep("LHCB2.2", 2),
          rep("CAB3", 2), rep("LHCB1.1", 2), rep("LHCB6", 2), rep("LHCB3", 2), 
          rep("LHCA4", 2))

# Build the table
fc.barplot.data <- data.frame(comparison, gene, fc.values )
head(fc.barplot.data)

## Plot the table
png(filename = "barplot_photoprotection_14d.png", width = 10, height = 4, units = "in", res = 300)
ggplot(fc.barplot.data, aes(x = gene, y=fc.values, fill=comparison)) +
  scale_x_discrete(limits=gene)+
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_manual(values=c("darkolivegreen3","darkorange2")) +
  labs(title="Photoprotection 14 DAG") +
  xlab("") + 
  ylab("Log2(fold-change)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5)
dev.off()
###############################################################################
########################### Oxidative stress 14 DAG ###########################
###############################################################################

cat2 <- c(x35SCO.wt["AT4G35090",]$logFC, co10.wt["AT4G35090",]$logFC)
gi <- c(x35SCO.wt["AT1G22770",]$logFC, co10.wt["AT1G22770",]$logFC)
atpal1 <- c(x35SCO.wt["AT2G37040",]$logFC, co10.wt["AT2G37040",]$logFC)
atpal2 <- c(x35SCO.wt["AT3G53260",]$logFC, co10.wt["AT3G53260",]$logFC)
hsp70.13 <- c(x35SCO.wt["AT1G09080",]$logFC, co10.wt["AT1G09080",]$logFC)
atosa1 <- c(x35SCO.wt["AT5G64940",]$logFC, co10.wt["AT5G64940",]$logFC)
pdr11 <- c(x35SCO.wt["AT1G66950",]$logFC, co10.wt["AT1G66950",]$logFC)
prxcb <- c(x35SCO.wt["AT3G49120",]$logFC, co10.wt["AT3G49120",]$logFC)
prx53 <- c(x35SCO.wt["AT5G06720",]$logFC, co10.wt["AT5G06720",]$logFC)
prx37 <- c(x35SCO.wt["AT4G08770",]$logFC, co10.wt["AT4G08770",]$logFC)
prx52 <- c(x35SCO.wt["AT5G05340",]$logFC, co10.wt["AT5G05340",]$logFC)
lhcb2.2 <- c(x35SCO.wt["AT2G05070",]$logFC, co10.wt["AT2G05070",]$logFC)


## merge fc values for the table
fc.values <- c(cat2, gi, atpal1, atpal2, hsp70.13, atosa1, 
               pdr11, prxcb, prx53, prx37, prx52, lhcb2.2)

## set comparison column for the table
number.of.genes.to.plot <- 12
comparison <- rep(c("35SCO", "co10"), number.of.genes.to.plot)

## set gene column for the table
gene <- c(rep("CAT2",2), rep("GI", 2), rep("ATPAL1", 2), rep("ATPAL2", 2), rep("HSP70-13", 2),
          rep("ATOSA1", 2), rep("PDR11", 2), rep("PRXCB", 2), 
          rep("PRX53", 2), rep("PRX37", 2), rep("PRX52", 2), rep("LHCB2.2", 2))

# Build the table
fc.barplot.data <- data.frame(comparison, gene, fc.values )
head(fc.barplot.data)

## Plot the table
png(filename = "barplot_oxidative_stress_14d.png", width = 10, height = 4, units = "in", res = 300)
ggplot(fc.barplot.data, aes(x = gene, y=fc.values, fill=comparison)) +
  scale_x_discrete(limits=gene)+
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_manual(values=c("darkolivegreen3","darkorange2")) +
  labs(title="Oxidative stress 14 DAG") +
  xlab("") + 
  ylab("Log2(fold-change)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5)
dev.off()
###############################################################################
####################### PIGMENT BIOSYNTHESIS 14 DAG ###########################
###############################################################################


mybd <- c(x35SCO.wt["AT1G70000",]$logFC, co10.wt["AT1G70000",]$logFC)
zep <- c(x35SCO.wt["AT5G67030",]$logFC, co10.wt["AT5G67030",]$logFC)
bch2 <- c(x35SCO.wt["AT5G52570",]$logFC, co10.wt["AT5G52570",]$logFC)
cao <- c(x35SCO.wt["AT1G44446",]$logFC, co10.wt["AT1G44446",]$logFC)
cdr1 <- c(x35SCO.wt["AT3G56940",]$logFC, co10.wt["AT3G56940",]$logFC)
elip2 <- c(x35SCO.wt["AT4G14690",]$logFC, co10.wt["AT4G14690",]$logFC)
atom1 <- c(x35SCO.wt["AT5G54160",]$logFC, co10.wt["AT5G54160",]$logFC)
pr5 <- c(x35SCO.wt["AT1G75040",]$logFC, co10.wt["AT1G75040",]$logFC)


## merge fc values for the table
fc.values <- c(mybd, zep, bch2, cao, cdr1, elip2, atom1, pr5)

## set comparison column for the table
number.of.genes.to.plot <- 8
comparison <- rep(c("35SCO", "co10"), number.of.genes.to.plot)

## set gene column for the table
gene <- c(rep("MYBD",2), rep("ZEP", 2), rep("BCH2", 2), rep("CAO", 2), rep("CDR1", 2),
          rep("ELIP2", 2), rep("ATOM1", 2), rep("PR5", 2))

# Build the table
fc.barplot.data <- data.frame(comparison, gene, fc.values )
head(fc.barplot.data)

## Plot the table
png(filename = "barplot_pigment_biosynthesis_14d.png", width = 10, height = 4, units = "in", res = 300)
ggplot(fc.barplot.data, aes(x = gene, y=fc.values, fill=comparison)) +
  scale_x_discrete(limits=gene)+
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_manual(values=c("darkolivegreen3","darkorange2")) +
  labs(title="Pigment biosynthesis 14 DAG") +
  xlab("") + 
  ylab("Log2(fold-change)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5)
dev.off()

