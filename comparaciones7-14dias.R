
# AUTHOR: Sara Braojos Juzgado

###########################################################################
################### GENES COMUNES A 7 y 14 días ###########################
###########################################################################
library(VennDiagram)

# Intersección activados 35SCO-wt a 7 y 14 días
act.35sco.wt.7 <- read.table(file="../7dias/tablas/activated_genes_35SCO_wt.txt")
act.35sco.wt.7 <- act.35sco.wt.7[,1]
act.35sco.wt.14 <- read.table(file="../14dias/tablas/activated_genes_35SCO_wt.txt")
act.35sco.wt.14 <- act.35sco.wt.14[,1]

act.35sco.7.14.dias <- intersect(act.35sco.wt.7, act.35sco.wt.14)
length(act.35sco.7.14.dias)

# write.table(act.35sco.7.14.dias, file="tablas/act_35sco_7_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

venn.diagram(
  x = list(act.35sco.wt.7, act.35sco.wt.14),
  category.names = c("Genes activados 35SCO 7 dias" , "Genes activados 35SCO 14 dias"),
  filename = 'images/venn_diagramm_act_35sco_7_14_dias.png',
  output=TRUE
)



# Intersección reprimidos 35SCO-wt a 7 y 14 días

rep.35sco.wt.7 <- read.table(file="../7dias/tablas/repressed_genes_35SCO_wt.txt")
rep.35sco.wt.7 <- rep.35sco.wt.7[,1]
rep.35sco.wt.14 <- read.table(file="../14dias/tablas/repressed_genes_35SCO_wt.txt")
rep.35sco.wt.14 <- rep.35sco.wt.14[,1]

rep.35sco.7.14.dias <- intersect(rep.35sco.wt.7, rep.35sco.wt.14)
length(rep.35sco.7.14.dias)

# write.table(rep.35sco.7.14.dias, file="tablas/rep_35sco_7_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

venn.diagram(
  x = list(rep.35sco.wt.7, rep.35sco.wt.14),
  category.names = c("Genes reprimidos 35SCO 7 dias" , "Genes reprimidos 35SCO 14 dias"),
  filename = 'images/venn_diagramm_rep_35sco_7_14_dias.png',
  output=TRUE
)


# Intersección activados co10-wt a 7 y 14 días

act.co10.wt.7 <- read.table(file="../7dias/tablas/activated_genes_co10_wt.txt")
act.co10.wt.7 <- act.co10.wt.7[,1]
act.co10.wt.14 <- read.table(file="../14dias/tablas/activated_genes_co10_wt.txt")
act.co10.wt.14 <- act.co10.wt.14[,1]

act.co10.7.14.dias <- intersect(act.co10.wt.7, act.co10.wt.14)
length(act.co10.7.14.dias)

# write.table(act.co10.7.14.dias, file="tablas/act_co10_7_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

venn.diagram(
  x = list(act.co10.wt.7, act.co10.wt.14),
  category.names = c("Genes activados co10 7 dias" , "Genes activados co10 14 dias"),
  filename = 'images/venn_diagramm_act_co10_7_14_dias.png',
  output=TRUE
)


# Intersección reprimidos co10-wt a 7 y 14 días

rep.co10.wt.7 <- read.table(file="../7dias/tablas/repressed_genes_co10_wt.txt")
rep.co10.wt.7 <- rep.co10.wt.7[,1]
rep.co10.wt.14 <- read.table(file="../14dias/tablas/repressed_genes_co10_wt.txt")
rep.co10.wt.14 <- rep.co10.wt.14[,1]

rep.co10.7.14.dias <- intersect(rep.co10.wt.7, rep.co10.wt.14)
length(rep.co10.7.14.dias)

# write.table(rep.co10.7.14.dias, file="tablas/rep_co10_7_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

venn.diagram(
  x = list(rep.co10.wt.7, rep.co10.wt.14),
  category.names = c("Genes reprimidos co10 7 dias" , "Genes reprimidos co10 14 dias"),
  filename = 'images/venn_diagramm_rep_co10_7_14_dias.png',
  output=TRUE
)

#######################################################################
############ GENES ESPECÍFICOS A 7 y 14 dias ##########################
#######################################################################

# función setdiff() para obtener genes específicos: The elements of setdiff(x,y) are those elements in x but not in y.

## Genes específicos activados 35sco-wt a 7 y 14 días

specific.act.35sco.wt.7 <- setdiff(act.35sco.wt.7, act.35sco.wt.14)
length(specific.act.35sco.wt.7)
# write.table(specific.act.35sco.wt.7, file="tablas/specific_act_35sco_7_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

specific.act.35sco.wt.14 <- setdiff(act.35sco.wt.14, act.35sco.wt.7)
length(specific.act.35sco.wt.14)
# write.table(specific.act.35sco.wt.14, file="tablas/specific_act_35sco_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)


## Genes específicos activados co10-wt a 7 y 14 días

specific.act.co10.wt.7 <- setdiff(act.co10.wt.7, act.co10.wt.14)
length(specific.act.co10.wt.7)
# write.table(specific.act.co10.wt.7, file="tablas/specific_act_co10_7_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

specific.act.co10.wt.14 <- setdiff(act.co10.wt.14, act.co10.wt.7)
length(specific.act.co10.wt.14)
# write.table(specific.act.co10.wt.14, file="tablas/specific_act_co10_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

## Genes específicos reprimidos 35sco-wt a 7 y 14 días

specific.rep.35sco.wt.7 <- setdiff(rep.35sco.wt.7, rep.35sco.wt.14)
length(specific.rep.35sco.wt.7)
# write.table(specific.rep.35sco.wt.7, file="tablas/specific_rep_35sco_7_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

specific.rep.35sco.wt.14 <- setdiff(rep.35sco.wt.14, rep.35sco.wt.7)
length(specific.rep.35sco.wt.14)
# write.table(specific.rep.35sco.wt.14, file="tablas/specific_rep_35sco_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)


## Genes específicos reprimidos co10-wt a 7 y 14 días

specific.rep.co10.wt.7 <- setdiff(rep.co10.wt.7, rep.co10.wt.14)
length(specific.rep.co10.wt.7)
# write.table(specific.rep.co10.wt.7, file="tablas/specific_rep_co10_7_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

specific.rep.co10.wt.14 <- setdiff(rep.co10.wt.14, rep.co10.wt.7)
length(specific.rep.co10.wt.14)
# write.table(specific.rep.co10.wt.14, file="tablas/specific_rep_co10_14_dias.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
