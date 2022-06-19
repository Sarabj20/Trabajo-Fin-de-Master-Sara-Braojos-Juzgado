## Go terms enrichment

##Authors:
# Fran Romero Campero
# Ana Belén Romero Losada
# Pedro de los Reyes Rodríguez
# Sara Braojos Juzgado

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.At.tair.db")
# BiocManager::install("rrvgo")
# install.packages("viridis")


library(clusterProfiler)
library(org.At.tair.db)
library(rrvgo)
library(ggplot2)
library(viridis)

## Auxiliary function to compute enrichments for GO table
compute.enrichments <- function(gene.ratios, bg.ratios)
{
  gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
  bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
  enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
  enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")
  
  return(enrichments.text)  
}

## Universo o background para el análisis. (Total de genes de Arabidopsis)
atha.universe <- unique(AnnotationDbi::select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR"))[["TAIR"]])
length(atha.universe)

## Leemos los genes que queremos analizar.
target.genes <- read.table(file="../../tablas/specific_act_35sco_7_dias.txt")
target.genes <- as.vector(target.genes[[1]])
length(target.genes)

## Calculamos el enriquecimiento
enrich.go <- enrichGO(gene          = target.genes ,
                      universe      = atha.universe,
                      OrgDb         = org.At.tair.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE,
                      keyType = "TAIR")

## Generate ouput table
enrich.go.result <- as.data.frame(enrich.go)

## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                           bg.ratios = enrich.go.result$BgRatio)

go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                              enrich.go.result$pvalue, enrich.go.result$qvalue,
                              go.term.enrichments, 
                              gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                              stringsAsFactors = FALSE)

colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                               "Enrichment (Target Ratio; BG Ration)","Genes")

#write.table(go.result.table, file = "go_terms_table_specific_act_35sco_7_dias.txt", sep = "\t", row.names = F, quote = F)

## Algunas representaciones gráficas
goplot(enrich.go,showCategory = 10)
barplot(enrich.go,drop=TRUE,showCategory = 40)  
emapplot(enrich.go)
cnetplot(enrich.go)
dotplot(enrich.go, showCategory =20)

png(filename = 'barplot_specific_act_35sco_7days.png', width =10, height = 15,units = 'in',res = 300)
barplot(enrich.go,drop=TRUE,showCategory = 40)

dev.off()

go.result.table[73,]
go.result.table[74,]

### Plotting specific GO terms #####
# GO:0009698 -> Phenylpropanoid metabolic process
# GO:0009699 -> phenylpropanoid biosynthetic process
# GO:0009808 -> lignin metabolic process
# GO:0009809 -> lignin biosynthetic process
# GO:0009813 -> flavonoid biosynthetic process
# GO:0009812 -> flavonoid metabolic process
# GO:0030026 -> cellular manganese ion homeostasis

key.terms <- c("GO:0009698", "GO:0009699", "GO:0009808", "GO:0009809",
               "GO:0009813", "GO:0009812", "GO:0030026")

keep <- go.result.table$`GO ID`[!go.result.table$`GO ID` %in% key.terms]
specific.go <- dropGO(enrich.go, term = keep )

png(filename = "barplot_activ_35sco_specific_7_selected.png", width = 10, height = 3, units = "in", res = 300)
barplot(specific.go,drop=TRUE,showCategory = 40, title = "Specific genes 7 DAG")
dev.off()

## Reducing the redundancy of GO using rrvgo package

# get the similarity matrix between terms
simMatrix <- calculateSimMatrix(go.result.table$`GO ID`,
                                orgdb="org.At.tair.db",
                                ont="BP",
                                method="Rel")

# group terms based on similarity
scores <- setNames(-log10(go.result.table$`q-value`), go.result.table$`GO ID`)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.At.tair.db")

# Plot similarity matrix as a heatmap
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

# Plot GO terms as scattered points
scatterPlot(simMatrix, reducedTerms)

# Treemaps are space-filling visualization of hierarchical structures
treemapPlot(reducedTerms)

png(filename = 'treemap_specific_act_35sco_7days.png', width =20, height = 10,units = 'in',res = 300)
treemapPlot(reducedTerms)

dev.off()

