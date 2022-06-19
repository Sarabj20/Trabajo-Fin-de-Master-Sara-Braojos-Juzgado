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
target.genes <- read.table(file="../tablas/repressed_genes_35SCO_wt_7d.txt")
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
# No hay enriquecimiento en nigún proceso biológico, no puedo continuar.

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

#write.table(go.result.table, file = "go_terms_table_repressed_35sco_wt.txt", sep = "\t", row.names = F, quote = F)

goplot(enrich.go,showCategory = 10)
barplot(enrich.go,drop=TRUE,showCategory = 20)  
emapplot(enrich.go)
cnetplot(enrich.go)
dotplot(enrich.go, showCategory =20)

