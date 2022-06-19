REPOSITORY: Trabajo de Fin de Máster - Sara Braojos Juzgado

Here you'll find all scripts used to do the data analysis of this work. The code is an adaptation of the one you can find at:
https://github.com/pedrodelosreyes/Thesis_Code
https://github.com/pedrodelosreyes/starryRNA
https://github.com/fran-romero-campero/miscellanomics

Next, all scripts and their content are detailed:

sample_quality_alignment7d.sh -> Quality control, alignment and assembly for 7 DAG samples (FastQC, STAR,samtools, bamCoverage, stringtie).
sample_quality_alignment14d.sh -> Quality control, alignment and assembly for 14 DAG samples (FastQC, STAR,samtools, bamCoverage, stringtie).

quantification7d.sh -> Quantification for 7 DAG samples (featureCounts)
quantification14d.sh -> Quantification for 14 DAG samples (featureCounts)

limma_voom_7.R -> Differential expression (voom) for 7 DAG samples
limma_voom_14.R -> Differential expression (voom) for 14 DAG samples

go_terms_enrichment_7_activated_35sco_wt.R -> GO terms enrichment for activated genes in 35S:CO-wt set at 7 DAG
go_terms_enrichment_7_activated_co10_wt.R -> GO terms enrichment for activated genes in co-10-wt set at 7 DAG
go_terms_enrichment_7_repressed_35sco_wt.R -> GO terms enrichment for repressed genes in 35S:CO-wt set at 7 DAG
go_terms_enrichment_7_repressed_co10_wt.R -> GO terms enrichment for repressed genes in co-10-wt set at 7 DAG

go_terms_enrichment_14_activated_35sco_wt.R -> GO terms enrichment for activated genes in 35S:CO-wt set at 14 DAG
go_terms_enrichment_14_activated_co10_wt.R -> GO terms enrichment for activated genes in co-10-wt set at 14 DAG
go_terms_enrichment_14_repressed_35sco_wt.R -> GO terms enrichment for repressed genes in 35S:CO-wt set at 14 DAG
go_terms_enrichment_14_repressed_co10_wt.R -> GO terms enrichment for repressed genes in co-10-wt set at 14 DAG

comparaciones7-14dias.R -> To get common 7 ad 14 DAG genes, 7DAG-specific genes and 14DAG-specific genes for activated 35S:CO-wt and repressed co-10-wt gene sets

go_terms_enrichment_act_35sco_7_14.R -> GO terms enrichment for common 7 and 14 DAG activated genes in 35S:CO-wt sets
go_terms_enrichment_specific_act_35sco_7dias.R -> GO terms enrichment for 7DAG-specific activated genes in 35S:CO-wt sets
go_terms_enrichment_specific_act_35sco_14dias.R -> GO terms enrichment for 14DAG-specific activated genes in 35S:CO-wt sets

go_terms_enrichment_rep_co10_7_14.R -> GO terms enrichment for common 7 and 14 DAG repressed genes in co-10-wt sets
go_terms_enrichment_specific_rep_co10_7dias.R -> GO terms enrichment for 7DAG-specific repressed genes in co-10-wt sets
go_terms_enrichment_specific_rep_co10_14dias.R -> GO terms enrichment for 14DAG-specific repressed genes in co-10-wt sets

limma_voom_global7-14.R -> Differential expression (voom) for 7 DAG and 14 DAG samples together to get global heatmap and PCA

fc_barplots_7.R -> Make barplots with fold change values for specific 7DAG genes
fc_barplots_14.R -> Make barplots with fold change values for specific 14DAG genes