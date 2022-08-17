library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(tidyverse)
library(RColorBrewer)

obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA_no_cell_cycle/snRNA_combo_Merged_is_Myeloid_lin_44_samples.rds")
dams_cell_state <- fread('dams_cell_state.tsv')
genes.toplot <- dams_cell_state %>% filter(avg_log2FC>2 & p_val_adj < 0.01 & cluster =='Mast') %>% pull(TF)

dir.create('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/DAM/v3.0/combo_only/fromRNA')
setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/DAM/v3.0/combo_only/fromRNA')


walk(genes.toplot, function(gene) {
  FeaturePlot(obj, features = gene, order = T)  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0('Mast','_FeaturePlot_',gene,'_dcis_idc.pdf'), width = 5, height = 4)
})


#############
genes.toplot <- dams_cell_state %>% filter(avg_log2FC>2 & p_val_adj < 0.01 & cluster =='DC') %>% pull(TF) %>% unique
if (sum(grepl('::', genes.toplot))>0) {
  genes.toplot <- c(str_split_fixed(genes.toplot, '::', 2)[,1],str_split_fixed(genes.toplot, '::', 2)[,2]) %>% unique()
  genes.toplot <- genes.toplot[genes.toplot!=""]
  
}
genes.toplot <- genes.toplot[!grepl('var|-', genes.toplot)]
walk(genes.toplot, function(gene) {
  FeaturePlot(obj, features = gene, order = T)  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0('cDC1','_FeaturePlot_',gene,'_TF_expresion.pdf'), width = 5, height = 4)
})

DefaultAssay(panc) <- 'ATAC_immune'
MotifPlot(
  object = panc,
  motifs = (motif.tf %>% filter(grepl('TBX', TF)) %>% pull(motif) %>% sort )
)

########
obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA_no_cell_cycle/snRNA_combo_Merged_is_Bcell_lin_44_samples.rds")

genes.toplot <- dams_cell_state %>% filter(avg_log2FC>2 & p_val_adj < 0.01 & cluster =='Plasma cells') %>% pull(TF) %>% unique
if (sum(grepl('::', genes.toplot))>0) {
  genes.toplot <- c(str_split_fixed(genes.toplot, '::', 2)[,1],str_split_fixed(genes.toplot, '::', 2)[,2]) %>% unique()
  genes.toplot <- genes.toplot[genes.toplot!=""]
  
}
genes.toplot <- genes.toplot[!grepl('var|-', genes.toplot)]
walk(genes.toplot, function(gene) {
  FeaturePlot(obj, features = gene, order = T)  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0('Plasma','_FeaturePlot_',gene,'_TF_expresion.pdf'), width = 5, height = 4)
})



###########
obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA_no_cell_cycle/snRNA_combo_Merged_is_Tcell_lin_44_samples.rds")
dams_cell_general$cluster %>% unique
genes.toplot <- dams_cell_general %>% filter(avg_log2FC>2 & p_val_adj < 0.01 & cluster =="NK cells") %>% pull(TF) %>% unique
if (sum(grepl('::', genes.toplot))>0) {
  genes.toplot <- c(str_split_fixed(genes.toplot, '::', 2)[,1],str_split_fixed(genes.toplot, '::', 2)[,2]) %>% unique()
  genes.toplot <- genes.toplot[genes.toplot!=""]
  
}
genes.toplot <- genes.toplot[!grepl('var|-', genes.toplot)]
walk(genes.toplot, function(gene) {
  FeaturePlot(obj, features = gene, order = T)  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("NK cells",'_FeaturePlot_',gene,'_TF_expresion.pdf'), width = 5, height = 4)
})

DimPlot(obj, group.by = 'cell_type_state', cols = 'Spectral')
