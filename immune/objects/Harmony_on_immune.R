library(Signac)
library(Seurat)
library(data.table)
library(tidyverse)
library(harmony)

setwd('/diskmnt/Projects/snATAC_analysis/immune/obj/v3.0')
panc <- readRDS('Reclustered_immune_snATAC_Merged_137_samples_v3.0_cluster_peaks.rds')
panc$tumor_type <- case_when(panc$Cancer=='MM' ~ 'Liquid',
                             TRUE ~ 'Solid')
combo.meta <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_is_Tcell_lin_44_samples.res2.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)

add_filename <- '137_samples_v3.0_cluster_peaks_harmony'
panc2 <- panc %>%
  RunHarmony( "tumor_type", reduction = 'lsi', assay.use = 'ATAC_immune',  project.dim = F) %>%
  FindNeighbors(reduction = "harmony", dims = 2:30) %>%
  FindClusters(verbose = TRUE, resolution = 2, algorithm = 3) %>%
  RunUMAP(reduction = "harmony", dims = 2:30)

p <- DimPlot(panc2, group.by = 'Cancer', pt.size = 0.1)
pdf(paste0("Dimplot_cancer_reclustered_immune_snATAC_Merged_", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p <- DimPlot(panc2, label = T, pt.size = 0.1)
pdf(paste0("Dimplot_clusters_reclustered_immune_snATAC_Merged_", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p3 <- DimPlot(panc2,group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_cell_type_reclustered_immune_snATAC_Merged", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

p1 <- DimPlot(panc,reduction = 'lsi',dims = 2:3, group.by = 'Cancer', pt.size = 0.1)
pdf(paste0("PCplot_cancer_reclustered_immune_snATAC_Merged_", '137_samples_v3.0_cluster_peaks', ".pdf"),height=10,width=12, useDingbats = F)
print(p1)
dev.off()

p2 <- DimPlot(panc2,reduction = 'harmony',dims = 2:3, group.by = 'Cancer', pt.size = 0.1)
pdf(paste0("PCplot_cancer_reclustered_immune_snATAC_Merged", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p2)
dev.off()

p2 <- DimPlot(panc2,reduction = 'harmony',dims = 2:3, group.by = 'data.type', pt.size = 0.1)
pdf(paste0("PCplot_data.type_reclustered_immune_snATAC_Merged", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p2)
dev.off()

combo.meta <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_is_Tcell_lin_44_samples.res2.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)

panc2 <- AddMetaData(panc2, combo.meta[,c('Cell_type_state', 'Cell_type_markers')])
p3 <- DimPlot(panc2,group.by = 'Cell_type_state', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_cell_type_reclustered_immune_snATAC_Merged", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

p3 <- DimPlot(panc2,group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_cell_type.harmonized.cancer_reclustered_immune_snATAC_Merged", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

p3 <- DimPlot(panc2,group.by = 'data.type', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_data_type_reclustered_immune_snATAC_Merged",add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()


##### use harmony to remove data type batch effect
panc3 <- panc %>%
  RunHarmony( "data.type", reduction = 'lsi', assay.use = 'ATAC_immune',  project.dim = F) %>%
  FindNeighbors(reduction = "harmony", dims = 2:30) %>%
  FindClusters(verbose = TRUE, resolution = 2) %>%
  RunUMAP(reduction = "harmony", dims = 2:30)
add_filename3 <- '_137_samples_v3.0_cluster_peaks_harmony_data.type'
p <- DimPlot(panc3, group.by = 'Cancer', pt.size = 0.1)
pdf(paste0("Dimplot_cancer_reclustered_immune_snATAC_Merged_", add_filename3, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p <- DimPlot(panc3, label = T, pt.size = 0.1)
pdf(paste0("Dimplot_clusters_reclustered_immune_snATAC_Merged_", add_filename3, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p3 <- DimPlot(panc3,group.by = 'data.type', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_data_type_reclustered_immune_snATAC_Merged",add_filename3, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

panc3 <- AddMetaData(panc3, combo.meta[,c('Cell_type_state', 'Cell_type_markers')])
p3 <- DimPlot(panc3,group.by = 'Cell_type_state', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_Cell_type_state_reclustered_immune_snATAC_Merged", add_filename3, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

p3 <- DimPlot(panc3,group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_cell_type.harmonized.cancer_reclustered_immune_snATAC_Merged", add_filename3, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

DepthCor(panc3, reduction = 'harmony')
ggsave(paste0("Corplot_", add_filename3, ".pdf"),height=5,width=6, useDingbats = F)

DepthCor(panc3, reduction = 'lsi')
ggsave(paste0("Corplot_", add_filename3, '_lsi', ".pdf"),height=5,width=6, useDingbats = F)

##################### now remove the effect of both MM and data.type
panc4 <- panc %>%
  RunHarmony( c("data.type", "tumor_type"), reduction = 'lsi', assay.use = 'ATAC_immune',  project.dim = F, max.iter.harmony = 20) %>%
  FindNeighbors(reduction = "harmony", dims = 2:30) %>%
  FindClusters(verbose = TRUE, resolution = 2) %>%
  RunUMAP(reduction = "harmony", dims = 2:30)

add_filename4 <- '_137_samples_v3.0_cluster_peaks_harmony_data.type_MM'
p <- DimPlot(panc4, group.by = 'Cancer', pt.size = 0.1)
pdf(paste0("Dimplot_cancer_reclustered_immune_snATAC_Merged_", add_filename4, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p <- DimPlot(panc4, label = T, pt.size = 0.1)
pdf(paste0("Dimplot_clusters_reclustered_immune_snATAC_Merged_", add_filename4, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p3 <- DimPlot(panc4,group.by = 'data.type', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_data_type_reclustered_immune_snATAC_Merged",add_filename4, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

combo.meta <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_immune_44_samples.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)
panc4 <- AddMetaData(panc4, combo.meta[,c('Cell_type_state', 'Cell_type_markers')])
p3 <- DimPlot(panc4,group.by = 'Cell_type_state', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_Cell_type_state_reclustered_immune_snATAC_Merged", add_filename4, ".pdf"),height=10,width=20, useDingbats = F)
print(p3)
dev.off()

p3 <- DimPlot(panc4,group.by = 'Cell_type_markers', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_Cell_type_markers_reclustered_immune_snATAC_Merged", add_filename4, ".pdf"),height=10,width=25, useDingbats = F)
print(p3)
dev.off()

p3 <- DimPlot(panc4,group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_cell_type.harmonized.cancer_reclustered_immune_snATAC_Merged", add_filename4, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

DepthCor(panc4, reduction = 'harmony')
ggsave(paste0("Corplot_", add_filename4, ".pdf"),height=5,width=6, useDingbats = F)

library(cowplot)
library(data.table)
library(ggplot2)
all.meta <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/data.freeze.v3.0/merged_by_cancer_50PCs/All_137_samples_metadata_data_freeze_v3.1.tsv')


ggplot(data = all.meta, aes(color=data.type, x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, shape = 16) +
  scale_color_brewer(palette = 'Paired')+
  theme_cowplot() +
  facet_wrap(~Cancer, ncol=5)
ggsave('Dimplot_all_v3.0_data.type.pdf',height=10,width=10, useDingbats=FALSE)

ggplot(data = all.meta, aes(color=cell_type.harmonized.cancer, x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, shape = 16) +
  scale_color_manual(values = colors$cell_type)+
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=6)))+
  facet_wrap(~Cancer, ncol=5)
ggsave('Dimplot_all_v3.0_data.type.pdf',height=10,width=10, useDingbats=FALSE)


saveRDS(panc2, 'Reclustered_immune_snATAC_Merged_137_samples_v3.0_cluster_peaks_harmony_MM.rds')
saveRDS(panc3, 'Reclustered_immune_snATAC_Merged_137_samples_v3.0_cluster_peaks_harmony_data.type.rds')
saveRDS(panc4, 'Reclustered_immune_snATAC_Merged_137_samples_v3.0_cluster_peaks_harmony_MM_data.type.rds')

colors <- readRDS('~/R_working_dir/scripts/snATAC/Colors_panatac_v1.0.rds')





panc$Barcodes_cancer <- paste(panc$Cancer, panc$Piece_ID, panc$Barcodes, sep = '_')
panc <- RenameCells(object = panc, new.names = panc$Barcodes_cancer)
panc <- AddMetaData(panc, meta)
add_filename4 <- '_137_samples_v3.0_cluster_peaks_harmony_data.type_MM'
p3 <- DimPlot(panc,group.by = 'Cell_type_state', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_Cell_type_state_reclustered_immune_snATAC_Merged", add_filename4, ".pdf"),height=10,width=12, useDingbats = F)
print(p3)
dev.off()

add_filename4 <- '_137_samples_v3.0_cluster_peaks_harmony_data.type_MM'
p3 <- DimPlot(panc,group.by = 'Cell_type_markers', pt.size = 0.1,label=T)
pdf(paste0("Dimplot_Cell_type_markers_reclustered_immune_snATAC_Merged", add_filename4, ".pdf"),height=10,width=20, useDingbats = F)
print(p3)
dev.off()
