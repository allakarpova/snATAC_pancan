#recluster cells on chromvar motifs
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))


panc <- readRDS('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/v3.0/combo_only/Reclustered_immune_snATAC_Merged_44_combo_samples_v3.0_cluster_peaks.type.chromvar.rds')
combo.meta <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_immune_44_samples.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)

panc <- AddMetaData(panc, combo.meta[,c('Cell_type_state', 'Cell_type_markers', 'cell_type_general')])
panc$Cell_type_state <- case_when(is.na(panc$Cell_type_state) | grepl('oublet', panc$Cell_type_state) ~ 'Doublet',
                                        TRUE ~ panc$Cell_type_state)

panc$cell_type_general <- case_when(is.na(panc$Cell_type_state) | grepl('oublet', panc$Cell_type_state) ~ 'Doublet',
                                          grepl('CD4/CD8', panc$Cell_type_state) ~ 'CD4/CD8− T−cells',
                                          grepl('CD4', panc$Cell_type_state) ~ 'CD4 T-cells',
                                          grepl('CD8', panc$Cell_type_state) ~ 'CD8 T-cells',
                                          grepl('NK', panc$Cell_type_state) ~ 'NK cells',
                                          TRUE ~ panc$Cell_type_state)

DefaultAssay(panc) <- 'chromvar'
rownames(panc)
hist(FetchData(panc, vars='MA1581.1')[,1])
panc <- ScaleData(panc, assay = 'chromvar')
panc <- FindVariableFeatures(panc, assay = 'chromvar')

panc <- panc %>% 
  RunPCA(assay = 'chromvar', ) %>%
  FindNeighbors(reduction = "pca", dims = 1:40) %>%
  FindClusters(verbose = FALSE, resolution = 1.1) %>%
  RunUMAP(reduction = "pca", dims = 1:40, 
          reduction.name = "chromvar.umap", 
          reduction.key = "chromvarUMAP_")

DimPlot(panc, reduction = "chromvar.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) +
  DimPlot(panc, reduction = "chromvar.umap", group.by = "Cell_type_state", label = TRUE, label.size = 2.5, repel = TRUE)

panc <- panc %>%
  RunUMAP(reduction = "pca", dims = 1:40, 
          reduction.name = "chromvar40.umap", 
          reduction.key = "chromvar40UMAP_")
ElbowPlot(panc, reduction = 'pca', ndims = 50)

DimPlot(panc, reduction = "chromvar.umap", group.by = "cell_type_general", label = TRUE, label.size = 2.5) +
  DimPlot(panc, reduction = "umap", group.by = "cell_type_general", label = TRUE, label.size = 2.5)
ggsave('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/motif_clustering/combo_only/Dimplot_chromvar_ATAC_umaps_by_cell_type_general.pdf', width = 12, height = 5, useDingbats = F)

DimPlot(panc, reduction = "chromvar.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5) +
  DimPlot(panc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5)
ggsave('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/motif_clustering/combo_only/Dimplot_chromvar_ATAC_umaps_by_seurat_clusters.pdf', width = 12, height = 5, useDingbats = F)







