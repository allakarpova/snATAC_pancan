# Recluster normal cells for better annotaion
###libraries
##################

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(future))


meta <- fread('/diskmnt/Projects/snATAC_primary/04_celltyped_rds/cell_type_snRNA_merged/v5.0_data_freeze/All_RNA_cell_type.harmonized.cancer.meta.data')%>% data.frame(row.names = 1)
meta <- meta %>% mutate(Cancer.subtype = case_when(grepl('HT029B1|HT035B1|1408|HT141B1|HT271B1|HT268B1', Piece_ID) ~'Basal',
                                                   grepl('HT214|HT265', Piece_ID) ~ 'Her2-enriched',
                                                   Cancer == 'BRCA' ~ 'Luminal',
                                                   TRUE ~ 'NA'))
meta$Cancer.new <- case_when(meta$Cancer.subtype=='Basal' ~ 'BRCA_Basal',
                             TRUE ~ meta$Cancer)


cat('PDAC \n')
cat('opening object \n')
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snRNA/Merged_objects/PDAC/PDAC_merge_obj_processed-with_doublets.rds')
cat('done \n')
panc <- AddMetaData(panc, meta)

colnames(panc@meta.data)

panc <-  subset(panc, cell_type.harmonized.cancer %in% c('Acinar', 'Normal epithelial cells'))
dim(panc)

panc <- panc %>% SCTransform(
  assay = 'RNA',
  vars.to.regress = c("nCount_RNA", "percent.mito"),
  conserve.memory = T,
  return.only.var.genes = T
) %>%
  RunPCA(assay = 'SCT', do.print = FALSE) %>%
  RunUMAP(dims = 1:30, assay = 'SCT') %>%
  NormalizeData(assay = 'RNA') %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(verbose = T, resolution = 2)

saveRDS(panc, '/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/PDAC_RNA_obj_normal_epith_acinar_v5.0_data_freeze.rds')

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/Dimplot_clusters_normal_cells_PDAC.pdf"), height=7,width=8,useDingbats=FALSE)

DimPlot(panc, group.by = 'Piece_ID', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/Dimplot_Piece_ID_normal_cells_PDAC.pdf"), height=7,width=9,useDingbats=FALSE)

#DimPlot(panc, group.by = 'cell_type', label = T)
#ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/Dimplot_cell_type_normal_cells_PDAC.pdf"), height=7,width=10,useDingbats=FALSE)

#####################################
cat('ccRCC \n')
cat('opening object \n')
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snRNA/Merged_objects/ccRCC/ccRCC_merge_obj_processed-with_doublets.rds')
cat('done \n')
panc <- AddMetaData(panc, meta)

panc <-  subset(panc, cell_type.harmonized.cancer %in% c('Normal epithelial cells'))

panc <- panc %>% SCTransform(
  assay = 'RNA',
  vars.to.regress = c("nCount_RNA", "percent.mito"),
  conserve.memory = T,
  return.only.var.genes = T
) %>%
  RunPCA(assay = 'SCT', do.print = FALSE) %>%
  RunUMAP(dims = 1:30, assay = 'SCT') %>%
  NormalizeData(assay = 'RNA') %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(verbose = T, resolution = 2)

saveRDS(panc, '/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/ccRCC/ccRCC_RNA_obj_normal_epith_v5.0_data_freeze.rds')

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/ccRCC/Dimplot_clusters_normal_cells_ccRCC.pdf"), height=7,width=8,useDingbats=FALSE)

DimPlot(panc, group.by = 'Piece_ID', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed//ccRCC/Dimplot_Piece_ID_normal_cells_ccRCC.pdf"), height=7,width=10,useDingbats=FALSE)

#DimPlot(panc, group.by = 'cell_type', label = T)
#ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed//ccRCC/Dimplot_cell_type_normal_cells_ccRCC.pdf"), height=7,width=10,useDingbats=FALSE)

#####################################
cat('CRC \n')
cat('opening object \n')
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snRNA/Merged_objects/CRC/CRC_merge_obj_processed-with_doublets.rds')
cat('done \n')
panc <- AddMetaData(panc, meta)

panc <-  subset(panc, cell_type.harmonized.cancer %in% c('Normal epithelial cells'))

panc <- panc %>% SCTransform(
  assay = 'RNA',
  vars.to.regress = c("nCount_RNA", "percent.mito"),
  conserve.memory = T,
  return.only.var.genes = T
) %>%
  RunPCA(assay = 'SCT', do.print = FALSE) %>%
  RunUMAP(dims = 1:30, assay = 'SCT') %>%
  NormalizeData(assay = 'RNA') %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(verbose = T, resolution = 2)

saveRDS(panc, '/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/CRC_RNA_obj_normal_epith_v5.0_data_freeze.rds')

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/Dimplot_clusters_normal_cells_CRC.pdf"), height=7,width=8,useDingbats=FALSE)

DimPlot(panc, group.by = 'Piece_ID', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/Dimplot_Piece_ID_normal_cells_CRC.pdf"), height=7,width=10,useDingbats=FALSE)

#DimPlot(panc, group.by = 'cell_type', label = T)
#ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/Dimplot_cell_type_normal_cells_CRC.pdf"), height=7,width=10,useDingbats=FALSE)


