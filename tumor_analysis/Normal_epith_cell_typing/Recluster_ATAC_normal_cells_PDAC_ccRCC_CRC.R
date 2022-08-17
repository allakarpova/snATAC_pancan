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

annotation <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

meta <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/All_159_samples_metadata_data_freeze_v5.0.tsv')%>% data.frame(row.names = 1)
meta <- meta %>% mutate(Cancer.subtype = case_when(grepl('HT029B1|HT035B1|1408|HT141B1|HT271B1|HT268B1', Piece_ID) ~'Basal',
                                                   grepl('HT214|HT265', Piece_ID) ~ 'Her2-enriched',
                                                   Cancer == 'BRCA' ~ 'Luminal',
                                                   TRUE ~ 'NA'))
meta$Cancer.new <- case_when(meta$Cancer.subtype=='Basal' ~ 'BRCA_Basal',
                             TRUE ~ meta$Cancer)


cat('PDAC \n')
cat('opening object \n')
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/PDAC_snATAC_Merged.PancanSet.noDoublets.20220209.rds')
cat('done \n')
Annotation(panc) <- annotation
panc <- AddMetaData(panc, meta)

colnames(panc@meta.data)

panc <-  subset(panc, cell_type.harmonized.cancer %in% c('Acinar', 'Normal epithelial cells'))
dim(panc)

panc <- panc %>% 
  RunTFIDF() %>%
  FindTopFeatures( min.cutoff = 20) %>%
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400
  ) %>%
  RunUMAP(dims = 2:30,reduction = 'lsi') %>%
  FindNeighbors(dims = 2:30, reduction = 'lsi') %>%
  FindClusters(verbose = T, resolution = 0.5, algorithm = 3)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities <- GeneActivity(panc)
panc[['RNA']] <- CreateAssayObject(counts = gene.activities)
panc <- NormalizeData(
  object = panc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(panc$nCount_RNA)
)

saveRDS(panc, '/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/ATAC/PDAC_ATAC_obj_normal_epith_acinar_v5.0_data_freeze.rds')

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/ATAC/Dimplot_clusters_normal_cells_PDAC.pdf"), height=7,width=8,useDingbats=FALSE)

DimPlot(panc, group.by = 'Piece_ID', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/ATAC/Dimplot_Piece_ID_normal_cells_PDAC.pdf"), height=7,width=9,useDingbats=FALSE)

#DimPlot(panc, group.by = 'cell_type', label = T)
#ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/PDAC/Dimplot_cell_type_normal_cells_PDAC.pdf"), height=7,width=10,useDingbats=FALSE)

#####################################
cat('ccRCC \n')
cat('opening object \n')
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/ccRCC_snATAC_Merged.PancanSet.noDoublets.20220209.rds')
cat('done \n')
Annotation(panc) <- annotation
panc <- AddMetaData(panc, meta)

panc <-  subset(panc, cell_type.harmonized.cancer %in% c('Normal epithelial cells'))

panc <- panc %>% 
  RunTFIDF() %>%
  FindTopFeatures( min.cutoff = 20) %>%
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400
  ) %>%
  RunUMAP(dims = 2:30,reduction = 'lsi') %>%
  FindNeighbors(dims = 2:30, reduction = 'lsi') %>%
  FindClusters(verbose = T, resolution = 0.5, algorithm = 3)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities <- GeneActivity(panc)
panc[['RNA']] <- CreateAssayObject(counts = gene.activities)
panc <- NormalizeData(
  object = panc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(panc$nCount_RNA)
)

saveRDS(panc, '/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/ccRCC/ATAC/ccRCC_ATAC_obj_normal_epith_v5.0_data_freeze.rds')

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/ccRCC/ATAC/Dimplot_clusters_normal_cells_ccRCC.pdf"), height=7,width=8,useDingbats=FALSE)

DimPlot(panc, group.by = 'Piece_ID', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed//ccRCC/ATAC/Dimplot_Piece_ID_normal_cells_ccRCC.pdf"), height=7,width=10,useDingbats=FALSE)

#DimPlot(panc, group.by = 'cell_type', label = T)
#ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed//ccRCC/Dimplot_cell_type_normal_cells_ccRCC.pdf"), height=7,width=10,useDingbats=FALSE)

#####################################
cat('CRC \n')
cat('opening object \n')
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/CRC_snATAC_Merged.PancanSet.noDoublets.20220209.rds')
cat('done \n')
Annotation(panc) <- annotation
panc <- AddMetaData(panc, meta)

panc <-  subset(panc, cell_type.harmonized.cancer %in% c('Normal epithelial cells'))

panc <- panc %>% 
  RunTFIDF() %>%
  FindTopFeatures( min.cutoff = 20) %>%
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400
  ) %>%
  RunUMAP(dims = 2:30,reduction = 'lsi') %>%
  FindNeighbors(dims = 2:30, reduction = 'lsi') %>%
  FindClusters(verbose = T, resolution = 0.5, algorithm = 3)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities <- GeneActivity(panc)
panc[['RNA']] <- CreateAssayObject(counts = gene.activities)
panc <- NormalizeData(
  object = panc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(panc$nCount_RNA)
)

saveRDS(panc, '/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/ATAC/CRC_ATAC_obj_normal_epith_v5.0_data_freeze.rds')

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/ATAC/Dimplot_clusters_normal_cells_CRC.pdf"), height=7,width=8,useDingbats=FALSE)

DimPlot(panc, group.by = 'Piece_ID', label = T)
ggsave(paste0( "/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/ATAC/Dimplot_Piece_ID_normal_cells_CRC.pdf"), height=7,width=10,useDingbats=FALSE)

#DimPlot(panc, group.by = 'cell_type', label = T)
#ggsave(paste0("/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/CRC/Dimplot_cell_type_normal_cells_CRC.pdf"), height=7,width=10,useDingbats=FALSE)


