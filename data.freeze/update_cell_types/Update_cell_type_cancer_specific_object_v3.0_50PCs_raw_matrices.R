library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(tidyverse)
library(data.table)
library(doParallel)
library(future)
library(RColorBrewer)

dir.create('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs_correct_peaks')
setwd('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs_correct_peaks')
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

objects <- c('BRCA_snATAC_Merged.PancanSet.noDoublets.20210929.rds', 'ccRCC_snATAC_Merged.PancanSet.noDoublets.20210929.rds', 'CESC_snATAC_Merged.PancanSet.noDoublets.20210929.rds',
  'CRC_snATAC_Merged.PancanSet.noDoublets.20210929.rds', 'GBM_snATAC_Merged.PancanSet.noDoublets.20210929.rds', 'HNSCC_snATAC_Merged.PancanSet.noDoublets.20210929.rds',
  'MM_snATAC_Merged.PancanSet.noDoublets.20210929.rds', 'OV_snATAC_Merged.PancanSet.noDoublets.20210929.rds', 'PDAC_snATAC_Merged.PancanSet.noDoublets.20210929.rds',
  'UCEC_snATAC_Merged.PancanSet.noDoublets.20210929.rds')
out.metas <- c( 'BRCA_19_samples_metadata_data_freeze_v3.0.tsv','ccRCC_29_samples_metadata_data_freeze_v3.0.tsv','CESC_5_samples_metadata_data_freeze_v3.0.tsv',
               'CRC_20_samples_metadata_data_freeze_v3.0.tsv', 'GBM_20_samples_metadata_data_freeze_v3.0.tsv', 'HNSCC_2_samples_metadata_data_freeze_v3.0.tsv',
               'MM_14_samples_metadata_data_freeze_v3.0.tsv', 'OV_4_samples_metadata_data_freeze_v3.0.tsv', 'PDAC_12_samples_metadata_data_freeze_v3.2.tsv',
               'UCEC_12_samples_metadata_data_freeze_v3.2.tsv')

meta.to.upd.name <- '/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v5.0/All_137_samples_metadata_data_freeze_v3.0.tsv'
previous.meta.name <- '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/All_137_samples_metadata_data_freeze_v3.1.tsv'
meta.to.upd <- fread(meta.to.upd.name)
meta.to.upd$Barcodes_cancer <- paste(meta.to.upd$Cancer, meta.to.upd$V1, sep='_')
colnames(meta.to.upd)[colnames(meta.to.upd)=='V1'] <- 'Barcodes'
previous.meta <- fread(previous.meta.name)
previous.meta <- previous.meta %>% 
  select(Cancer, Barcodes_cancer, cell_type.harmonized.cancer) 
#colnames(previous.meta)[1] <- 'Barcodes_merged'

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/Cell_type_colors_panatac_v1.0.rds')

walk2(objects, out.metas, function(input.obj, output.meta.name) {
  tryCatch( {
    print(input.obj)
    #input.obj <- 'BRCA_snATAC_Merged.PancanSet.noDoublets.20210929.rds'
    #output.meta.name <- 'BRCA_19_samples_metadata_data_freeze_v3.0.tsv'
    
    output.obj <- paste0(str_split_fixed(input.obj, '2021', 2)[,1], '20211001.rds')
    cancer.type <- str_split_fixed(input.obj, '_', 2)[1]

    cat('opening object\n')
    panc <- readRDS(paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/', input.obj))
    
    #sample_barcodes <- meta.to.upd[[1]]
    cells.in.object <- colnames(panc)
    previous.meta.cancer <- previous.meta %>%
      filter(Cancer == cancer.type) %>% 
      #data.frame(row.names = 'Barcodes_merged') %>%
      select(-Cancer) 
    
    meta.to.upd.cancer <- meta.to.upd %>%
      filter(Cancer == cancer.type) %>%
      select(-seurat_clusters)
    
    meta.merged <- left_join(meta.to.upd.cancer, 
                             previous.meta.cancer,
                             by = 'Barcodes_cancer')
    
    if (cancer.type=='BRCA') {
      obj.meta <- panc@meta.data %>% 
        mutate(Barcodes_merged = rownames(.)) %>% 
        mutate (Barcodes.change = gsub('_', '-', Barcodes_merged),
                Piece_ID = str_split_fixed(string = Piece_ID, pattern = '_', 2)[,2],
                Piece_ID = gsub('_', '-', Piece_ID)) %>%
        mutate (Barcodes.change = str_split_fixed(Barcodes.change, '-', 4)[,4]) %>%
        mutate(Barcodes_cancer = paste(cancer.type, Piece_ID, Barcodes.change, sep='_')) %>%
        select(Barcodes_merged, Barcodes_cancer)
      
     
      meta.merged <- right_join(meta.merged, obj.meta, by = 'Barcodes_cancer')
    } else {
      meta.merged$Barcodes_merged <- meta.merged$Barcodes_cancer
    }
    
    meta.merged <- meta.merged %>% 
      data.frame(row.names = 'Barcodes_merged')
    
    total_fragments_cell <- panc$passed_filters
    peak.counts <- colSums(x = GetAssayData(panc, slot='counts'))
    frip <- peak.counts *100 / total_fragments_cell
    panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_pancan')
    panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan')
    
    panc@meta.data <- panc@meta.data[,c('blacklist_ratio', "orig.ident",'nCount_pancan','nFeature_pancan', 'seurat_clusters', 'pancan_snn_res.0.8', 
                                        'pancan_snn_res.2', 'pct_read_in_peaks_pancan', 'peak_RF_pancan')]
    cat('adding meta.merged\n')
    panc <- AddMetaData(panc, meta.merged)
    #panc@meta.data <- (panc@meta.data %>% select(-Barcodes))
    #cat('adding previous.meta.cancer\n')
    #panc <-AddMetaData(panc, (previous.meta.cancer %>% select(cell_type.harmonized.cancer)))
    
    n <- length(unique(panc$cell_type.harmonized))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cat('plotting\n')
    ps <- ifelse(ncol(panc)>50000, 0.1, 1)
    p3 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized', cols = col_vector)
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.pdf"), plot = p3, height=12,width=15, useDingbats = F)
    p1 <- DimPlot(panc, group.by = 'Piece_ID', pt.size =ps,label = T)
    ggsave(paste0(cancer.type,"_snATAC_Merged_Piece_ID.pdf"),plot = p1,height=12,width=13, useDingbats = F)
    
    p1 <- DimPlot(panc, split.by = 'Piece_ID', pt.size =ps,label = T, ncol = 5)
    ggsave(paste0(cancer.type,"_snATAC_Merged_Piece_ID_split.pdf"),plot = p1,height=40,width=40, useDingbats = F)
    
    p2 <- DimPlot(panc, group.by = 'seurat_clusters', pt.size = ps,label=T)
    ggsave(paste0(cancer.type,"_snATAC_Merged_clusters0.8.pdf"), plot = p2, height=12,width=14, useDingbats = F)
    
    p3 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized', cols = col_vector)
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized_with_doublets.pdf"), plot = p3, height=12,width=15, useDingbats = F)
    
    p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = col_vector)
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer_old.pdf"), plot = p4, height=12,width=15, useDingbats = F)
    
    tb <- table(panc$seurat_clusters, panc$cell_type.harmonized.cancer)
    cluster.match.celltype <- apply(tb, 1, function(x) {
      colnames(tb)[which.max (x)]
    })
    panc$cell_type.harmonized.cancer <- cluster.match.celltype[as.character(panc$seurat_clusters)]
    
    p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = colors[names(colors) %in% (panc$cell_type.harmonized.cancer %>% unique)])
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer_new.pdf"), plot = p4, height=12,width=15, useDingbats = F)
    
    # Add annotation
    Annotation(panc) <- annotations
    fwrite(cbind(Embeddings(panc, reduction = "umap"), panc@meta.data), output.meta.name, row.names = T, sep = '\t')
    cat('saving\n')
    saveRDS(panc, output.obj)
  }, 
  error=function(cond) {
    message(paste("cancer type failed", input.obj))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  })
})






############## END #################################

#check gigantic object
panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/137Samples_PanCan_merged_obj/137_snATAC_113K_peaks_diffPCs.motifsAdded.chromvar.20210826.rds.gz')
meta <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs_correct_peaks/All_137_samples_metadata_data_freeze_v3.0.tsv') %>% 
  data.frame(row.names = 1, check.names = F, check.rows = F) %>%
  dplyr::rename(seurat_clusters.cancer = seurat_clusters)
cancer.type <- 'Pancan'
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/Cell_type_colors_panatac_v1.0.rds')

panc <- FindNeighbors(
  object = panc,
  reduction = 'lsi',
  dims = 2:30
)

panc <- FindClusters( 
  object = panc,
  algorithm = 3,
  resolution = 0.8,
  verbose = FALSE
)

p2 <- DimPlot(panc, group.by = 'seurat_clusters', pt.size = ps,label=T)
ggsave(paste0(cancer.type,"_snATAC_Merged_clusters0.8.pdf"), plot = p2, height=12,width=14, useDingbats = F)

total_fragments_cell <- panc$passed_filters
peak.counts <- colSums(x = GetAssayData(panc, slot='counts'))
frip <- peak.counts *100 / total_fragments_cell
panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_pancan_s')
panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan_s')

panc@meta.data <- panc@meta.data[,c('nCount_pancan_s','nFeature_pancan_s', 'seurat_clusters', 'pct_read_in_peaks_pancan_s', 'peak_RF_pancan_s')]
panc <- AddMetaData(panc, meta)
p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = colors)
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer.pdf"), plot = p4, height=12,width=15, useDingbats = F)




if (cancer.type=='BRCA') {
  obj.meta <- panc@meta.data %>% 
    mutate(Barcodes_merged = rownames(.)) %>% 
    mutate (Barcodes.change = gsub('_', '-', Barcodes_merged),
            Piece_ID = gsub('_', '-', Piece_ID)) %>%
    mutate (Barcodes.change = str_split_fixed(Barcodes.change, '-', 4)[,4]) %>%
    mutate(Barcodes_cancer = paste(Piece_ID, Barcodes.change, sep='_')) %>%
    select(Barcodes_merged, Barcodes_cancer)
  head(obj.meta)
} else if (cancer.type %in% c('ccRCC', 'CESC', 'CRC', 'GBM')) {
  obj.meta <- panc@meta.data %>% 
    mutate(Barcodes_merged = rownames(.)) %>% 
    mutate (Barcodes_merged = str_split_fixed(Barcodes, '_', 3)[,3]) %>%
    mutate(Barcodes_cancer = paste(Piece_ID, Barcodes, sep='_')) %>%
    select(Barcodes_merged, Barcodes_cancer)
  setequal(obj.meta$Sample_barcodes, sample_barcodes)
}


if (cancer.type=='BRCA' | cancer.type=='ccRCC') {
  panc <- RenameCells(panc, new.names = rownames(meta.scrub))
}








