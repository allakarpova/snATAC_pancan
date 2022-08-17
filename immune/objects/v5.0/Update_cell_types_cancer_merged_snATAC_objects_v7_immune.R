library(Signac)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(doParallel)
library(future)
library(RColorBrewer)

dir.create('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze_for_immune')
setwd('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze_for_immune')
#annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

objects <- c('CESC_9_snATAC_Merged_new_peaks_normalized_v5_samples.rds', 
             'PBMC_3_snATAC_Merged_new_peaks_normalized_v5_samples.rds',
             'SKCM_4_snATAC_Merged_new_peaks_normalized_v5_samples.rds'
             )
out.metas <- c('CESC_9_samples_metadata_v7_data_freeze_v5.0.tsv',
                'PBMC_3_samples_metadata_v7_data_freeze_v5.0.tsv', 
                'SKCM_4_samples_metadata_v7_data_freeze_v5.0.tsv')

meta.to.upd.name <- '/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v7.0_for_immune/output/All_164_samples_metadata_data_freeze_v5.0.tsv'
previous.meta.name <- '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/All_159_samples_metadata_data_freeze_v5.0.tsv'
meta.to.upd <- fread(meta.to.upd.name)
meta.to.upd$Barcodes_cancer <- meta.to.upd$V1
#colnames(meta.to.upd)[colnames(meta.to.upd)=='V1'] <- 'Barcodes'
previous.meta <- fread(previous.meta.name) %>% rename(Barcodes_cancer=V1)
previous.meta <- previous.meta %>% 
  select(Cancer, Barcodes_cancer, cell_type.harmonized.cancer) 
#colnames(previous.meta)[1] <- 'Barcodes_merged'

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/Colors_panatac_v1.0.rds')
colors$cell_type <- c(colors$cell_type, c( 'Monocytes' = '#b5179e', 'MAIT' = '#4cc9f0'))
new.cancer.colors <- c(colorRampPalette(colors = c("#006ba6","#0496ff","#ffbc42","#d81159","#8f2d56"))(6),
                       colorRampPalette(colors = c("#a8d5e2","#f9a620","#ffd449","#548c2f","#104911"))(6))

names(new.cancer.colors) <- c('BRCA', 'ccRCC', 'CESC',  'CRC', 'GBM', 'HNSCC', 'MM', 'OV', 'PBMC', 'PDAC','SKCM', 'UCEC')
colors$Cancer <- new.cancer.colors
saveRDS(colors, '/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v5.0/Colors_panatac_v1.0.rds')

walk2(objects, out.metas, function(input.obj, output.meta.name) {
  tryCatch( {
    print(input.obj)
    input.obj <- 'CESC_9_snATAC_Merged_new_peaks_normalized_v5_samples.rds'
    output.meta.name <- 'CESC_9_samples_metadata_v7_data_freeze_v5.0.tsv'
    
    cancer.type <- str_split_fixed(input.obj, '_', 2)[1]
    
    cat('opening object\n')
    panc <- readRDS(paste0('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v5.0/snATAC/Merged_objects_intermediate_with_tumor/cancer_level/', input.obj))
  
    cells.in.object <- colnames(panc)
    previous.meta.cancer <- previous.meta %>%
      filter(Cancer == cancer.type) %>% 
      select(-Cancer) 
    
    head(previous.meta.cancer)
    
    meta.to.upd.cancer <- meta.to.upd %>%
      filter(Cancer == cancer.type) %>%
      select(-seurat_clusters)
    
    meta.merged <- left_join(meta.to.upd.cancer, 
                             previous.meta.cancer,
                             by = 'Barcodes_cancer')
    
    meta.merged <- meta.merged %>% 
      data.frame(row.names = 'Barcodes_cancer') %>% 
      select(-V1)
    
    head(meta.merged)
      
    panc@meta.data <- panc@meta.data[,c('blacklist_ratio', "orig.ident",'nCount_pancan','nFeature_pancan', 'seurat_clusters', 
                                        'pancan_snn_res.2', 'doublet_score_rna' ,'predicted_doublet_rna' , 'doublet_score_atac', 'predicted_doublet_atac', 'Doublet_final')]
    
    cat('adding meta.merged\n')
    panc <- AddMetaData(panc, meta.merged)
    
    total_fragments_cell <- panc$passed_filters
    peak.counts <- colSums(x = GetAssayData(panc, slot='counts'))
    frip <- peak.counts *100 / total_fragments_cell
    panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_pancan')
    panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan')
    
    head(panc@meta.data)
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
    
    p2 <- DimPlot(panc, group.by = 'pancan_snn_res.2', pt.size = ps,label=T)
    ggsave(paste0(cancer.type,"_snATAC_Merged_pancan_snn_res.2.pdf"), plot = p2, height=12,width=14, useDingbats = F)
    
    p3 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized', cols = col_vector)
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized_with_doublets.pdf"), plot = p3, height=12,width=15, useDingbats = F)
    
    p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = col_vector)
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer_old.pdf"), plot = p4, height=12,width=15, useDingbats = F)
    
    p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'data.type')
    ggsave(paste0(cancer.type,"_snATAC_Merged_data.type.pdf"), plot = p4, height=12,width=13, useDingbats = F)
    
    if (cancer.type %in% c('PBMC', 'SKCM')) {
      tb <- table(panc$pancan_snn_res.2, panc$cell_type.harmonized)
      cluster.match.celltype <- apply(tb, 1, function(x) {
        colnames(tb)[which.max (x)]
      })
      panc$cell_type.harmonized.cancer <- cluster.match.celltype[as.character(panc$pancan_snn_res.2)]
    } else {
      tb <- table(panc$pancan_snn_res.2, panc$cell_type.harmonized.cancer)
      cluster.match.celltype <- apply(tb, 1, function(x) {
        colnames(tb)[which.max (x)]
      })
      panc$cell_type.harmonized.cancer <- cluster.match.celltype[as.character(panc$pancan_snn_res.2)]
    }
    
    
    p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = colors$cell_type[names(colors$cell_type) %in% (panc$cell_type.harmonized.cancer %>% unique)])
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer_new.pdf"), plot = p4, height=12,width=15, useDingbats = F)
    
    p1 <- DimPlot(panc, split.by = 'Piece_ID', group.by = 'cell_type.harmonized.cancer',  pt.size =ps,label = T, ncol = 5)
    ggsave(paste0(cancer.type,"_snATAC_Merged_Piece_ID_split2.pdf"),plot = p1,height=40,width=40, useDingbats = F)
    
    panc$cell_type.harmonized.cancer %>% unique
    # Add annotation
    #Annotation(panc) <- annotations
    towrite <- cbind(Embeddings(panc, reduction = "umap"), panc@meta.data)
    fwrite(towrite, output.meta.name, row.names = T, sep = '\t')
    cat('saving\n')
    #saveRDS(panc, output.obj)
  }, 
  error=function(cond) {
    message(paste("cancer type failed", input.obj))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  })
})


