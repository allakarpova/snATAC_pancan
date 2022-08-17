library(Signac)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(doParallel)
library(future)
library(RColorBrewer)

dir.create('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v4.0_data_freeze')
setwd('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v4.0_data_freeze')
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

objects <- c('BRCA_snATAC_Merged.PancanSet.noDoublets.20211022.rds', 'ccRCC_snATAC_Merged.PancanSet.noDoublets.20211022.rds', 'CESC_snATAC_Merged.PancanSet.noDoublets.20211022.rds',
  'CRC_snATAC_Merged.PancanSet.noDoublets.20211022.rds', 'GBM_snATAC_Merged.PancanSet.noDoublets.20211022.rds', 'HNSCC_snATAC_Merged.PancanSet.noDoublets.20211022.rds',
  'MM_snATAC_Merged.PancanSet.noDoublets.20211022.rds', 'OV_snATAC_Merged.PancanSet.noDoublets.20211022.rds', 'PDAC_snATAC_Merged.PancanSet.noDoublets.20211022.rds',
  'UCEC_snATAC_Merged.PancanSet.noDoublets.20211022.rds')
out.metas <- c( 'BRCA_19_samples_metadata_v6_data_freeze_v4.1.tsv','ccRCC_29_samples_metadata_v6_data_freeze_v4.1.tsv','CESC_5_samples_metadata_v6_data_freeze_v4.1.tsv',
               'CRC_22_samples_metadata_v6_data_freeze_v4.1.tsv', 'GBM_22_samples_metadata_v6_data_freeze_v4.1.tsv', 'HNSCC_7_samples_metadata_v6_data_freeze_v4.1.tsv',
               'MM_14_samples_metadata_v6_data_freeze_v4.1.tsv', 'OV_8_samples_metadata_v6_data_freeze_v4.1.tsv', 'PDAC_12_samples_metadata_v6_data_freeze_v4.1.tsv',
               'UCEC_11_samples_metadata_v6_data_freeze_v4.1.tsv')

meta.to.upd.name <- '/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v6.0/output/All_149_samples_metadata_data_freeze_v4.0.tsv'
previous.meta.name <- '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs_correct_peaks/All_137_samples_metadata_data_freeze_v3.0.tsv'
meta.to.upd <- fread(meta.to.upd.name)
meta.to.upd$Barcodes_cancer <- meta.to.upd$V1
#colnames(meta.to.upd)[colnames(meta.to.upd)=='V1'] <- 'Barcodes'
previous.meta <- fread(previous.meta.name)
previous.meta <- previous.meta %>% 
  select(Cancer, Barcodes_cancer, cell_type.harmonized.cancer) 
#colnames(previous.meta)[1] <- 'Barcodes_merged'

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/Colors_panatac_v1.0.rds')
colors$cell_type <- c(colors$cell_type, c(Unknown = 'grey', 'Cholangiocytes' = '#ffd92f', 'Astrocytes' = '#e78ac3', Microglia= '#8da0cb'))
saveRDS(colors, '/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/Colors_panatac_v1.0.rds')
walk2(objects, out.metas, function(input.obj, output.meta.name) {
  tryCatch( {
    print(input.obj)
    input.obj <- 'BRCA_snATAC_Merged.PancanSet.noDoublets.20211022.rds'
    output.meta.name <- 'BRCA_19_samples_metadata_v6_data_freeze_v4.1.tsv'
    
    output.obj <- paste0(str_split_fixed(input.obj, '2021', 2)[,1], '20211022.rds')
    cancer.type <- str_split_fixed(input.obj, '_', 2)[1]

    cat('opening object\n')
    panc <- readRDS(paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/', input.obj))
    
    #sample_barcodes <- meta.to.upd[[1]]
    cells.in.object <- colnames(panc)
    if (cancer.type=='PDAC') {
      previous.meta.cancer <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/PDAC_12_samples_metadata_data_freeze_v3.0.tsv') %>%
        select(Barcodes_cancer, cell_type.harmonized.cancer) 
    } else if (cancer.type=='UCEC'){
      previous.meta.cancer <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs_correct_peaks/UCEC_12_samples_metadata_data_freeze_v3.0.tsv') %>%
        select(-Barcodes_cancer) %>%
        dplyr::rename(Barcodes_cancer=V1) %>%
        select(Barcodes_cancer, cell_type.harmonized.cancer) 
    } else {
      previous.meta.cancer <- previous.meta %>%
        filter(Cancer == cancer.type) %>% 
        #data.frame(row.names = 'Barcodes_merged') %>%
        select(-Cancer) 
    }
    head(previous.meta.cancer)
    
    meta.to.upd.cancer <- meta.to.upd %>%
      filter(Cancer == cancer.type) %>%
      select(-seurat_clusters)
    
    meta.merged <- left_join(meta.to.upd.cancer, 
                             previous.meta.cancer,
                             by = 'Barcodes_cancer')
    if (cancer.type=='ccRCC') {
      meta.merged$Barcodes_merged <- meta.merged$Barcodes_cancer
      meta.merged$Barcodes_merged <- case_when(grepl('K1900070', meta.merged$Barcodes_merged) ~ str_replace(string = meta.merged$Barcodes_merged, pattern = 'ccRCC', replacement = 'PKD'),
                                               TRUE ~ meta.merged$Barcodes_merged)
    } else {
      meta.merged$Barcodes_merged <- meta.merged$Barcodes_cancer
    }
    
    
    
    meta.merged <- meta.merged %>% 
      data.frame(row.names = 'Barcodes_merged')
    
    x <- meta.merged$cell_type.harmonized.cancer=='B-cells/Plasma' & !is.na(meta.merged$cell_type.harmonized.cancer=='B-cells/Plasma')
    meta.merged$cell_type.harmonized.cancer[x] <- meta.merged$cell_type.harmonized[x]
    
    panc@meta.data <- panc@meta.data[,c('blacklist_ratio', "orig.ident",'nCount_pancan','nFeature_pancan', 'seurat_clusters', 'pancan_snn_res.0.8', 
                                        'pancan_snn_res.2')]
    
    cat('adding meta.merged\n')
    panc <- AddMetaData(panc, meta.merged)
    
    total_fragments_cell <- panc$passed_filters
    peak.counts <- colSums(x = GetAssayData(panc, slot='counts'))
    frip <- peak.counts *100 / total_fragments_cell
    panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_pancan')
    panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan')
    
    
    
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
    
    if (cancer.type %in% c('HNSCC', 'OV', 'GBM', 'UCEC')) {
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
    
    if (cancer.type=='CESC') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2==22] <- 'Plasma'
    }
    
    if (cancer.type=='CRC') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(47,46,50)] <- 'Unknown'
    }
    
    if (cancer.type=='GBM') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(38)] <- 'Astrocytes'
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(54,52)] <- 'Unknown'
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(59)] <- 'Tumor'
    }
    
    if (cancer.type=='MM') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(32,33)] <- 'Unknown'
    }
    if (cancer.type=='UCEC') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(29)] <- 'Plasma'
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(10)] <- 'Macrophages'
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(21)] <- 'T-cells'
      
    }
    if (cancer.type=='PDAC') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(35)] <- 'Unknown'
    }
    
    p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = colors$cell_type[names(colors$cell_type) %in% (panc$cell_type.harmonized.cancer %>% unique)])
    ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer_new.pdf"), plot = p4, height=12,width=15, useDingbats = F)
    
    p1 <- DimPlot(panc, split.by = 'Piece_ID', group.by = 'cell_type.harmonized.cancer',  pt.size =ps,label = T, ncol = 5)
    ggsave(paste0(cancer.type,"_snATAC_Merged_Piece_ID_split2.pdf"),plot = p1,height=40,width=40, useDingbats = F)
    
    panc$cell_type.harmonized.cancer %>% unique
    # Add annotation
    #Annotation(panc) <- annotations
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


ccrcc <- fread('ccRCC_29_samples_metadata_v6_data_freeze_v4.1.tsv')
ccrcc$cell_type.harmonized.cancer <- case_when(ccrcc$pancan_snn_res.2 %in% c(3, 57) ~ 'Unknown',
                                               TRUE ~ ccrcc$cell_type.harmonized.cancer)
head(ccrcc, n=1)
my.col <- colors$cell_type[names(colors$cell_type) %in% (ccrcc$cell_type.harmonized.cancer  %>% unique)]
ggplot(data=ccrcc[,c('UMAP_1', 'UMAP_2', 'cell_type.harmonized.cancer')], aes(x = UMAP_1, y=UMAP_2, color = cell_type.harmonized.cancer)) +
  geom_point(size = 0.5, shape=16) +
  scale_color_manual(values = my.col) +
  cowplot::theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=6)))
ggsave('ccRCC_snATAC_Merged_cell_type.harmonized.cancer_new2.pdf', height=20,width=20, useDingbats = F)


fwrite(ccrcc, 'ccRCC_29_samples_metadata_v6_data_freeze_v4.2.tsv', sep = '\t')

####
ucec <- fread('ccRCC_29_samples_metadata_v6_data_freeze_v4.1.tsv')
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/Colors_panatac_v1.0.rds')
ucec$cell_type.harmonized.cancer <- case_when(ucec$pancan_snn_res.2 %in% c(27) ~ 'Secretory Endometrial epithelial cells',
                                              ucec$pancan_snn_res.2 %in% c(33) ~ 'Endothelial',
                                               TRUE ~ ucec$cell_type.harmonized.cancer)
head(ucec, n=1)
my.col <- colors$cell_type[names(colors$cell_type) %in% (ucec$cell_type.harmonized.cancer  %>% unique)]
ggplot(data=ucec[,c('UMAP_1', 'UMAP_2', 'cell_type.harmonized.cancer')], aes(x = UMAP_1, y=UMAP_2, color = cell_type.harmonized.cancer)) +
  geom_point(size = 0.5, shape=16) +
  scale_color_manual(values = my.col) +
  cowplot::theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=6)))
ggsave('ccRCC_snATAC_Merged_cell_type.harmonized.cancer_new2.pdf', height=20,width=20, useDingbats = F)


fwrite(ucec, 'ccRCC_29_samples_metadata_v6_data_freeze_v4.2.tsv', sep = '\t')


colors$cell_type <- c(colors$cell_type, c('Secretory Endometrial epithelial cells' = '#a6761d', 'Ciliated Endometrial epithelial cells' = '#e6ab02'))
saveRDS(colors, '/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/Colors_panatac_v1.0.rds')

####
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/Colors_panatac_v1.0.rds')

brca <- fread('BRCA_19_samples_metadata_v6_data_freeze_v4.1.tsv')
brca <- brca[,-11]
head(brca, n = 1)
brca <- brca %>% mutate (
  cell_type.harmonized.cancer = case_when(Cancer == 'BRCA' & pancan_snn_res.2 == 27 ~ 'Luminal mature',
                                              Cancer == 'BRCA' & pancan_snn_res.2== 32 ~ 'Luminal progenitor',
                                              Cancer == 'BRCA' & pancan_snn_res.2== 18 ~ 'Basal progenitor',
                                              TRUE ~ cell_type.harmonized.cancer))
colors$cell_type <- c(colors$cell_type, 
                      c( 'Luminal mature' = '#005824', 'Luminal progenitor' = '#41ae76', 'Basal progenitor' = '#99d8c9'))
saveRDS(colors, '/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/Colors_panatac_v1.0.rds')
my.col <- colors$cell_type[names(colors$cell_type) %in% (brca$cell_type.harmonized.cancer  %>% unique)]
ggplot(data=brca[,c('UMAP_1', 'UMAP_2', 'cell_type.harmonized.cancer')], aes(x = UMAP_1, y=UMAP_2, color = cell_type.harmonized.cancer)) +
  geom_point(size = 0.5, shape=16) +
  scale_color_manual(values = my.col) +
  cowplot::theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=6)))
ggsave('BRCA_snATAC_Merged_cell_type.harmonized.cancer_new2.pdf', height=20,width=20, useDingbats = F)

fwrite(brca, 'BRCA_19_samples_metadata_v6_data_freeze_v4.2.tsv', sep = '\t')


