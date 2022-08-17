library(Signac)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(doParallel)
library(future)
library(RColorBrewer)

dir.create('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze')
setwd('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze')
#annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

objects <- c('BRCA_snATAC_Merged.PancanSet.noDoublets.20220209.rds', 'ccRCC_snATAC_Merged.PancanSet.noDoublets.20220209.rds', 'CESC_snATAC_Merged.PancanSet.noDoublets.20220209.rds',
             'CRC_snATAC_Merged.PancanSet.noDoublets.20220209.rds', 'GBM_snATAC_Merged.PancanSet.noDoublets.20220209.rds', 'HNSCC_snATAC_Merged.PancanSet.noDoublets.20220209.rds',
             'MM_snATAC_Merged.PancanSet.noDoublets.20220209.rds', 'OV_snATAC_Merged.PancanSet.noDoublets.20220209.rds', 'PDAC_snATAC_Merged.PancanSet.noDoublets.20220209.rds',
             'UCEC_snATAC_Merged.PancanSet.noDoublets.20220209.rds')
out.metas <- c( 'BRCA_19_samples_metadata_v7_data_freeze_v5.0.tsv','ccRCC_29_samples_metadata_v6_data_freeze_v5.0.tsv','CESC_11_samples_metadata_v6_data_freeze_v5.0.tsv',
                'CRC_22_samples_metadata_v6_data_freeze_v5.0.tsv', 'GBM_20_samples_metadata_v6_data_freeze_v5.0.tsv', 'HNSCC_11_samples_metadata_v6_data_freeze_v5.0.tsv',
                'MM_14_samples_metadata_v6_data_freeze_v5.0.tsv', 'OV_10_samples_metadata_v6_data_freeze_v5.0.tsv', 'PDAC_12_samples_metadata_v6_data_freeze_v5.0.tsv',
                'UCEC_11_samples_metadata_v6_data_freeze_v5.0.tsv')

meta.to.upd.name <- '/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v7.0/output/All_159_samples_metadata_data_freeze_v5.0.tsv'
previous.meta.name <- '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v4.0_data_freeze/All_149_samples_metadata_data_freeze_v4.1.tsv'
meta.to.upd <- fread(meta.to.upd.name)
meta.to.upd$Barcodes_cancer <- meta.to.upd$V1
#colnames(meta.to.upd)[colnames(meta.to.upd)=='V1'] <- 'Barcodes'
previous.meta <- fread(previous.meta.name)
previous.meta <- previous.meta %>% 
  select(Cancer, Barcodes_cancer, cell_type.harmonized.cancer) 
#colnames(previous.meta)[1] <- 'Barcodes_merged'

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v4.0/Colors_panatac_v1.0.rds')
colors$cell_type <- c(colors$cell_type, c(Unknown = 'grey', 'Cholangiocytes' = '#ffd92f', 'Astrocytes' = '#e78ac3', Microglia= '#8da0cb'))
colors$Cancer <- c(BRCA = "#fb9a99", BRCA_Basal = "#e31a1c",
                   CESC = "#fdbf6f",CRC = '#ff7f00',
                   ccRCC = "#cab2d6", 
                   GBM = "#6a3d9a" ,
                   MM = "#b15928" ,
                   HNSCC = "#b2df8a",
                   OV = "#33a02c", 
                   PDAC = "#a6cee3", 
                   UCEC = "#1f78b4" )
saveRDS(colors, '/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/Colors_panatac_v1.0.rds')
walk2(objects, out.metas, function(input.obj, output.meta.name) {
  tryCatch( {
    print(input.obj)
    input.obj <- 'PDAC_snATAC_Merged.PancanSet.noDoublets.20220209.rds'
    output.meta.name <- 'PDAC_12_samples_metadata_v6_data_freeze_v5.0.tsv'
    
    output.obj <- paste0(str_split_fixed(input.obj, '2022', 2)[,1], '20220215.rds')
    cancer.type <- str_split_fixed(input.obj, '_', 2)[1]
    
    cat('opening object\n')
    panc <- readRDS(paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/', input.obj))
    
    #sample_barcodes <- meta.to.upd[[1]]
    cells.in.object <- colnames(panc)
    # if (cancer.type=='PDAC') {
    #   previous.meta.cancer <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/PDAC_12_samples_metadata_data_freeze_v3.0.tsv') %>%
    #     select(Barcodes_cancer, cell_type.harmonized.cancer) 
    # } else if (cancer.type=='UCEC'){
    #   previous.meta.cancer <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs_correct_peaks/UCEC_12_samples_metadata_data_freeze_v3.0.tsv') %>%
    #     select(-Barcodes_cancer) %>%
    #     dplyr::rename(Barcodes_cancer=V1) %>%
    #     select(Barcodes_cancer, cell_type.harmonized.cancer) 
    # } else {
      previous.meta.cancer <- previous.meta %>%
        filter(Cancer == cancer.type) %>% 
        #data.frame(row.names = 'Barcodes_merged') %>%
        select(-Cancer) 
    #}
    head(previous.meta.cancer)
    
    meta.to.upd.cancer <- meta.to.upd %>%
      filter(Cancer == cancer.type) %>%
      select(-seurat_clusters)
    
    meta.merged <- left_join(meta.to.upd.cancer, 
                             previous.meta.cancer,
                             by = 'Barcodes_cancer')
    
    # if (cancer.type=='ccRCC') {
    #   meta.merged$Barcodes_merged <- meta.merged$Barcodes_cancer
    #   meta.merged$Barcodes_merged <- case_when(grepl('K1900070', meta.merged$Barcodes_merged) ~ str_replace(string = meta.merged$Barcodes_merged, pattern = 'ccRCC', replacement = 'PKD'),
    #                                            TRUE ~ meta.merged$Barcodes_merged)
    # } else {
    #   meta.merged$Barcodes_merged <- meta.merged$Barcodes_cancer
    # }
    
    meta.merged <- meta.merged %>% 
      data.frame(row.names = 'Barcodes_cancer')
    
    #x <- meta.merged$cell_type.harmonized.cancer=='B-cells/Plasma' & !is.na(meta.merged$cell_type.harmonized.cancer=='B-cells/Plasma')
    #meta.merged$cell_type.harmonized.cancer[x] <- meta.merged$cell_type.harmonized[x]
    
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
    
    if (cancer.type %in% c('HNSCC', 'OV', 'CESC')) {
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
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2==31] <- 'Plasma'
    }
    
    if (cancer.type=='OV') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2==29] <- 'Unknown'
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2==31] <- 'Unknown'
    }
    if (cancer.type=='PDAC') {
         panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(29)] <- 'Acinar'
    }
    
    if (cancer.type=='BRCA') {
      panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(44)] <- 'Unknown'
    }
    # 
    # if (cancer.type=='CRC') {
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(47,46,50)] <- 'Unknown'
    # }
    # 
    # if (cancer.type=='GBM') {
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(38)] <- 'Astrocytes'
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(54,52)] <- 'Unknown'
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(59)] <- 'Tumor'
    # }
    # 
    # if (cancer.type=='MM') {
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(32,33)] <- 'Unknown'
    # }
    # if (cancer.type=='UCEC') {
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(29)] <- 'Plasma'
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(10)] <- 'Macrophages'
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(21)] <- 'T-cells'
    #   
    # }
    # if (cancer.type=='PDAC') {
    #   panc$cell_type.harmonized.cancer[panc$pancan_snn_res.2 %in% c(35)] <- 'Unknown'
    # }
    
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


## 
## fix BRCA cohort
library(ggplot2)
library(tidyverse)
library(data.table)

meta <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types//v5.0_data_freeze/BRCA_19_samples_metadata_v7_data_freeze_v5.0.tsv')
meta$cell_type.harmonized.cancer[meta$pancan_snn_res.2 %in% c(44)] <- 'Unknown'
meta$cell_type.harmonized.cancer[meta$pancan_snn_res.2 %in% c(37)] <- 'Pericytes'
meta$cell_type.harmonized.cancer[meta$Piece_ID == 'HT1408-06' & meta$cell_type.harmonized.cancer == 'Tumor' & !(meta$pancan_snn_res.2 %in% c(23,14))] <- 'Unknown'
meta$cell_type.harmonized.cancer[meta$Piece_ID == 'HT1408-06' & meta$cell_type.harmonized.cancer == 'Endothelial' & meta$UMAP_1 > -5] <- 'Unknown'


fwrite(meta, '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/BRCA_19_samples_metadata_v7_data_freeze_v5.1.tsv', sep='\t', row.names = FALSE)


### fix ccRCC cohort
library(ggplot2)
library(tidyverse)
library(data.table)

meta <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/ccRCC_29_samples_metadata_v6_data_freeze_v5.0.tsv')

meta$cell_type.harmonized.cancer[meta$cell_type.harmonized.cancer %in% c('Endothelial', 'Fibroblasts') & meta$UMAP_2 <= 3.44117] <- 'Unknown'
meta[meta$Piece_ID=='C3L-00790-T1' & meta$UMAP_1 > 5  & meta$UMAP_2 <= 3.44117 & meta$cell_type.harmonized.cancer != 'Tumor', ]
fwrite(meta, '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/ccRCC_29_samples_metadata_v6_data_freeze_v5.1.tsv', sep='\t', row.names = FALSE)


meta[meta$Piece_ID == 'HT1408-06' & meta$cell_type.harmonized.cancer == 'Endothelial' & UMAP_1 > -5, ] <- 'Unknown'

###fix UCEC cohort
meta <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/UCEC_11_samples_metadata_v6_data_freeze_v5.0.tsv')

meta$cell_type.harmonized.cancer[meta$pancan_snn_res.2 == 3] <- 'Unknown'
fwrite(meta, '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/UCEC_11_samples_metadata_v6_data_freeze_v5.0.tsv', sep='\t', row.names = FALSE)


############## END #################################

#### add snRNA based CNV annotation

meta <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v5.0_data_freeze/BRCA_19_samples_metadata_v7_data_freeze_v5.0.tsv')
cnv <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/paired_sn_rna/114_snRNAsamples_CNV_count_annotation.20220214.tsv')
cnv.cancer <- cnv %>% filter(Disease == 'BRCA')
intersect(meta$V1, cnv$Barcode)

meta$CNV_count <- cnv$CNV_count[match(meta$V1, cnv$Barcode)]
meta <- meta[,-11]
sample = 'HT263B1-S1H1'


ggplot(filter(meta, Piece_ID == sample), aes (x = cell_type.harmonized.cancer, y = log10(CNV_count))) +
  geom_boxplot() +
  geom_jitter(height = 0)

meta$CNV_present <- case_when(log10(meta$CNV_count) > 2 ~ 'Yes',
                              log10(meta$CNV_count) < 2 ~ 'No',
                              TRUE ~ 'NA')


ggplot(data = filter(meta, Piece_ID == sample), aes(x = UMAP_1, y = UMAP_2, color = CNV_present)) +
  geom_point(size = 0.001)+
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  facet_wrap(~cell_type.harmonized.cancer +seurat_clusters)
ggsave(paste0('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v5.0_data_freeze/BRCA_CNV_count_', sample,'.pdf'), width = 12, height = 12, limitsize = F)




ggplot(data = filter(meta, seurat_clusters %in% c(43,37)), aes(x = UMAP_1, y = UMAP_2, color = CNV_present)) +
  geom_point(size = 0.001)+
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  facet_wrap(~cell_type.harmonized.cancer +Piece_ID)
ggsave(paste0('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v5.0_data_freeze/PDAC_cluster43_37_CNV_yes_no.pdf'), width = 6, height = 6, limitsize = F)

ggplot(data = filter(meta, seurat_clusters %in% c(43,37)), aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.001)+
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  facet_wrap(~cell_type.harmonized.cancer +Piece_ID)
ggsave(paste0('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v5.0_data_freeze/PDAC_cluster43_37_CNV_yes_no.pdf'), width = 6, height = 6, limitsize = F)

ggplot(data = filter(meta, seurat_clusters %in% c(43,37)), aes(x = UMAP_1, y = UMAP_2, color = Piece_ID)) +
  geom_point(size = 0.001)+
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  facet_wrap(~cell_type.harmonized.cancer)



