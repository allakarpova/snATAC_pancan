suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(require(magrittr))
suppressMessages(require(readr))
suppressMessages(library(Matrix))
suppressMessages(library(tidyr))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(data.table))
library(stringr)

setwd('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types')
input.obj <- 'GBM_snATAC_Merged.PancanSet.20210818.rds'
output.obj <- 'GBM_snATAC_Merged.PancanSet.20210824.rds'
cancer.type <- str_split_fixed(input.obj, '_', 2)[1]
input.meta.name <- 'GBM_20_samples_metadata_data_freeze_v2.0.tsv'
output.meta.name <- 'GBM_20_samples_metadata_data_freeze_v2.0.tsv'


panc <- readRDS(paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/', input.obj))
scrub <- fread('/diskmnt/Projects/snATAC_analysis/scrublet/Scrublet_allATACadjusted.cutoff.tsv')
meta <- fread(paste0('/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v4.0/',input.meta.name))
sample_barcodes <- meta[[1]]

panc$Piece_ID%>% unique
meta$Piece_ID %>% unique
meta <- meta %>%
  data.frame(row.names = 1, check.rows = F, check.names = F) %>%
  dplyr::rename(seurat_clusters_indiv = seurat_clusters) %>%
  mutate(Piece_barcode = paste(Cancer, Piece_ID, V1, sep = '_')) %>%
  data.frame(row.names = 'Piece_barcode', check.rows = F, check.names = F)

scrub.sub <- scrub  %>%
  mutate(Cancer = str_split_fixed(Sample,'_', 2)[,1]) %>%
  mutate(Cancer = gsub('PKD', 'ccRCC', Cancer)) %>%
  filter(Cancer == cancer.type)  %>%
  mutate(Piece_barcode = paste(Cancer, Piece_ID, Barcodes, sep = '_'))
scrub.sub <- scrub.sub[Piece_barcode %in% rownames(meta), ]
scrub.sub <- scrub.sub %>% data.frame(row.names = 'Piece_barcode', check.rows = F, check.names = F)
scrub.sub <- scrub.sub[rownames(meta), ]
scrub.sub <- scrub.sub[,c("doublet_score.rna","predicted_doublet.rna",
        "Cutoff.rna","predicted.doublet.upd.rna",
        "doublet_score.atac",
        "predicted_doublet.atac","Cutoff.atac",
        "predicted.doublet.upd.atac","predicted.doublet.upd.rna.atac", "final.doublet")]

dim(scrub.sub)
dim(meta)
dim(panc)

setequal(colnames(panc), rownames(meta))

setequal(rownames(scrub.sub), rownames(meta))
meta.scrub <- cbind(meta, scrub.sub)

colnames(meta.scrub)[colnames(meta.scrub)=='V1'] <- 'Barcodes'
fwrite(meta.scrub, output.meta.name, row.names = T, sep = '\t')

if (cancer.type=='BRCA') {
  obj.meta <- panc@meta.data %>% 
    mutate(Barcodes = rownames(.)) %>% 
    mutate (Barcodes = gsub('_', '-', Barcodes)) %>%
    mutate (Barcodes = str_split_fixed(Barcodes, '-', 4)[,4]) %>%
    mutate(Sample_barcodes = paste(dataset, Barcodes, sep='_'))
  setequal(obj.meta$Sample_barcodes, sample_barcodes)
} else if (cancer.type %in% c('ccRCC', 'CESC', 'CRC', 'GBM')) {
  obj.meta <- panc@meta.data %>% 
    mutate(Barcodes = rownames(.)) %>% 
    mutate (Barcodes = str_split_fixed(Barcodes, '_', 3)[,3]) %>%
    mutate(Sample_barcodes = paste(dataset, Barcodes, sep='_'))
  setequal(obj.meta$Sample_barcodes, sample_barcodes)
}


if (cancer.type=='BRCA' | cancer.type=='ccRCC') {
  panc <- RenameCells(panc, new.names = rownames(meta.scrub))
}

setequal(colnames(panc), rownames(meta.scrub))
panc@meta.data <- panc@meta.data[,c("orig.ident",'nCount_pancan','nFeature_pancan', 'seurat_clusters')]
panc <- AddMetaData(object = panc, metadata = meta.scrub)

n <- length(unique(panc$cell_type.harmonized))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ps <- ifelse(ncol(panc)>50000, 0.1, 1)

p5 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'final.doublet', cols = c('black', 'yellow'))
ggsave(paste0(cancer.type,"_snATAC_Merged_final.doublet.pdf"), plot = p5, height=12,width=13, useDingbats = F)
p3 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized', cols = col_vector)
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.pdf"), plot = p3, height=12,width=15, useDingbats = F)
p1 <- DimPlot(panc, group.by = 'Piece_ID', pt.size =ps,label = T)
ggsave(paste0(cancer.type,"_snATAC_Merged_Piece_ID.pdf"),plot = p1,height=12,width=13, useDingbats = F)

p1 <- DimPlot(panc, split.by = 'Piece_ID', pt.size =ps,label = T, ncol = 5)
ggsave(paste0(cancer.type,"_snATAC_Merged_Piece_ID_split.pdf"),plot = p1,height=40,width=40, useDingbats = F)

p2 <- DimPlot(panc, group.by = 'seurat_clusters', pt.size = ps,label=T)
ggsave(paste0(cancer.type,"_snATAC_Merged_clusters0.8.pdf"), plot = p2, height=12,width=14, useDingbats = F)

total_fragments_cell <- panc$passed_filters
peak.counts <- colSums(x = GetAssayData(panc, slot='counts'))
frip <- peak.counts *100 / total_fragments_cell
panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_pancan')
panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan')



panc <- FindClusters( 
  object = panc,
  algorithm = 3,
  resolution = 2,
  verbose = FALSE
)
p2 <- DimPlot(panc, group.by = 'seurat_clusters', pt.size = ps,label=T)
ggsave(paste0(cancer.type,"_snATAC_Merged_clusters2.pdf"), plot = p2, height=12,width=14, useDingbats = F)

#panc@commands
panc$cell_type.harmonized.doublets <- ifelse(panc$final.doublet, 'Doublet', panc$cell_type.harmonized)

p3 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.doublets', cols = col_vector)
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.doublets.pdf"), plot = p3, height=12,width=15, useDingbats = F)


tb <- table(panc$seurat_clusters, panc$cell_type.harmonized.doublets)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
panc$cell_type.harmonized.cancer <- cluster.match.celltype[as.character(panc$seurat_clusters)]
panc$cell_type.harmonized.cancer <- ifelse(panc$final.doublet, 'Doublet', panc$cell_type.harmonized.cancer)

p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = col_vector)
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer2.pdf"), plot = p4, height=12,width=15, useDingbats = F)

if (cancer.type=='ccRCC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(48,49), 'Doublet', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='CESC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(17,23), 'Doublet', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='CRC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(42,46,45,47,48,52, 54), 'Doublet', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(13,29), 'Epithelial', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(21), 'Goblet', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='GBM') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(43), 'OPC', panc$cell_type.harmonized.cancer)
}
panc@meta.data %>% select(-cell_type.harmonized.doublets)




p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = col_vector)
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer2.pdf"), plot = p4, height=12,width=15, useDingbats = F)


fwrite((panc@meta.data %>% select(-cell_type.harmonized.doublets)), output.meta.name, row.names = T, sep = '\t')
saveRDS(panc, output.obj)


