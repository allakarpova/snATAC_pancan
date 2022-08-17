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

setwd('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs')
input.obj <- 'UCEC_snATAC_Merged.PancanSet.20210822.rds'
output.obj <- 'UCEC_snATAC_Merged.PancanSet.20210825.rds'
cancer.type <- str_split_fixed(input.obj, '_', 2)[1]
input.meta.name <- 'UCEC_12_samples_metadata_data_freeze_v2.0.tsv'
output.meta.name <- 'UCEC_12_samples_metadata_data_freeze_v3.0.tsv'


panc <- readRDS(paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/RDS.50PCs/', input.obj))
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
scrub.sub$Piece_ID %>% unique
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
  verbose = T,
  group.singletons = T
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
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer.pdf"), plot = p4, height=12,width=15, useDingbats = F)

if (cancer.type=='ccRCC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(51, 57), 'Doublet', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='CESC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(21), 'Doublet', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='CRC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(45), 'Doublet', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(16), 'B-cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(11,43), 'Plasma', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(8,14,22,30), 'Normal epithelial cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(21), 'Goblet', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='GBM') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$cell_type.harmonized=='OPC', 'OPC', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='BRCA') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(14, 30, 27), 'Normal epithelial cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(24), 'B-cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(11), 'Plasma', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='MM') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(7,29), 'T-cells', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='PDAC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(31), 'Normal epithelial cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(25), 'B-cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(5,47), 'Plasma', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='UCEC') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(18,32,40), 'B-cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(19,44), 'Plasma', panc$cell_type.harmonized.cancer)
} else if (cancer.type=='OV') {
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(14,17), 'B-cells', panc$cell_type.harmonized.cancer)
  panc$cell_type.harmonized.cancer <- ifelse(panc$seurat_clusters %in% c(4,19,32), 'Plasma', panc$cell_type.harmonized.cancer)
}

panc$cell_type.harmonized.cancer <- ifelse(panc$final.doublet, 'Doublet', panc$cell_type.harmonized.cancer)
p4 <- DimPlot(panc, pt.size = ps, label=T, group.by = 'cell_type.harmonized.cancer', cols = col_vector)
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer2.pdf"), plot = p4, height=12,width=15, useDingbats = F)


# Add annotation
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
Annotation(panc) <- annotations
fwrite(panc@meta.data, output.meta.name, row.names = T, sep = '\t')
saveRDS(panc, output.obj)



### this is for breast only
brca.meta <- fread('/diskmnt/Projects/HTAN_analysis_2/BRCA/Analyses/Alla/ST_interactions/normal_meta/snRNA_combo_metadata_v1.1.tsv')
obj.meta <- panc$meta.data
obj.meta$Piece_ID <- gsub('-','_', obj.meta$Piece_ID)
obj.meta <- merge(obj.meta,brca.meta[,c('original_barcode', 'Piece_ID', 'cell_type_specific'), with = F], by.x = c('Barcodes', 'Piece_ID'), by.y = c('original_barcode', 'Piece_ID'), all.x = T)
table(obj.meta$seurat_clusters, obj.meta$cell_type_specific)[,c('B', 'Plasma')]
panc <- readRDS('BRCA_snATAC_Merged.PancanSet.20210825.rds')

## this is for PDAC
pdac.meta <- fread(output.meta.name) %>% data.frame(row.names = 1, check.rows = F, check.names = F)
panc@meta.data <- panc@meta.data[,c("orig.ident",'nCount_pancan','nFeature_pancan', 'seurat_clusters')]
panc <- AddMetaData(object = panc, metadata = pdac.meta)
Annotation(panc) <- annotations




# save UMAP embeddings
library(purrr)
library(tidyverse)
rds.files <- list.files(path = '/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/RDS.50PCs/', pattern = 'rds', full.names = T)
meta.files <- list.files(path = '.', pattern = 'tsv', full.names = T)
meta.files <- meta.files[!grepl('137', meta.files)]
#p <- './UCEC_12_samples_metadata_data_freeze_v3.0.tsv'
walk2(rds.files, meta.files, function(r, p) {
  panc <- readRDS(r)
  panc@meta.data <- panc@meta.data[,c("orig.ident",'nCount_pancan','nFeature_pancan', 'seurat_clusters')]
  s <- str_replace(string = p,pattern = '[.]/', replacement = '')
  s <- str_split_fixed(s,pattern = '_', 2)[,1]
  meta <- fread(p) %>%
    data.frame(row.names = 1, check.rows = F, check.names = F)
  panc <- AddMetaData(object = panc, metadata = meta)
  fwrite(cbind(Embeddings(panc, reduction = "umap"), meta), 
         glue::glue("{s}_UMAP_embeddings_metadata_data_freeze_v3.0.tsv"), row.names = T, sep = '\t')
  
  peak.counts <- colSums(x = GetAssayData(panc, slot= 'data'))
  panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan_data')
  VlnPlot(panc, features =  c('nCount_pancan', 'nFeature_pancan', 'peak_RF_pancan_data'), group.by = 'Piece_ID', ncol = 1)
  ggsave(glue::glue("{s}_vlnplot_nCount.pdf"), width = 7, height = 15)
})



#map new scrublet calls
setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/data.freeze.v3.0/merged_by_cancer_50PCs/update_cell_type/')
scrub <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_allATAC__15p_rate_adjusted.cutoff.tsv')
cancer.type <- 'HNSCC'
meta <- fread(paste0(cancer.type, '_UMAP_embeddings_metadata_data_freeze_v3.0.tsv'))
meta <- meta %>%
  select(-c("doublet_score.rna","predicted_doublet.rna",
            "Cutoff.rna","predicted.doublet.upd.rna",
            "doublet_score.atac",
            "predicted_doublet.atac","Cutoff.atac",
            "predicted.doublet.upd.atac","predicted.doublet.upd.rna.atac", "final.doublet"))
sample_barcodes <- meta[[1]]

meta$Piece_ID %>% unique
meta <- meta %>%
  data.frame(row.names = 'Barcodes_cancer', check.rows = F, check.names = F)

scrub.sub <- scrub  %>%
  mutate(Cancer = str_split_fixed(Sample,'_', 2)[,1]) %>%
  mutate(Cancer = gsub('PKD', 'ccRCC', Cancer)) %>%
  filter(Cancer == cancer.type)  %>%
  mutate(Piece_barcode = paste(Cancer, Piece_ID, Barcodes, sep = '_'))

scrub.sub$Piece_ID %>% unique
scrub.sub <- scrub.sub[Piece_barcode %in% rownames(meta), ]
scrub.sub <- scrub.sub %>% 
  data.frame(row.names = 'Piece_barcode', check.rows = F, check.names = F)
scrub.sub <- scrub.sub[rownames(meta), ]
scrub.sub <- scrub.sub[,c("doublet_score.rna","predicted_doublet.rna",
                          "Cutoff.rna","predicted.doublet.upd.rna",
                          "doublet_score.atac",
                          "predicted_doublet.atac","Cutoff.atac",
                          "predicted.doublet.upd.atac","predicted.doublet.upd.rna.atac", "final.doublet")]

meta.scrub <- cbind(meta, scrub.sub)
colnames(meta.scrub)[colnames(meta.scrub)=='V1'] <- 'Barcodes_merged'
ggplot(data = meta.scrub, aes(x = UMAP_1, y = UMAP_2, color = final.doublet)) +
  geom_point (shape = 16, size = 1) + 
  scale_color_manual(values = c('black', 'yellow')) +
  theme_cowplot() + guides(color = guide_legend(override.aes = list(size=6)))
ggsave(paste0(cancer.type,"_snATAC_Merged_final.doublet_new.pdf"), height=12,width=13, useDingbats = F)

meta.scrub$cell_type.harmonized.doublets <- ifelse(meta.scrub$final.doublet, 'Doublet', meta.scrub$cell_type.harmonized)
tb <- table(meta.scrub$seurat_clusters, meta.scrub$cell_type.harmonized.doublets)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
meta.scrub$cell_type.harmonized.cancer <- cluster.match.celltype[as.character(meta.scrub$seurat_clusters)]
meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$final.doublet, 'Doublet', meta.scrub$cell_type.harmonized.cancer)

#n <- length(unique(meta.scrub$cell_type.harmonized))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
unique(meta.scrub$cell_type.harmonized.cancer)
col_vector <- c('Doublet' = 'grey', 
                "Tumor" = '#e7298a', 
                "Fibroblasts" =  '#d95f02',
                "Pericytes" = '#ff7f00',
                "Endothelial" = '#66a61e',
                "Normal epithelial cells"= '#1b9e77', 
                "Endothelial" = "#e6ab02",
                "Macrophages" = '#7570b3', 
                "Mast" = '#6a3d9a',
                "T-cells" = '#1f78b4',
                "NK" = '#a6cee3',
                "Tregs" = '#b3cde3',
                "Plasma" = '#fb9a99',
                "B-cells" = '#e31a1c',
                "B-cells/Plasma" = '#cb181d',
                "DC" = '#cab2d6',
                "Oligodendrocytes" = '#33a02c',
                "Neurons" = '#b2df8a',
                "OPC" = '#fccde5',
                "Erythrocytes" = '#fbb4ae', 
                "Hepatocytes" = '#a65628',
                "Epithelial"= '#1b9e77',
                "Alveolar" = '#b3cde3',
                "Goblet" = '#fed9a6',
                "Acinar" = '#e78ac3',
                "Islets" = '#7fc97f')
cancer.col <- c("BRCA"= "#c70039", "CESC"="#ff5733", "CRC"="#ff8d1a", "GBM"="#ffc300","HNSCC"="#eddd53", "MM"="#add45c", "OV"="#57c785", "PDAC"="#00baad","UCEC"="#2a7b9b", "ccRCC"="#3d3d6b")
saveRDS(col_vector, '~/R_working_dir/scripts/snATAC/Cell_type_colors_panatac_v1.0.rds')

readRDS('~/R_working_dir/scripts/snATAC/Cell_type_colors_panatac_v1.0.rds')
ggplot(data = meta.scrub, aes(x = UMAP_1, y = UMAP_2, color = cell_type.harmonized.cancer)) +
  geom_point (shape = 16, size = 1) + 
  scale_color_manual(values = col_vector[names(col_vector) %in% unique(meta.scrub$cell_type.harmonized.cancer)]) +
  theme_cowplot() + guides(color = guide_legend(override.aes = list(size=6)))
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer_new.pdf"), height=12,width=13, useDingbats = F)

if (cancer.type=='ccRCC') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(51, 57), 'Doublet', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='CESC') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(21), 'Doublet', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='CRC') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(45), 'Doublet', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(16), 'B-cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(11,43), 'Plasma', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(8,14,22,30), 'Normal epithelial cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(21), 'Goblet', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='GBM') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$cell_type.harmonized=='OPC', 'OPC', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='BRCA') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(14, 30, 27), 'Normal epithelial cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(24), 'B-cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(11), 'Plasma', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(48), 'Doublet', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='MM') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(7,29), 'T-cells', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='PDAC') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(31), 'Normal epithelial cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(25), 'B-cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(5,47), 'Plasma', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='UCEC') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(18,32,40), 'B-cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(19,44), 'Plasma', meta.scrub$cell_type.harmonized.cancer)
} else if (cancer.type=='OV') {
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(14,17), 'B-cells', meta.scrub$cell_type.harmonized.cancer)
  meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$seurat_clusters %in% c(4,19,32), 'Plasma', meta.scrub$cell_type.harmonized.cancer)
}

meta.scrub$cell_type.harmonized.cancer <- ifelse(meta.scrub$final.doublet, 'Doublet', meta.scrub$cell_type.harmonized.cancer)
ggplot(data = meta.scrub, aes(x = UMAP_1, y = UMAP_2, color = cell_type.harmonized.cancer)) +
  geom_point (shape = 16, size = 0.5) + 
  scale_color_manual(values = col_vector[names(col_vector) %in% unique(meta.scrub$cell_type.harmonized.cancer)]) +
  theme_cowplot() + guides(color = guide_legend(override.aes = list(size=6)))
ggsave(paste0(cancer.type,"_snATAC_Merged_cell_type.harmonized.cancer2_new.pdf"), height=12,width=14, useDingbats = F)

fwrite(meta.scrub, paste0(cancer.type,'_UMAP_embeddings_metadata_data_freeze_v3.1.tsv'), row.names = T, sep = '\t')

getwd()
total.meta <- lapply (list('BRCA', 'ccRCC', 'CESC','HNSCC','CRC','GBM','MM', 'OV', 'PDAC', 'UCEC'), function(cancer.type) {
  x <- fread(paste0(cancer.type,'_UMAP_embeddings_metadata_data_freeze_v3.1.tsv'))
  colnames(x)[1] <- 'Barcodes_cancer'
  x <- x %>% data.frame(row.names = 'Barcodes_merged', check.rows = F, check.names = F)
  fwrite(x, paste0(cancer.type,'_UMAP_embeddings_metadata_data_freeze_v3.1.tsv'), row.names = T, sep= '\t')
})

total.meta <- lapply (list('HNSCC'), function(cancer.type) {
  x <- fread(paste0(cancer.type,'_UMAP_embeddings_metadata_data_freeze_v3.1.tsv'))
  colnames(x)[1] <- 'Barcodes_cancer'
  x <- x %>% data.frame(row.names = 'Barcodes_merged', check.rows = F, check.names = F)
  fwrite(x, paste0(cancer.type,'_UMAP_embeddings_metadata_data_freeze_v3.1.tsv'), row.names = T, sep= '\t')
})

total.meta <- lapply (list('BRCA', 'ccRCC', 'CESC','HNSCC','CRC','GBM','MM', 'OV', 'PDAC', 'UCEC'), function(cancer.type) {
  x <- fread(paste0(cancer.type,'_UMAP_embeddings_metadata_data_freeze_v3.1.tsv'))
  })
total.meta <- rbindlist(total.meta)


setdiff(colnames(total.meta[[2]]), colnames(total.meta[[1]]))

