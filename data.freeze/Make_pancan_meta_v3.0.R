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

panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/137Samples_PanCan_merged_obj/137_snATAC_113K_peaks_diffPCs.motifsAdded.chromvar.20210826.rds.gz')
meta <- fread('/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/All_137_samples_metadata_data_freeze_v3.0.tsv')

meta <- data.frame(meta, row.names = 1, check.rows = F, check.names = F)

panc@meta.data <- panc@meta.data[,c('nCount_pancan_s','nFeature_pancan_s')]
panc <- AddMetaData(object = panc, metadata = meta)
total_fragments_cell <- panc$passed_filters
peak.counts <- colSums(x = GetAssayData(panc, slot='counts'))
frip <- peak.counts *100 / total_fragments_cell
panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_pancan_s')
panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_RF_pancan_s')
fwrite(panc@meta.data, '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/All_137_samples_metadata_data_freeze_v3.1.tsv', row.names = T, sep = '\t')



DimPlot(panc, group.by = 'data.type', cols = 'Paired')
ggsave('Dimplot_all_v3.0_data.type.pdf', width = 10, height = 10)
DimPlot(panc, group.by = 'data.type',split.by ='data.type', cols = 'Paired')
ggsave('Dimplot_all_v3.0_data.type_split.pdf', width = 20, height = 10)

DimPlot(panc, group.by = 'cell_type.harmonized.cancer', cols = colors)
ggsave('Dimplot_all_v3.0_cell_type.harmonized.cancer.pdf', width = 14, height = 10)
DimPlot(panc, group.by = 'Cancer', cols = colors$Cancer)
ggsave('Dimplot_all_v3.0_Cancer.pdf', width = 11, height = 10)


panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/RDS.50PCs/UCEC_snATAC_Merged.PancanSet.20210822.rds')
meta.file <- '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/UCEC_12_samples_metadata_data_freeze_v3.0.tsv'
meta <- fread(meta.file)
colnames(meta)[1] <- 'Barcodes_cancer'
meta <- cbind(Barcodes_merged = colnames(panc), meta)
meta[,1:4]
fwrite(meta, meta.file, sep = '\t')