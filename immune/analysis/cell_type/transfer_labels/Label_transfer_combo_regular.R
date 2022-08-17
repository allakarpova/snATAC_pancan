suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))

add_filename <- 'snATAC_scATAC_label_transferred_from_comboATAC'
output <- '/diskmnt/Projects/snATAC_analysis/immune/cell_typing/label_transfer'
setwd(output)
combo.obj <- readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v3.0/combo_only/Reclustered_immune_snATAC_Merged_44_combo_samples_v3.0_cluster_peaks.rds')
combo.meta <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_immune_44_samples.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)
combo.obj <- AddMetaData(combo.obj, combo.meta[,c('Cancer', 'Cell_type_state', 'Cell_type_markers')])
combo.obj$Cell_type_state[grepl(pattern = 'Doublet',combo.obj$Cell_type_state)] <- 'Doublet'
combo.obj$Cell_type_markers[grepl(pattern = 'Doublet',combo.obj$Cell_type_markers)] <- 'Doublet'

combo.obj$Cell_type_state[is.na(combo.obj$Cell_type_state)] <- 'Doublet'
combo.obj$Cell_type_markers[is.na(combo.obj$Cell_type_markers)] <- 'Doublet'

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/Colors_panatac_v1.0.rds')
int.obj <- readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v3.0_integrated/137_samples_v3.0_integrated_data_type.rds')

int.obj <- subset(int.obj, data.type != '10x_SC_Multi_ATAC_SEQ')

# compute UMAP and store the UMAP model
combo.obj <- RunUMAP(combo.obj, reduction = "lsi", dims = 2:50, return.model = TRUE)

plan("multicore", workers = 20)
options(future.globals.maxSize = 500 * 1024^3) # for 500 Gb RAM

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = combo.obj,
  query = int.obj, 
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:50
)

# map query onto the reference dataset
int.obj <- MapQuery(
  anchorset = transfer.anchors,
  reference = combo.obj,
  query = int.obj,
  refdata = combo.obj$Cell_type_state,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)


p1 <- DimPlot(combo.obj, reduction = "umap", group.by = "Cell_type_state", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p2 <- DimPlot(int.obj, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query")

p1 | p2
ggsave(paste0(add_filename, '_cell_type_predicted.pdf'), width = 12, height = 4.5)

p1 <- DimPlot(combo.obj, reduction = "umap", group.by = "Cancer",cols = colors$Cancer, label = TRUE, repel = TRUE, pt.size = 0.0005)  + ggtitle("Reference")
p2 <- DimPlot(int.obj, reduction = "ref.umap", group.by = "Cancer", cols = colors$Cancer, label = TRUE, repel = TRUE)  + ggtitle("Query")
p1 | p2
ggsave(paste0(add_filename, '_cancer.pdf'), width = 12, height = 4.5)



saveRDS(int.obj, paste0(add_filename, '_10142021.rds'))



