library(future)

plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(doParallel))


setwd('/diskmnt/Projects/snATAC_analysis/immune/obj/v4.0_script/mark_doublets')
combined <- readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v4.0_script/immune_cells_integrated_cluster_peaks_chemistry.rds')
combined@reductions

metas.combo <- c('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/v4.0_comboRNA/B-cell/Celltypedata_v3.0_doublets_snRNA_combo_Merged_is_Bcell_lin.tsv',
                 '/diskmnt/Projects/snATAC_analysis/immune/cell_typing/v4.0_comboRNA/T-cell/Celltypedata_v3.0_doublets_snRNA_combo_Merged_is_Tcell_lin.tsv',
                 '/diskmnt/Projects/snATAC_analysis/immune/cell_typing/v4.0_comboRNA/Myeloid/Celltypedata_v3.0_doublets_snRNA_combo_Merged_is_Myeloid_lin.tsv')

metas <- map (metas.combo, ~fread(.x, header = T))
metas <- rbindlist(metas) %>% data.frame(row.names = 1)

combined <- AddMetaData(combined, metas)
DimPlot(combined, group.by = 'Cell_type_state', label = T) + DimPlot(combined, label = T)
ggsave('Dimplot_Cell_type_state_from_combo_immune_cells_integrated.pdf', width = 14, height = 6)

DimPlot(combined, group.by = 'data.type') + DimPlot(combined, label = T)
ggsave('Dimplot_data.type_from_combo_immune_cells_integrated.pdf', width = 14, height = 6)

combined$Cell_type_combo_reg_doublets <- case_when(combined$data.type == '10x_SC_Multi_ATAC_SEQ' & is.na(combined$Cell_type_state) ~ 'Doublet',
                                      TRUE ~ combined$Cell_type_state)

tb <- table(combined$seurat_clusters, combined$Cell_type_state)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
combined$Cell_type_combo_reg_doublets <- cluster.match.celltype[as.character(combined$seurat_clusters)]
combined$Cell_type_combo_reg_doublets[combined$seurat_clusters %in% c(34,18,43)] <- 'Microglia'
combined$Cell_type_combo_reg_doublets[combined$seurat_clusters %in% c(14)] <- 'CD8 T-cells'
combined$Cell_type_combo_reg_doublets[combined$Cell_type_state=='Doublet'] <- 'Doublet'

DimPlot(combined, group.by = 'Cell_type_combo_reg_doublets', label = T) + DimPlot(combined, label = T)
ggsave('Dimplot_Cell_type_combo_reg_doublets_from_combo_immune_cells_integrated.pdf', width = 14, height = 6)


combo.cells <- rownames(filter(combined@meta.data, data.type == '10x_SC_Multi_ATAC_SEQ'))

genes.toplot <- c('PDCD1', 'TMEM119', 'CTLA4', 'HAVCR2', 'CD163')
genes.toplot %>% walk (function(g) {
  exons <- Annotation(combined) %>% data.frame() %>% filter(gene_name == g & type =='exon') 
  chr <- exons$seqnames[1]
  st <- exons$start[1]
  en <- exons$end[nrow(exons)]
  region <- paste(chr, st-1000, en+1000, sep ='-')

  CoveragePlot(combined, region =  region,peaks = T, group.by = 'Cell_type_combo_reg_doublets')
  ggsave(glue::glue('CoveragePlot_{g}_Cell_type_combo_reg_doublets_from_combo_immune_cells_integrated.pdf'), width = 14, height = 14)
}) 


exons <- Annotation(combined) %>% data.frame() %>% filter(gene_name == g & type =='exon') 
chr <- exons$seqnames[1]
st <- exons$start[1]
en <- exons$end[nrow(exons)]
region <- paste(chr, st-1000, en+1000)

CoveragePlot(combined, region =  region,peaks = T, group.by = 'Cell_type_combo_reg_doublets')
ggsave(glue::glue('CoveragePlot_{g}_Cell_type_combo_reg_doublets_from_combo_immune_cells_integrated.pdf'), width = 14, height = 14)


fwrite((combined@meta.data %>% dplyr::select(Cell_type_state, Cell_type_combo_reg_doublets)), 
       'Celltypedata_v3.0_doublets_immune_cells_integrated_cluster_peaks_chemistry.tsv',sep='\t', row.names = T)




