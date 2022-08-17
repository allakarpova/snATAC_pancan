library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(tidyverse)
library(RColorBrewer)

###Load in seurat object
obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA_no_cell_cycle/snRNA_combo_Merged_is_Myeloid_lin_44_samples.rds")
obj <- obj %>% NormalizeData()
add_filename <- 'is_Myeloid_lin' 
dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/')
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/')

dim(obj)

# plot dotplot markers MINE
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v09012021.txt', data.table = F, header = T)
#myeloid.genes <- myeloid.genes %>% filter(grepl('Activ|T-cells|CD4|CD8|Treg|Cost|Coin|Cytoki|emory|Trm|xhau|NK|naive|Th|Tfh|Term|Effect|Prog', Gene_set))
genes2plot <- myeloid.genes$Gene %>% unique() %>% sort

p <- DotPlot(object = obj, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA')
p$data <- merge(p$data, myeloid.genes[1:2], by.x = 'features.plot', by.y = 'Gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Gene_set , scales = "free",  drop = T, ncol = 7)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)

ggsave(paste0( "Dotplot_marker_gene_expression_", add_filename, "_RNA.pdf"),plot=p,height=70,width=50,useDingbats=FALSE,limitsize = FALSE)

FindMarkers(obj, ident.1 = 27, test.use = 'LR')

walk(c('PRDM16', 'ANKFN1', 'ST6GALNAC5', 'IDO2', 'FLT3', 'BACH2', 'TANC1', 'ITGAE', 'ITGAX'), function(gene) {
  FeaturePlot(obj, features = gene, order = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0('Feature_plot_cluster_27_markers_', gene,'_', add_filename, '.pdf'),height=5,width=6,useDingbats=FALSE )
  
})

walk(c('CD14', 'FCGR3A', 'CD68'), function(gene) {
  FeaturePlot(obj, features = gene, order = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0('Feature_plot_', gene,'_', add_filename, '.pdf'),height=5,width=6,useDingbats=FALSE )
  
})

FeaturePlot(obj, features = 'CD4', order = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

cell.typed <- fread('Cell type_M.txt')
obj$Cell_type_state <- cell.typed$Cell_type_state[match(obj$seurat_clusters, cell.typed$Cluster)]
obj$Cell_type_markers <- cell.typed$Cell_type_markers[match(obj$seurat_clusters, cell.typed$Cluster)]

DimPlot(obj, group.by = 'Cell_type_state', label = T, cols = 'Spectral')
ggsave(paste0( "Dimplot_Cell_type_state_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)


DimPlot(obj, group.by = 'Cell_type_markers', label = T, cols = colorRampPalette(brewer.pal(n = 11, 'Spectral'))(18))
ggsave(paste0( "Dimplot_Cell_type_markers_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

obj$Barcodes <- map2_chr(rownames(obj@meta.data), paste0(obj$Sample, '_'), ~ str_replace(.x,.y,''))
obj$Barcodes_cancer <- paste(obj$Cancer, obj$Piece_ID, obj$Barcodes, sep='_')
fwrite(obj@meta.data, paste0('snRNA_combo_Merged_', add_filename,'_44_samples.metadata.tsv'), sep = '\t', row.names = T)

