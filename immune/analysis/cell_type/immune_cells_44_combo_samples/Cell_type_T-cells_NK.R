library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(tidyverse)
library(RColorBrewer)

###Load in seurat object
obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA/snRNA_combo_Merged_is_Tcell_lin_44_samples.rds")
DefaultAssay(obj) <- 'RNA'
add_filename <- 'is_Tcell_lin' 
dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells/singleR')
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells/singleR')

dim(obj)

# plot dotplot markers MINE
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v09012021.txt', data.table = F, header = T)
myeloid.genes <- myeloid.genes %>% filter(grepl('Activ|T-cells|CD4|CD8|Treg|Cost|Coin|Cytoki|emory|Trm|xhau|NK|naive|Th|Tfh|Term|Effect|Prog', Gene_set))
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

getwd()

DefaultAssay(obj) <- 'SCT'
obj <- obj %>% FindNeighbors(dims = 1:50)
obj <- obj %>% FindClusters(verbose = T, resolution = 2)
obj <- obj %>% RunUMAP(dims = 1:50, assay = 'SCT') 

Graphs(obj)
obj@commands
DimPlot(obj, label = T)
ggsave(paste0( "Dimplot_clusters_res2_pc50_", add_filename, ".pdf"),height=7,width=8,useDingbats=FALSE)

DimPlot(obj, group.by = 'BlueprintEncodeData', label = F, cols = col_vecor)
ggsave(paste0( "Dimplot_BlueprintEncodeData_res2_pc50_", add_filename, ".pdf"),height=7,width=11,useDingbats=FALSE)

VlnPlot(obj, features = c('nCount_RNA', 'percent.mt'))

p <- DotPlot(object = obj, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'SCT')
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

ggsave(paste0( "Dotplot_res2_pc50_marker_gene_expression_", add_filename, "_SCT.pdf"),plot=p,height=25,width=50,useDingbats=FALSE,limitsize = FALSE)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
FeaturePlot(obj, features = c('TOX','TBX21', 'PDCD1'), order = T)

avg.obj <- AverageExpression(object = obj, assays = 'RNA',return.seurat = T)
avg.obj <- ScaleData(avg.obj)

toplot <- FetchData(avg.obj,slot = 'scale.data',vars = genes2plot) %>% t
p <- pheatmap::pheatmap (toplot,
                         cellwidth = 6, 
                         cellheight =6.5,
                         cutree_cols = 12,
                         cutree_rows = 6,
                         #color=color.palette,
                         #gaps_row = c(3,12, 16),
                         #gaps_row = c(38),
                         #gaps_col = c(53,60,85, 100, 125, 131),
                         na_col = 'white',
                         cluster_rows = T,
                         cluster_cols = T, 
                         clustering_distance_rows = "euclidean", 
                         clustering_distance_cols = "euclidean", 
                         clustering_method = "ward.D2", 
                         #annotation_row = a[3:4],
                         #annotation_colors = annotation_colors, 
                         border_color = NA, 
                         #breaks=breaks, 
                         fontsize_col  = 7,
                         fontsize_row  = 7,
                         show_rownames= T,
                         show_colnames = T,
                         #display_numbers = heatmap.mut.table,
                         fontsize_number = 5)
pdf ('Heatmap_average_expr_T-cell_genes.pdf', width = 8, height = 20)
print(p)
dev.off()

cell.typed <- fread('Cell_types_with_cycling_cells_50PC_res2.txt')
obj$cell_type_upd <- cell.typed$`Cell-type`[match(obj$seurat_clusters, cell.typed$Cluster)]
DimPlot(obj, group.by = 'cell_type_upd', label = T, cols = col_vecor)
ggsave(paste0( "Dimplot_cell_type_upd_res2_pc50_", add_filename, ".pdf"),height=7,width=12,useDingbats=FALSE)

save.it <- data.frame(Barcodes = colnames(obj), cell_type_upd = obj$cell_type_upd)
rownames(save.it) <- save.it$Barcodes


### use object where cell cycle is regressed out
obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA_no_cell_cycle/snRNA_combo_Merged_is_Tcell_lin_44_samples.rds")
DefaultAssay(obj) <- 'RNA'
add_filename <- 'is_Tcell_lin_no_cell_cycle' 
dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/', recursive = T)
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/')

colnames(obj)
save.it <- save.it[colnames(obj),]
obj <- AddMetaData(obj, save.it$cell_type_upd, col.name = 'cell_type_upd')
DimPlot(obj, group.by = 'cell_type_upd', label = T, cols = col_vecor, label.size = 2)
ggsave(paste0( "Dimplot_cell_type_upd_res2_pc30_", add_filename, ".pdf"),height=7,width=12,useDingbats=FALSE)

DefaultAssay(obj) <- 'SCT'
obj <- obj %>% FindNeighbors(dims = 1:50)
obj <- obj %>% FindClusters(verbose = T, resolution = 2)
obj <- obj %>% RunUMAP(dims = 1:50, assay = 'SCT') 

DimPlot(obj, group.by = 'cell_type_upd', label = T, cols = col_vecor, label.size = 2)
ggsave(paste0( "Dimplot_cell_type_upd_res2_pc50_", add_filename, ".pdf"),height=7,width=12,useDingbats=FALSE)

tb <- table(obj$seurat_clusters, obj$cell_type_upd)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
obj$cell_type_upd_cluster <- cluster.match.celltype[as.character(obj$seurat_clusters)]

DimPlot(obj, group.by = 'cell_type_upd_cluster', label = T, cols = col_vecor, label.size = 2)
ggsave(paste0( "Dimplot_cell_type_upd_cluster_res2_pc50_", add_filename, ".pdf"),height=7,width=12,useDingbats=FALSE)

DimPlot(obj, label = T)
ggsave(paste0( "Dimplot_cluster_res2_pc50_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

#add my classification
cell.typed <- fread('Cell_types_T_with_cycling_cells_50PC_res2.txt')
obj$Cell_type_state <- cell.typed$Cell_type_state[match(obj$seurat_clusters, cell.typed$Cluster)]
obj$Cell_type_markers <- cell.typed$Cell_type_markers[match(obj$seurat_clusters, cell.typed$Cluster)]

DimPlot(obj, group.by = 'Cell_type_state', label = T)
ggsave(paste0( "Dimplot_Cell_type_state_res2_pc50_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

DimPlot(obj, group.by = 'Cancer', label = T, cols = col_vecor)
ggsave(paste0( "Dimplot_Cancer_res2_pc50_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

toplot <- data.frame(Cancer = obj$Cancer, Cell_type_state = obj$Cell_type_state)
ggplot(data = toplot, aes(x=Cancer, fill = Cell_type_state)) +
  geom_bar(position = 'fill', stat = 'count') +
  scale_fill_manual(values = col_vecor) +
  theme_cowplot()
ggsave('Barplot_fill_Cell_type_state_in_cancer.pdf',height=7,width=10,useDingbats=FALSE)

ggplot(data = toplot, aes(x=Cancer, fill = Cell_type_state)) +
  geom_bar(position = 'stack', stat = 'count') +
  scale_fill_manual(values = col_vecor)+
  theme_cowplot()
ggsave('Barplot_stack_Cell_type_state_in_cancer.pdf',height=7,width=10,useDingbats=FALSE)


ggplot(data = toplot, aes(fill=Cancer, x = Cell_type_state)) +
  geom_bar(position = 'fill', stat = 'count') +
  scale_fill_brewer(palette = 'Set1')+
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Barplot_fill_cancer_in_Cell_type_state.pdf',height=10,width=10,useDingbats=FALSE)

obj$Barcodes <-  map2_chr(rownames(obj@meta.data), paste0(obj$Sample, '_'), ~ str_replace(.x,.y,''))
obj$Barcodes_cancer <- paste(obj$Cancer, obj$Piece_ID, obj$Barcodes, sep='_')
fwrite(obj@meta.data, 'snRNA_combo_Merged_is_Tcell_lin_44_samples.res2.metadata.tsv', sep = '\t', row.names = T)

markers <- c('FOXP3','IL2RA','CTLA4',
             'CD8A', 'IFNG','FASLG', 'GZMA', 'GZMK','KLRG1', 'CCR5',
             'GNLY','PRF1','NKG7', 'NCAM1','FCGR3A', 'KLRB1', 'KLRF1',
             'TCF7', 'CXCR5', 'SLAMF6',
             'TOX', 'TBX21', 'EOMES', 'HAVCR2', 
             'PDCD1', 'CD38',
             'ITGAE', 'RUNX3', 'ENTPD1',
             'CD4','BTLA','CXCL13', 'CD200', 'BCL6',
             'IL7R', 'SELL', 'CCR7')
VlnPlot(obj, group.by = 'Cell_type_state', features=markers,stack = T,flip = T,assay = "RNA")

DoHeatmap(obj, group.by = 'Cell_type_state', features = markers)

obj <- RenameCells(object = obj, new.names = paste(colnames(obj), obj$Cell_type_state, sep ='_'))
Idents(obj) <- 'Cell_type_state'
avg.obj <- AverageExpression(object = obj, assays = 'RNA',slot = 'data', return.seurat = T)
avg.obj <- ScaleData(avg.obj)
DoHeatmap(avg.obj, draw.lines = F, features = markers)


toplot <- FetchData(avg.obj,slot = 'scale.data',vars = markers) %>% t
p <- pheatmap::pheatmap (toplot,
                         cellwidth = 12, 
                         cellheight =12,
                         cutree_cols = 1,
                         cutree_rows = 1,
                         #color=color.palette,
                         #gaps_row = c(3,12, 16),
                         #gaps_row = c(38),
                         #gaps_col = c(53,60,85, 100, 125, 131),
                         na_col = 'white',
                         cluster_rows = T,
                         cluster_cols = T, 
                         clustering_distance_rows = "euclidean", 
                         clustering_distance_cols = "euclidean", 
                         clustering_method = "ward.D2", 
                         #annotation_row = a[3:4],
                         #annotation_colors = annotation_colors, 
                         border_color = NA, 
                         #breaks=breaks, 
                         fontsize_col  = 11,
                         fontsize_row  = 11,
                         show_rownames= T,
                         show_colnames = T,
                         #display_numbers = heatmap.mut.table,
                         fontsize_number = 5)
p
pdf ('Heatmap_average_expr_T-cell_genes_markers.pdf', width = 10, height = 20)
print(p)
dev.off()

getwd()
