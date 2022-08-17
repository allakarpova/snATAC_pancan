# manually annotate PDAC objects
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)

rds.files <- list.files(path = '~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typing', pattern = '.rds', full.names = T, recursive = T)
f = rds.files
panc <- readRDS(f)
#sample <- str_split_fixed(f, pattern = '/PDAC_', 3)[,2]
sample <- str_split_fixed(f, pattern = 'cell_typing/|[.]', 3)[,2]
sample
DefaultAssay(panc) <- 'RNA'
dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typing')
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typing')

DimPlot(panc, label = T, group.by = 'seurat_clusters')
DimPlot(panc, label = T, group.by = 'seurat_clusters')

panc@meta.data
# plot dotplot markers MINE
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v04222021.txt', data.table = F, header = T)


genes2plot <- myeloid.genes$Gene %>% unique()

p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA')

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
p <- p + scale_color_viridis_c("viridis", direction = 1)
p
ggsave(paste0( "Dotplot_marker_gene_expression_", sample, "_RNA.pdf"),height=40,width=40,useDingbats=FALSE,limitsize = FALSE)


FeaturePlot(panc, features = subset(myeloid.genes, Gene_set =='NK cells')$Gene, order = T)
ggsave(paste0( "FeaturePlot_NK_marker_gene_expression_", sample, "_RNA.pdf"),height=20,width=20,useDingbats=FALSE,limitsize = FALSE)

FeaturePlot(panc, features = subset(myeloid.genes, Gene_set =='CD8+ cytotoxic T cells')$Gene, order = T)
ggsave(paste0( "FeaturePlot_CD8+ cytotoxic T cells_marker_gene_expression_", sample, "_RNA.pdf"),height=20,width=20,useDingbats=FALSE,limitsize = FALSE)

FeaturePlot(panc, features = subset(myeloid.genes, Gene_set =='CD4 cells')$Gene, order = T)
ggsave(paste0( "FeaturePlot_CD4 cells_marker_gene_expression_", sample, "_RNA.pdf"),height=10,width=10,useDingbats=FALSE,limitsize = FALSE)


metadata <- fread ('Everything_merge_v5_annotated_metadata.tsv', data.table = F)

toplot <- metadata %>% group_by(orig.ident, cell_type) %>% tally ()
ggplot(data=toplot, aes(x=cell_type, y= n) ) +
  geom_bar(stat = 'identity') +
  facet_wrap(~orig.ident) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


