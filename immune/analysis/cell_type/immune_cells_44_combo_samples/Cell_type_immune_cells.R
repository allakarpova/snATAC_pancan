# cell type immune cells based on RNA from 44 combo samples
library(Seurat)
library(data.table)
library(dplyr)
library(scatterpie)
library(paletteer)
library(ggplot2)
library(readxl)

panc <- readRDS ('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA/44_snRNA_combo_Merged_immune_RNA_immune_cells.rds')

add_filename <- 'comboRNA_immune_cells'
wd <- paste0('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/', add_filename)
dir.create(wd, recursive = T)
setwd(wd)

DimPlot(panc,label = T)
ggsave(paste0( "Dimplot_clusters.pdf"), height=5,width=6,useDingbats=FALSE,limitsize = FALSE)

# plot dotplot markers MINE
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v06222021.txt', data.table = F, header = T)
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
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)

ggsave(paste0( "Dotplot_marker_gene_expression_", add_filename, "_RNA.pdf"),plot=p,height=70,width=50,useDingbats=FALSE,limitsize = FALSE)


panc$cell_type_upd <- case_when(panc$seurat_clusters == 3 ~ 'B-cells', 
                            
                            panc$seurat_clusters %in% c(32,37) ~ 'IgG- IgA+ plasma cells',
                            panc$seurat_clusters %in% c(17,18,20,23,28) ~ 'IgG+ plasma cells',
                            panc$seurat_clusters %in% c(7) ~ 'Plasma cells',
                            panc$seurat_clusters == 12 ~ 'NK cells',
                            panc$seurat_clusters == 5 ~ 'Tregs',
                            panc$seurat_clusters %in% c(2,13,14)~ 'CD4 T-cells',
                            panc$seurat_clusters == 1 ~ 'CD8 T-cells',
                            panc$seurat_clusters == 34 ~ 'T-cells likely',
                            panc$seurat_clusters %in% c(0,4,8,9,10,11) ~ 'Macrophages',
                            panc$seurat_clusters == 26 ~ 'pDC',
                            panc$seurat_clusters == 16 ~ 'Doublet/Fibroblasts',
                            panc$seurat_clusters %in% c(6,15,19,25,27,30,35,36, 31,24, 22) ~ 'Doublet/Epithelial',
                            panc$seurat_clusters %in% c(29) ~ 'Doublet/Hepatocytes',
                            panc$seurat_clusters %in% c(33) ~ 'Doublet/ADM',
                            panc$seurat_clusters == 21 ~ 'Mast',
                            TRUE ~ 'Unknown')

panc@meta.data %>% subset(cell_type_upd == 'IgG+ plasma cells' & Cancer =='OV') 
panc@meta.data$Piece_ID %>% unique()
library(RColorBrewer)
n <- length(unique(panc$cell_type))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


DimPlot(panc, group.by = 'cell_type_upd', label = T, cols = col_vector) 
ggsave(paste0( "Dimplot_cell_type_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

DimPlot(panc, group.by = 'cell_type_upd', label = T, cols = col_vector, cells = rownames(subset(panc@meta.data, !grepl('Doublet', cell_type_upd)))) 
ggsave(paste0( "Dimplot_cell_type_", add_filename, "_no_doublet.pdf"),height=7,width=10,useDingbats=FALSE)

length(rownames(subset(panc@meta.data, !grepl('Doublet', cell_type_upd))))
####### human liver single cell paper
degs <- read_excel('~/lab_Ding/work/single_cell/senescence/snATAC_combo/papers/41467_2018_6318_MOESM5_ESM_formatted.xlsx', col_names = T)

genes2plot <- degs$Gene %>% unique()
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA')
p$data <- merge(p$data, degs[c(1,3)], by.x = 'features.plot', by.y = 'Gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell_type , scales = "free",  drop = T, ncol = 3)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_viridis_c("viridis", direction = 1)

ggsave(paste0( "Dotplot_20human_liver_paper_markers_gene_expression_", add_filename, "_many_RNA.pdf"),plot = p, height=45,width=20,useDingbats=FALSE,limitsize = FALSE)


## add doublets
colnames(panc)
doublets <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_comboATAC_comboRNA_scores_adjusted.cutoff.tsv') %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)

db <- doublets
db <- db[colnames(panc),]
panc <- AddMetaData(object = panc, metadata = db$predicted.doublet.upd.rna.atac, col.name = 'predicted.doublet.upd.rna.atac')
panc <- AddMetaData(object = panc, metadata = db$predicted.doublet.upd.rna, col.name = 'predicted.doublet.upd.rna')
panc <- AddMetaData(object = panc, metadata = db$predicted.doublet.upd.atac, col.name = 'predicted.doublet.upd.atac')

panc$predicted.doublet.upd.rna.atac <- factor(panc$predicted.doublet.upd.rna.atac, levels = c(FALSE,TRUE))
DimPlot(panc, group.by = 'predicted.doublet.upd.rna.atac', label = F, cols = c('yellow', 'black'))
ggsave(paste0( "Dimplot_predicted.doublet.upd.rna.atac_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

DimPlot(panc, group.by = 'predicted.doublet.upd.rna', label = F, cols = c('black', 'yellow')) 
ggsave(paste0( "Dimplot_predicted.doublet.upd.rna_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

DimPlot(panc, group.by = 'predicted.doublet.upd.atac', label = F, cols = c('black', 'yellow')) 
ggsave(paste0( "Dimplot_predicted.doublet.upd.atac_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

table(panc$seurat_clusters,  panc$predicted.doublet.upd)
table(panc$seurat_clusters,  panc$predicted_doublet)
getwd()

ggplot(data = panc@meta.data[,c('seurat_clusters','predicted.doublet.upd.rna.atac')], aes(x=as.character(seurat_clusters), fill =  predicted.doublet.upd.rna.atac)) +
  geom_bar(stat = 'count', position = 'fill') +
  scale_fill_manual(values = c( 'black','yellow'))
ggsave(paste0( "Barplot_predicted.doublet.upd.rna.atac_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

ggplot(data = panc@meta.data[,c('seurat_clusters','predicted.doublet.upd.rna')], aes(x=as.character(seurat_clusters), fill =  predicted.doublet.upd.rna)) +
  geom_bar(stat = 'count', position = 'fill') +
  scale_fill_manual(values = c('black', 'yellow'))
ggsave(paste0( "Barplot_predicted.doublet.upd.rna_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

ggplot(data = panc@meta.data[,c('seurat_clusters','predicted.doublet.upd.atac')], aes(x=as.character(seurat_clusters), fill =  predicted.doublet.upd.atac)) +
  geom_bar(stat = 'count', position = 'fill') +
  scale_fill_manual(values = c('black', 'yellow'))
ggsave(paste0( "Barplot_predicted.doublet.upd.atac_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

panc@meta.data %>% select (seurat_clusters, predicted.doublet.upd.rna.atac) %>%
  group_by (seurat_clusters, predicted.doublet.upd.rna.atac) %>% tally
doublet.cluster <- rownames(table(panc$seurat_clusters, panc$predicted.doublet.upd.rna.atac) %>% as.data.frame.matrix() %>% 
  mutate(pct = `TRUE`/(`TRUE`+`FALSE`)) %>% filter(pct > 0.2)) %>% as.numeric

panc$doublet.cluster <- case_when(panc$seurat_clusters %in% doublet.cluster ~ '20%doublets', TRUE ~ 'few_doublets')
DimPlot(panc, group.by = 'doublet.cluster', label = F, cols = c('yellow', 'black')) 
ggsave(paste0( "Dimplot_doublet.cluster0.2_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)


panc$cell_type_upd <- case_when(panc$doublet.cluster == '20%doublets' ~ 'Doublet', 
                                panc$predicted.doublet.upd.rna.atac == 'TRUE' ~ 'Doublet',
                                TRUE ~ panc$cell_type_upd)
DimPlot(panc, group.by = 'cell_type_upd', label = T) 
fwrite(panc@meta.data, '44_snRNA_combo_Merged_immune_RNA_immune_cells_metadata.tsv', row.names = T, sep = '\t')

table(panc$cell_type_upd)




