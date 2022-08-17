suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))

panc <- readRDS('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/v1.2/Reclustered_immune_snATAC_Merged_94_sample_obj.v1.2_new_macs2.rds')

wd <- '/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/v1.2'
dir.create(wd)
setwd(wd)
add_filename <- 'snATAC_Merged_94_sample_obj.v1.2_new_macs2'
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
ggsave(paste0( "Dotplot_marker_gene_expression_", add_filename, "_RNA.pdf"),plot = p, height=70,width=50,useDingbats=FALSE,limitsize = FALSE)

DimPlot(panc, label = T) + DimPlot(panc, label = T, group.by = 'cell_type')
ggsave(paste0( "Dimplot_clusters_cell_type.pdf"), height=4,width=12,useDingbats=FALSE,limitsize = FALSE)

panc$cell_type_upd <- case_when(panc$seurat_clusters %in% c(0,2,7, 19) ~ 'Macrophages', 
                                panc$seurat_clusters %in% c(8) ~ 'Microglia', 
                                panc$seurat_clusters %in% c(4, 26, 22) ~ 'Monocytes', 
                            panc$seurat_clusters %in% c(11,25,15) ~ 'B-cells',
                            panc$seurat_clusters %in% c(17, 28) ~ 'NK cells',
                            panc$seurat_clusters %in% c(14) ~ 'CD8 CTL',
                            panc$seurat_clusters %in% c(23) ~ 'Exhausted CD8 T-cells',
                            panc$seurat_clusters %in% c(3) ~ 'CD8 T-cells',
                            panc$seurat_clusters %in% c(6) ~ 'CD4 T-cells',
                            panc$seurat_clusters %in% c(12,1) ~ 'T-cells',
                            panc$seurat_clusters %in% c(5, 24)~ 'Plasma',
                            panc$seurat_clusters == 9 ~ 'GBM crap',
                            panc$seurat_clusters %in% c(16,20,21,27,29) ~ 'Epithelial', 
                            panc$seurat_clusters %in% c(13) ~ 'Fibroblasts', 
                            TRUE ~ 'Unknown')

library(RColorBrewer)
n <- length(unique(panc$cell_type))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


DimPlot(panc,  group.by = 'cell_type_upd', label = T, cols = col_vector) 
ggsave(paste0( "Dimplot_cell_type_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

fwrite(panc@meta.data, 'Reclustered_immune_snATAC_Merged_94_sample_obj.v1.2_new_macs2_upd.metaData', row.names = T, sep = '\t')




