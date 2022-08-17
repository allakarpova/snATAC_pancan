library(Signac)
library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(paletteer)
library(readxl)
library(ggplot2)

### FUNCTIONS 
Caps <- function(x) {
  s <- x
  toreturn = sapply(s, function(x) paste(toupper(substring(x, 1,1)), substring(x, 2),
                                         sep="") )
  return(as.character(toreturn))
}
##########

panc <- readRDS ('~/lab_Ding/work/single_cell/senescence/snRNA_combo/objects/mouse_liver/2_snRNA_Merged_mouse_liver_young_old.rds')

wd <- '~/lab_Ding/work/single_cell/senescence/snRNA_combo/cell-annotation'
dir.create(wd)
setwd(wd)
add_filename <- 'mouse_liver'

# plot dotplot markers MINE
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v12182020.txt', data.table = F, header = T)
myeloid.genes$Gene <- myeloid.genes$Gene %>% tolower() %>% Caps()

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
p
ggsave(paste0( "Dotplot_marker_gene_expression_", add_filename, "_chromatin.pdf"),height=50,width=50,useDingbats=FALSE,limitsize = FALSE)


# ### use markers from mouse liver paper https://www.liebertpub.com/doi/10.1089/hum.2019.366 ######
#######
read_excel('~/lab_Ding/work/single_cell/senescence/snATAC_combo/papers/Supp_Table-S1-converted.xlsx')
degs <- sapply(seq(1,11), function (x) {
  read_excel('~/lab_Ding/work/single_cell/senescence/snATAC_combo/papers/Supp_Table-S1-converted.xlsx', sheet = x,trim_ws = T, col_names = T) })

names(degs) <- c('Hep1', 'Hep2', 'Hep2too', 'Hep3', 'Hep3too', 'Endothelial', 'hepatic stellate cells' ,'Kupffer cells', 'B cells', 'T cells', 'cholangiocytes')
degs$Hep2 <- degs$Hep2[-2:-3]
degs$Hep3 <- degs$Hep3[-2]
degs <- lapply(names(degs), function (x) {
  y <- degs[[x]]
  colnames(y) <- c("Gene symbol", "P value",  "Average difference", "cells.pct")
  y$Cell.type <- x
  return(y)
})

degs.long <- do.call ('rbind', degs)
degs.long <- degs.long %>% dplyr::filter (!is.na(`P value`))
degs.long.set <- degs.long %>% 
  dplyr::filter(!grepl('too', Cell.type)) %>%
  group_by(Cell.type) %>%
  top_n(n = 20, wt = `Average difference`)

fwrite(degs.long.set, '~/lab_Ding/work/single_cell/Liver_mouse_markers.txt', sep = '\t')
degs.long.set %>% group_by(Cell.type) %>% tally

genes2plot <- degs.long.set$`Gene symbol` %>% unique()
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA' )
p$data <- merge(p$data, degs.long.set[c(1,5)], by.x = 'features.plot', by.y = 'Gene symbol', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell.type , scales = "free",  drop = T, ncol = 3)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_viridis_c("viridis", direction = 1)
p
ggsave(paste0( "Dotplot_20mouse_liver_paper_markers_gene_expression_", add_filename, "_many_fromRNA.pdf"),height=20,width=20,useDingbats=FALSE,limitsize = FALSE)

#### REYKA markers
Hepatocytes=c("ALB","CYP3A7","CYP2D6","CYP2A7","CYP2A6","ARG1","GPC3","SCD","HMGCS1","ACSS2","TM7SF2","GSTA2","AKR1C1")
sinusoidal=c("F8","PECAM1","LYVE1","STAB2","CD34","ENG","MCAM","ICAM1","VWF","STAB1","COL1A1","CD40","CD80","AOC3")
Cholangiocytes=c("SOX9","EPCAM","KRT19","KRT7","CFTR","NPBF3","PKD2")
stellate=c("ACTA2","COL1A1","COL1A2","COL3A1","RBP1","TAGLN","SPARC","DES","SMN1","GFAP","NTRK3","RBP1")
portal_endothelial=c("RAMP3","INMT","DNASEIL3","LIFR")
LSEC=c("MGP","SPARCL1","TM4SF1","CLEC14A","CCL14","CLEC1B","FCN2","S100A13")
inflamm_macs=c("S100A8","S100A9","LYZ","HLA-DBP1")
noninflamm_macs=c("CD5L","MARCO","VSIG4")
fasiculata=c("FDX1","CYP11B1","CYP11A1","MEG3","MGARP1")

markers <- rbind(data.frame(Gene = Hepatocytes, Cell.type = 'Hepatocytes'),
                 data.frame(Gene = sinusoidal, Cell.type = 'Sinusoidal'),
                 data.frame(Gene = Cholangiocytes, Cell.type = 'Cholangiocytes'),
                 data.frame(Gene = stellate, Cell.type = 'Stellate'),
                 data.frame(Gene = portal_endothelial, Cell.type = 'Portal endothelial'),
                 data.frame(Gene = LSEC, Cell.type = 'LSEC'),
                 data.frame(Gene = fasiculata, Cell.type = 'Fasiculata'),
                 data.frame(Gene = inflamm_macs, Cell.type = 'inflamm_macs'),
                 data.frame(Gene = noninflamm_macs, Cell.type = 'noninflamm_macs'))
markers$Gene <- markers$Gene %>%tolower() %>% Caps()
  
genes2plot <- markers$Gene %>% unique()
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA' )
p$data <- merge(p$data, markers, by.x = 'features.plot', by.y = 'Gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell.type , scales = "free",  drop = T, ncol = 3)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0( "Dotplot_reyka_markers_gene_expression_", add_filename, "_many_fromRNA.pdf"),height=20,width=20,useDingbats=FALSE,limitsize = FALSE)


FeaturePlot(panc, 'Gzmk')

# change back to working with peaks instead of gene activities
DefaultAssay(panc) <- 'peaksinters'

degs.0 <- FindMarkers(
  object = panc,
  ident.1 = "0",
  min.pct = 0.1,
  test.use = 'LR'
)

VlnPlot(
  object = panc,
  features = rownames(subset(degs.0, avg_logFC > 0))[1:10],
  pt.size = 0.1,
  #idents = c("CD4 Memory","CD14 Mono")
)

FeaturePlot(
  object = panc,
  features = rownames(da_peaks)[4],
  pt.size = 0.1
)

## add doublet info
scrublet01 <- fread('~/lab_Ding/work/single_cell/senescence/snATAC_combo/Scrublet/ATAC.only/SM001H1-Md/SM001H1-Md_scrublet_output_table.csv', data.table = F)
rownames(scrublet01) <- paste('SM001H1-Md', scrublet01$Barcodes, sep = '_')
scrublet02 <- fread('~/lab_Ding/work/single_cell/senescence/snATAC_combo/Scrublet/ATAC.only/SM002H1-Md/SM002H1-Md_scrublet_output_table.csv', data.table = F)
rownames(scrublet02) <- paste('SM002H1-Md', scrublet02$Barcodes, sep = '_')

panc$Predicted_doublet <- case_when(rownames(panc@meta.data) %in% rownames(subset(scrublet01, predicted_doublet)) | rownames(panc@meta.data) %in% rownames(subset(scrublet02, predicted_doublet)) ~ 'Doublet',
                                    TRUE ~ 'Singlet')
DimPlot(panc, group.by = 'Predicted_doublet', label = T)
rownames(panc@meta.data)
# ### DEFINE CELL TYPES ATAC 
# panc$Cell_type <- case_when( panc$seurat_clusters %in% c('3','7','6', '15') ~ 'Endothelial cells',
#                              panc$seurat_clusters %in% c('5','2', '8','13') ~ 'Hepathocytes',
#                              panc$seurat_clusters %in% c('4') ~ 'Hepatic stellate cells',
#                              panc$seurat_clusters %in% c('1','9', '16') ~ 'Kuppfer cells',
#                              panc$seurat_clusters %in% c('12') ~ 'CD8 T-cells',
#                              panc$seurat_clusters %in% c( '10') ~ 'B-cells',
#                              panc$seurat_clusters %in% c('11') ~ 'Immune cells',
#                              panc$seurat_clusters %in% c('14') ~ 'Cholangiocytes',
#                              panc$seurat_clusters %in% c('0') ~ 'Hepathocytes unusual')

### DEFINE CELL TYPES RNA object
panc$Cell_type <- case_when( panc$seurat_clusters %in% c('2','3','17') ~ 'Endothelial cells',
                             panc$seurat_clusters %in% c('1','4', '5','6','7','9','10','14','13','20') ~ 'Hepathocytes',
                             panc$seurat_clusters %in% c('8', '11') ~ 'Hepatic stellate cells',
                             panc$seurat_clusters %in% c('0','19', '16') ~ 'Kuppfer cells',
                             panc$seurat_clusters %in% c('15') ~ 'CD8 T-cells',
                             panc$seurat_clusters %in% c( '12') ~ 'B-cells',
                             panc$seurat_clusters %in% c('18') ~ 'Cholangiocytes')
panc$Age <- ifelse (panc$orig.ident=='SM001H1-Md', 'Old', 'Young')
DefaultAssay(panc) <- 'peaksinters'


DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave('Dimplot_clusters_macs2.pdf', width = 6, height =5, useDingbats = F)

DimPlot(panc,split.by = 'Age', group.by = 'Cell_type', label = T, label.size = 2)
ggsave('Dimplot_cell_type_split_by_age.pdf', width = 12, height =5, useDingbats = F)

DimPlot(panc,group.by = 'Age', label = T, label.size = 2, cols = 'Paired')
ggsave('Dimplot_age.pdf', width = 6, height =5, useDingbats = F)

DimPlot(panc, group.by = 'Cell_type', label = T)
ggsave('Dimplot_cell_type.pdf', width = 6, height =5, useDingbats = F)

panc$Cell_type_macs2 <- case_when( panc$peaksMACS2_snn_res.0.8 %in% c('1','4','8') ~ 'Endothelial cells',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('5','6', '9','7') ~ 'Hepathocytes',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('3') ~ 'Hepatic stellate cells',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('2','15', '11') ~ 'Kuppfer cells',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('13') ~ 'CD8 T-cells',
                                   panc$peaksMACS2_snn_res.0.8 %in% c( '10') ~ 'B-cells',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('12') ~ 'Myeloid cells',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('14') ~ 'Cholangiocytes',
                                   panc$peaksMACS2_snn_res.0.8 %in% c('0') ~ 'Hepathocytes unusual')


DimPlot(panc,split.by = 'Age', group.by = 'Cell_type_macs2', label = T, label.size = 2)
ggsave('Dimplot_cell_type_macs2_split_by_age.pdf', width = 12, height =5, useDingbats = F)
DimPlot(panc,group.by = 'Cell_type_macs2', label = T, label.size = 2)
ggsave('Dimplot_cell_type_macs2.pdf', width = 6, height =5, useDingbats = F)

panc@assays$chromvar@data

