# annotate CRC objects cell types
library(Seurat)
library(ggplot2)
library(viridis)
library(data.table)
library(dplyr)

sample <- 'CM618'
obj <- paste0('/diskmnt/Projects/MetNet_analysis/Colorectal/Seurat/Integration/',sample,'.rds')
out.dir <- paste0('/diskmnt/Projects/snATAC_primary/03_paired_sn_rna/CRC/cell_typing/', sample)
dir.create(out.dir)
setwd(out.dir)
panc <- readRDS(obj)

DimPlot(panc, group.by = 'seurat_clusters', label = T)
ggsave(paste0('Dimplot_clusters_', sample,'.pdf'), width = 6, height = 5, useDingbats = F)

#panc <- FindClusters(panc, resolution = 1)
#DimPlot(panc, group.by = 'SCT_snn_res.1', label = T)
#ggsave(paste0('Dimplot_SCT_snn_res.1_', sample,'.pdf'), width = 6, height = 5, useDingbats = F)

DimPlot(panc, group.by = 'cell_type', label = T)
ggsave(paste0('Dimplot_cell_type_', sample,'.pdf'), width = 6.5, height = 5, useDingbats = F)

degs.long.set <- fread('/diskmnt/Projects/Users/allakarpova/Data/markers/Liver_mouse_markers.txt', data.table = F)
degs.long.set$`Gene symbol` <- toupper(degs.long.set$`Gene symbol`)
genes2plot <- degs.long.set$`Gene symbol` %>% unique() %>% toupper()
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
ggsave(paste0( "Dotplot_20mouse_liver_paper_markers_gene_expression_", sample, "_many_fromRNA.pdf"),height=20,width=20,useDingbats=FALSE,limitsize = FALSE)


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
p <- p + scale_color_viridis_c("viridis", direction = 1)
p
ggsave(paste0( "Dotplot_reyka_markers_gene_expression_", sample, "_many_fromRNA.pdf"),height=20,width=20,useDingbats=FALSE,limitsize = FALSE)


########### MY MARKERS
myeloid.genes <- fread('/diskmnt/Projects/Users/allakarpova/Data/markers/Cell_state_markers_v12182020.txt', data.table = F, header = T)
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
ggsave(paste0( "Dotplot_marker_gene_expression_", sample, "_chromatin.pdf"),height=50,width=50,useDingbats=FALSE,limitsize = FALSE)

##### for CM556
Idents(panc) <- 'SCT_snn_res.0.5'
panc <- RenameIdents(object = panc, 
                   `0` = "Tumor", 
                   `1` = "Macrophage", 
                   `2` = "Tumor", #high TRPM3, CASC15
                   `3` = "T",
                   `4` = "T",
                   `5` = "Fibroblast",
                   `6` = "Macrophage",
                   `7` = "Plasma",
                   `8` = "T",
                   `9` = "Plasma",  
                   `10` = "Portal_Endothelial", 
                   `11` = "Tumor",
                   `12` = "B", #Subset is mregDC
                   `13` = "Fibroblast", #specific to adrenal
                   `14` = "Macrophage",
                   `15` = "Oligodendrocyte",
                   `16` = "Fibroblast",
                   `17` = "Unknown", #liver specific
                   `18` = "Tumor",
                   `19` = "Hepatocyte",
                   `20` = "Cholangiocyte",
                   `21` = "Neuron",
                   `22` = "Fasiculata"
)
DimPlot(panc, label = T)
ggsave(paste0('Dimplot_reyka_cell_types_', sample,'.pdf'), width = 6, height = 5, useDingbats = F)

panc$cell_type = Idents(panc)
saveRDS(panc, paste0(sample, '_cell_typed.rds'))

### for CM618
Idents(panc) <- 'SCT_snn_res.0.5'
s1 <- RenameIdents(object = s1, 
                   `0` = "Tumor", 
                   `1` = "Tumor", 
                   `2` = "T",
                   `3` = "Kupffer",
                   `4` = "Tumor",
                   `5` = "Fibroblast",
                   `6` = "Tumor",
                   `7` = "Tumor",
                   `8` = "Tumor",
                   `9` = "Zone_1_Periportal_Hepatocyte",  #subset of 9  is cholangio
                   `10` = "Tumor", 
                   `11` = "Endothelial",
                   `12` = "Intestinal_Goblet",
                   `13` = "Unknown",
                   #`14` = "Periportal_LSEC",
                   `15` = "Plasma",
                   `16` = "U3",
                   `17` = "pDC",
                   `18` = "Fibroblast",
                   `LSEC1` = "periportal_LSEC",
                   `LSEC2` = "portal_endothelial",
                   `LSEC3` = "central_venous_LSEC"
)

s1 <- RenameIdents(object = s1, 
                   `DC` = "mregDC") 
s1 <- RenameIdents(object = s1, 
                   `U1` = "Intestinal_Goblet") 
s1 <- RenameIdents(object = s1, 
                   `U2` = "Other",
                   `U4`="pDC",
                   `U5`="Fibroblast"
) 
