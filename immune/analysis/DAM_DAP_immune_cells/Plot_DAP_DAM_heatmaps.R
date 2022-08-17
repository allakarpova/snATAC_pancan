suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))
library(JASPAR2020)
library(TFBSTools)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

panc <- readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v3.0/combo_only/chromvar/Reclustered_immune_snATAC_Merged_44_combo_samples_v3.0_cluster_peaks.type.chromvar.rds')
panc@assays$chromvar@scale.data <- panc@assays$chromvar@data

combo.meta <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_immune_44_samples.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)
combo.meta$Cell_type_state <- case_when(is.na(combo.meta$Cell_type_state) | grepl('oublet', combo.meta$Cell_type_state) ~ 'Doublet',
                                        TRUE ~ combo.meta$Cell_type_state)

combo.meta$cell_type_general <- case_when(is.na(combo.meta$Cell_type_state) | grepl('oublet', combo.meta$Cell_type_state) ~ 'Doublet',
                                    grepl('CD4/CD8', combo.meta$Cell_type_state) ~ 'CD4/CD8− T−cells',
                                    grepl('CD4', combo.meta$Cell_type_state) ~ 'CD4 T-cells',
                                    grepl('CD8', combo.meta$Cell_type_state) ~ 'CD8 T-cells',
                                    grepl('NK', combo.meta$Cell_type_state) ~ 'NK cells',
                                    TRUE ~ combo.meta$Cell_type_state)

panc <- AddMetaData(panc, combo.meta[,c('Cell_type_state', 'Cell_type_markers', 'cell_type_general')])

tb <- table(panc$seurat_clusters, panc$cell_type)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
panc$cell_type1 <- cluster.match.celltype[panc$seurat_clusters]


#starts
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/v3.0/combo_only/')
panc <- readRDS('Reclustered_immune_snATAC_Merged_44_combo_samples_v3.0_cluster_peaks.type.chromvar.rds')
add_filename <- 'combo_snATAC_42_samples'
panc$cell_type.harmonized.cancer
DimPlot(panc, group.by = 'Cell_type_state', label = T) 
ggsave(paste0( "Dimplot_cell_type_", add_filename, ".pdf"),height=7,width=10,useDingbats=FALSE)

DimPlot(panc, group.by = 'cell_type', label = T) 
DimPlot(panc, label = T)

Motifs(panc)
panc@meta.data

library(future)
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 100 * 1024 ^ 3)


#DAMs for cell state 
DefaultAssay(panc) <- 'chromvar'
Idents(panc) <- 'Cell_type_state'
dams_cell_state <- FindAllMarkers(panc, test.use = 'LR', assay = 'chromvar', random.seed = 123, 
                       latent.vars = 'nCount_ATAC_immune')
setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/DAM/v3.0/combo_only')
fwrite(dams_cell_state, 'dams_cell_state.tsv', row.names = T, sep ='\t')

#DAMs for cell type general on katmai
panc <- readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v3.0/combo_only/chromvar/Reclustered_immune_snATAC_Merged_44_combo_samples_v3.0_cluster_peaks.type.chromvar.rds')
panc@assays$chromvar@scale.data <- panc@assays$chromvar@data
combo.meta <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/comboRNA_immune_cells_no_cell_cycle/manual/snRNA_combo_Merged_immune_44_samples.metadata.tsv') %>% 
  data.frame(row.names = 'Barcodes_cancer', check.rows=F, check.names=F)
panc <- AddMetaData(panc, combo.meta[,c('Cell_type_state', 'Cell_type_markers')])
panc$Cell_type_state <- case_when(is.na(panc$Cell_type_state) | grepl('oublet', panc$Cell_type_state) ~ 'Doublet',
                                  TRUE ~ panc$Cell_type_state)
panc$cell_type_general <- case_when(is.na(panc$Cell_type_state) | grepl('oublet', panc$Cell_type_state) ~ 'Doublet',
                                    grepl('CD4/CD8', panc$Cell_type_state) ~ 'CD4/CD8− T−cells',
                                    grepl('CD4', panc$Cell_type_state) ~ 'CD4 T-cells',
                                    grepl('CD8', panc$Cell_type_state) ~ 'CD8 T-cells',
                                    grepl('NK', panc$Cell_type_state) ~ 'NK cells',
                                    TRUE ~ panc$Cell_type_state)

Idents(panc) <- 'cell_type_general'
dams_cell_general <- FindAllMarkers(panc, test.use = 'LR', assay = 'chromvar', random.seed = 123, 
                                    latent.vars = 'nCount_ATAC_immune')
fwrite(dams_cell_general, '/diskmnt/Projects/snATAC_analysis/immune/DAMs/v3.0/combo_only/dams_cell_general.tsv', row.names = T, sep ='\t')

#DAPs for cell state
DefaultAssay(panc) <- 'ATAC_immune'
Idents(panc) <- 'cell_type_general'
daps_cell_general <- FindAllMarkers(panc, test.use = 'LR', assay = 'ATAC_immune', random.seed = 123, min.pct = 0.2,
                                  latent.vars = 'nCount_ATAC_immune')


fwrite(daps_cell_general, '/diskmnt/Projects/snATAC_analysis/immune/DAPs/v3.0/combo_only/daps_cell_general_min_pct_0.2.tsv', row.names = F, sep ='\t')

DefaultAssay(panc) <- 'ATAC_immune'
Idents(panc) <- 'Cell_type_state'
daps_cell_state <- FindAllMarkers(panc, test.use = 'LR', assay = 'ATAC_immune', random.seed = 123, min.pct = 0.2,
                                  latent.vars = 'nCount_ATAC_immune')
fwrite(daps_cell_state, '/diskmnt/Projects/snATAC_analysis/immune/DAPs/v3.0/combo_only/daps_cell_state_min_pct_0.2.tsv', row.names = F, sep ='\t')



############ PLOOOOOOOT 
Idents(panc) <- 'Cell_type_state'

panc@assays$chromvar@counts <- panc@assays$chromvar@data
panc.aver.state <- AverageExpression(panc, slot = 'counts', assays = c('chromvar'),return.seurat = T)
Idents(panc) <- 'cell_type_general'
panc.aver.general <- AverageExpression(panc, slot = 'counts', assays = c( 'chromvar'), return.seurat = T)

dams_cell_state


setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/DAM/v3.0/combo_only')
dams_cell_state <- fread('dams_cell_state.tsv')
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
motif.tf <- data.frame(motif = names(pfm@listData),
           TF = map_chr(names(pfm@listData), ~pfm@listData[[.x]]@name))

rm(pfm)
dams_cell_state$TF <- motif.tf$TF[match(dams_cell_state$gene, motif.tf$motif)]

genes.oi <- (dams_cell_state %>% filter(avg_log2FC>2 & p_val_adj < 0.01) %>% pull(gene)) %>% unique
dont.show <- dams_cell_state %>% filter(avg_log2FC>2 & p_val_adj < 0.01) %>% 
  group_by(gene) %>% 
  tally %>% 
  arrange(-n) %>% 
  filter(n>4) %>%
  pull(gene)


DefaultAssay(panc.aver.state) <- 'chromvar'
toplot <- FetchData(panc.aver.state, vars = genes.oi, slot = 'counts')
colnames(toplot) <- motif.tf$TF[match(colnames(toplot), motif.tf$motif)]

p <- pheatmap::pheatmap (toplot,
                   cellwidth = 6, 
                   cellheight =6.5,
                   cutree_cols = 12,
                   cutree_rows = 4,
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
pdf ('Heatmap_average_motif_scores_dams_cell_state.pdf', width = 30, height = 10)
print(p)
dev.off()

DefaultAssay(panc.aver.state) <- 'chromvar'
toplot <- FetchData(panc.aver.state, vars = setdiff(genes.oi, dont.show), slot = 'counts')
colnames(toplot) <- motif.tf$TF[match(colnames(toplot), motif.tf$motif)]
p <- pheatmap::pheatmap (toplot,
                         cellwidth = 6, 
                         cellheight =6.5,
                         cutree_cols = 12,
                         cutree_rows = 4,
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
pdf ('Heatmap_average_motif_scores_dams_cell_state_no_seen_in4.pdf', width = 20, height = 7)
print(p)
dev.off()

dams <- dams_cell_state %>% filter(p_val_adj < 0.01) %>% pull(gene) %>% unique
DefaultAssay(panc)

tocor <- FetchData(panc.aver.state, vars = dams, slot = 'counts') %>% t
cor.result <- psych::corr.test(tocor, method = 'spearman',adjust = 'fdr')

p <- pheatmap::pheatmap (cor.result$r,
                    cellwidth = 6, 
                    cellheight = 6,
                    #cutree_cols = 12,
                    #cutree_rows = 4,
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
pdf ('Cor_Heatmap_average_motif_scores_dams_cell_state.pdf', width = 5, height = 5)
print(p)
dev.off()


# now do DAMs on general cell types
dams_cell_general <- fread('dams_cell_general.tsv')
dams_cell_general$TF <- motif.tf$TF[match(dams_cell_general$gene, motif.tf$motif)]
genes.oi <- (dams_cell_general %>% filter(avg_log2FC>2 & p_val_adj < 0.01) %>% pull(gene)) %>% unique
DefaultAssay(panc.aver.general) <- 'chromvar'
toplot <- FetchData(panc.aver.general, vars = genes.oi, slot = 'counts')
colnames(toplot) <- motif.tf$TF[match(colnames(toplot), motif.tf$motif)]
p <- pheatmap::pheatmap (toplot,
                         cellwidth = 6, 
                         cellheight =6.5,
                         cutree_cols = 12,
                         cutree_rows = 4,
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
pdf ('Heatmap_average_motif_scores_dams_cell_general.pdf', width = 30, height = 10)
print(p)
dev.off()

genes.oi <- (dams_cell_general %>% filter(avg_log2FC>2 & p_val_adj < 0.01) %>% pull(gene)) %>% unique
dont.show <- dams_cell_general %>% filter(avg_log2FC>2 & p_val_adj < 0.01) %>% 
  group_by(gene) %>% 
  tally %>% 
  arrange(-n) %>% 
  filter(n>4) %>%
  pull(gene)
DefaultAssay(panc.aver.general) <- 'chromvar'
toplot <- FetchData(panc.aver.general, vars = genes.oi, slot = 'counts')
colnames(toplot) <- motif.tf$TF[match(colnames(toplot), motif.tf$motif)]
p <- pheatmap::pheatmap (toplot,
                         cellwidth = 6, 
                         cellheight =6.5,
                         cutree_cols = 12,
                         cutree_rows = 4,
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
pdf ('Heatmap_average_motif_scores_dams_cell_general_no_seen_in4.pdf', width = 20, height = 6)
print(p)
dev.off()
getwd()
(dams_cell_general %>% filter(avg_log2FC>2 & p_val_adj < 0.01) %>% group_by(gene)) %>% tally %>% arrange(-n) %>% filter(n>2)

### now DO DAPs
panc$Cell_type_state <- case_when(is.na(panc$Cell_type_state) | grepl('oublet', panc$Cell_type_state) ~ 'Doublet',
                                  TRUE ~ panc$Cell_type_state)
panc$cell_type_general <- case_when(is.na(panc$Cell_type_state) | grepl('oublet', panc$Cell_type_state) ~ 'Doublet',
                                    grepl('CD4/CD8', panc$Cell_type_state) ~ 'CD4/CD8− T−cells',
                                    grepl('CD4', panc$Cell_type_state) ~ 'CD4 T-cells',
                                    grepl('CD8', panc$Cell_type_state) ~ 'CD8 T-cells',
                                    grepl('NK', panc$Cell_type_state) ~ 'NK cells',
                                    TRUE ~ panc$Cell_type_state)
setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/DAP/v3.0/combo_only')
daps_cell_state <- fread('daps_cell_general_min_pct_0.2.tsv')
DefaultAssay(panc) <- 'ATAC_immune'
Annotation(panc)

genes.oi <- daps_cell_state %>% filter(p_val_adj < 0.01) %>% pull(gene) %>% unique
genes.oi
Idents(panc) <- 'cell_type_general'
panc.aver.general <- AverageExpression(panc, slot = 'data', assays = c( 'ATAC_immune'), return.seurat = T, features =genes.oi)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Signac)
peakAnno.important <- annotatePeak(StringToGRanges(genes.oi), tssRegion=c(-1500, 500),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')
peakAnno.important <- as.data.table(peakAnno.important@anno)
peakAnno.important$peak <- paste(peakAnno.important$seqnames, peakAnno.important$start, peakAnno.important$end, sep = '-')
peakAnno.important[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]


rm(genes.oi2)
genes.oi2 <- daps_cell_state %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.01) %>% 
  group_by(cluster) %>%
  top_n(100, wt = avg_log2FC) %>% 
  pull(gene) %>% unique
toplot <- FetchData(panc.aver.general, vars = genes.oi2, slot = 'scale.data')

p <- pheatmap::pheatmap (toplot, scale = 'none',
                         cellwidth = 0.5, 
                         cellheight =6.5,
                         #cutree_cols = 12,
                         #cutree_rows = 4,
                         #color=color.palette,
                         na_col = 'white',
                         cluster_rows = T,
                         cluster_cols = T, 
                         clustering_distance_rows = "euclidean", 
                         clustering_distance_cols = "euclidean", 
                         clustering_method = "ward.D2", 
                         border_color = NA, 
                         #breaks=breaks, 
                         fontsize_col  = 7,
                         fontsize_row  = 7,
                         show_rownames= T,
                         show_colnames = F,
                         fontsize_number = 5)

pdf ('Heatmap_average_DAPs_first100_FC0.5_clusters.pdf', width = 10, height = 10)
print(p)
dev.off()

df <- peakAnno.important %>% filter(peak %in% genes.oi2)

ggplot(data = df,  aes(x=peak.position, fill=peak.position)) +
  geom_bar(stat = 'count', show.legend = F) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('DAPs between immune cell types')
ggsave('Bar_stack_important_nadya_in_my_ccRCC_nonover400.pdf' , useDingbats = F, width = 4.5, height = 6)


ct <- daps_cell_state$cluster%>% unique
ct %>% walk(function(c) {
  daps.Macrophages <- daps_cell_state %>% filter(cluster == c & avg_log2FC>0  & p_val_adj < 0.01) %>% pull(gene) %>% unique
  peakAnno.important.TAMS <- peakAnno.important %>% filter(peak %in% daps.Macrophages)
  
  ggplot(data = peakAnno.important.TAMS,  aes(x=peak.position, fill=peak.position)) +
    geom_bar(stat = 'count', show.legend = F) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(color = 'grey')) +
    scale_fill_brewer(palette = 'Set1') +
    ggtitle(paste('DAPs in',c, 'vs others'))
  ggsave(paste0('Bar_DAP_position_', make.names(c),'.pdf') , useDingbats = F, width = 4.5, height = 6)
  
})

ct
toplot <- daps_cell_state %>% filter(cluster == 'mregDC' & avg_log2FC>0.5  & p_val_adj < 0.01)
toplot <- merge(toplot, peakAnno.important[,c('peak', 'SYMBOL', 'peak.position')], by.x = 'gene', by.y = 'peak', all.x = T)

ggplot(data = toplot,  aes(x=peak.position,y=SYMBOL, color=avg_log2FC)) +
  geom_point(size = 4) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  viridis::scale_color_viridis( direction = 1, limits = c(0,4)) +
  ggtitle(paste('DAPs in','mregDC', 'vs others'))
ggsave(paste0('Dotplot_DAP_position_', make.names(c),'.pdf') , useDingbats = F, width = 4.5, height = 6)



t.cell.peaks <- daps_cell_state %>% filter(cluster %in% c('CD4 T-cells', 'CD8 T-cells', 'NK cells') & avg_log2FC>0  & p_val_adj < 0.01) %>% group_by(gene) %>% tally %>% arrange(-n) %>% filter (n==3) %>% pull(gene)
toplot <- daps_cell_state %>% filter(cluster %in% c('CD4 T-cells', 'CD8 T-cells', 'NK cells') & avg_log2FC>0  & p_val_adj < 0.01 & gene %in% t.cell.peaks)
toplot <- merge(toplot, peakAnno.important[,c('peak', 'SYMBOL', 'peak.position')], by.x = 'gene', by.y = 'peak', all.x = T)
ggplot(data = toplot,  aes(x=peak.position,y=SYMBOL, color=avg_log2FC)) +
  geom_point(size = 4) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  viridis::scale_color_viridis( direction = 1, limits = c(0,1.5)) +
  facet_wrap(~cluster, scales = 'free') +
  ggtitle(paste('DAPs commons for T-cells and NK cells vs others')) 
ggsave(paste0('Dotplot_DAP_position_', make.names(c),'.pdf') , useDingbats = F, width = 4.5, height = 6)

