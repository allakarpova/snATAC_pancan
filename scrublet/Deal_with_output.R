
library(readxl)
library(plyr)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)
library(paletteer)
library(cowplot)


setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_v2/')
df.res <- list.files(path = './', pattern = 'output*', recursive = T, full.names = F)


df.meta.total <- lapply(df.res, FUN = function (x) {
  sample <- str_split_fixed(x, '_scrub', 2)[,1]
  df.metadata <- fread(x)
  df.metadata$Sample <- sample
  df.metadata$merged_barcode <- paste(sample, df.metadata$Barcodes, sep ='_')
  
  df.metadata
  return (df.metadata)
})
df.meta.total <- rbindlist(df.meta.total)

df.meta.total$Sample %>% unique %>% sort
df.meta.total <- df.meta.total %>% 
  mutate(predicted.doublet.upd = case_when(Sample %in% c('CRC_CM354C1-T1Y2N1Bma1_1', 'CRC_CM354C2-T1Y2N1Bma1_1', 'CRC_HT250C1-Th1K1Fc2A2N1_1Bma1_1', 'PDAC_HT270P1-S1H2Fc2A2N1Bma1_1') & doublet_score >= 0.2 ~ 'True',
                                           Sample %in% c('OV_VF027V1-S1Y1N1_1Bma1_1') & doublet_score >= 0.22 ~ 'True',
                                           Sample %in% c('BRCA_HT243B1-S1H4Fc2A2N1Bma1_1', 'PDAC_HT288P1-S1H4Fc2A2N1Bma1_1') & doublet_score >= 0.24 ~ 'True',
                                           Sample %in% c('BRCA_HT268B1-Th1H3Fc2A2N1Bma1_1', 'BRCA_HT271B1-S1H3Fc2A5N1Bma1_1', 'BRCA_HT297B1-S1H1Fc2A2N1_1Bma1', 'CRC_HT307C1-Th1K1Fc2A2N1Bma1_1', 'PDAC_HT242P1-S1H1Fc2A2N1Bma1_1', 'PDAC_HT264P1-S1H2Fc2A2N1Bma1_1', 'UCEC_CPT3936DU-T1N1_1Bma1_1') & doublet_score >= 0.26 ~ 'True',
                                           Sample %in% c('CRC_CM478C2-T1Y2N1Bma1_1', 'CRC_CM1563C1-T1Y1N1Bma1_1','CRC_HT225C1-Th1Fc1A2N1Bma1_1', 'CRC_HT253C1-Th1K1Fc2A2N1Bma1_1', 'HNSCC_P5539-N1_1Bma1_1', 'PDAC_HT224P1-S1Fc2A2N1Bma1_1') & doublet_score >= 0.28 ~ 'True',
                                           Sample %in% c('BRCA_HT235B1-S1H1Fc2A2N1Bma1_1','BRCA_HT263B1-S1H1A3N1Ba1_1', 'BRCA_HT305B1-S1H1Fc2A2_1N1Bma1_1', 'HNSCC_P5504-N1_1Bma1_1', 'PDAC_HT231P1-S1H3Fc2A2N1Bma1_1', 'PDAC_HT232P1-S1H1Fc2A2_1N1Bma1_1', 'PDAC_HT306P1-S1H1Fc2A2N1Bma1_1', 'UCEC_CPT704DU-T1N1_1Bma1_1', 'UCEC_CPT1541DU-S1N1Bma1_1') & doublet_score >= 0.3 ~ 'True',
                                           Sample %in% c('UCEC_CPT1541DU-T1N1_1Bma1_1', 'UCEC_CPT2552DU-S1N1Z1_1Bma1_1', 'UCEC_CPT3936DU-S1N1_1Bma1_1')  & doublet_score >= 0.32 ~ 'True',
                                           Sample %in% c('CRC_HT260C1-Th1K1Fc2A2N1Bma1_1', 'CRC_SP369H1-Mc1N1Bma1_1', 'CRC_SP819H1-Mc1_1N1Bma1_1', 'OV_VF034V1-T1Y1N1_1Bma1_1', 'PDAC_HT259P1-S1H1Fc2A2_1N1Bma1_1', 'UCEC_CPT2552DU-T1N1Z1_1Bma1_1') & doublet_score >= 0.34 ~ 'True',
                                           Sample %in% c('CRC_CM663C1-T1Y1N1Bma1_1') & doublet_score >= 0.36 ~ 'True',
                                           Sample %in% c('CRC_CM478C1-T1Y2N1Bma1_1') & doublet_score >= 0.38 ~ 'True',
                                           Sample %in% c('CRC_CM1563C1-S1Y1N1Bma1_1', 'CRC_HT230C1-Th1Fc2A2N1Bma1_1', 'UCEC_CPT704DU-S1N1_1Bma1_1') & doublet_score >= 0.4 ~ 'True',
                                           Sample %in% c('CRC_HT291C1-M1A3N1Bma1_1') & doublet_score >= 0.42 ~ 'True',
                                           TRUE ~ 'False'))

fwrite(df.meta.total, '/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_v2/Scrublet_scores_adjusted.cutoff.tsv', sep = '\t')
fwrite(df.meta.total %>% 
         filter(grepl('BRCA', Sample)) %>% 
         mutate(Piece_ID = str_split_fixed(string = gsub(pattern = 'BRCA_', replacement = '', Sample), pattern = 'Fc', 2)[,1]) %>%
         mutate(Piece_ID = str_split_fixed( Piece_ID,'A3', 2)[,1]), '~/lab_Ding/work/single_cell/htan_brca/revision/metadata/Scrublet_BRCA_scores_adjusted.cutoff.tsv', sep = '\t')

#compare scrublet runs with different parameters
rna.atac.all <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_v2/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_scrublet_output_table.csv')
rna <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_RNAonly/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_scrublet_output_table.csv')
rna.log <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_RNAonly_withlog/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_scrublet_output_table.csv')
rna.log.mean <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_RNAonly_withlog_mean/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_scrublet_output_table.csv')
rna.all <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_RNAonly_withlog_mean_var/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_scrublet_output_table.csv')
atac.log <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboATAC_withlog/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_scrublet_output_table.csv')

table(rna.atac.all$predicted_doublet, rna.all$predicted_doublet)
table(rna.all$predicted_doublet, atac.log$predicted_doublet)


panc <- readRDS('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/BRCA/signac/v4.0/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1/BRCA_HT235B1-S1H1Fc2A2N1Bma1_1_processed_multiomic.rds')
panc <- NormalizeData(panc)
DimPlot(panc, group_by = 'predicted_doublet')
rna.atac.all <- data.frame(rna.atac.all,  row.names = 1, check.rows = F, check.names = F) %>% 
  mutate(predicted_doublet.upd = case_when(doublet_score >= 0.3 ~ 'True',
                                            TRUE ~ 'False'))

panc <- AddMetaData(panc, rna.atac.all)
rna.all <- data.frame(rna.all,  row.names = 1, check.rows = F, check.names = F) %>% 
  mutate(predicted_doublet.upd = case_when(doublet_score >= 0.3 ~ 'True',
                                           TRUE ~ 'False'))
colnames(rna.all) <- paste(colnames(rna.all) , 'rna.only', sep='-')
panc <- AddMetaData(panc, rna.all)
p <- DimPlot(panc, group.by = 'predicted_doublet', cols = c('black', 'yellow'), reduction = 'umap.rna') + DimPlot(panc, group.by = 'predicted_doublet.upd', cols = c('black', 'yellow'), reduction = 'umap.rna')

p1 <- DimPlot(panc, group.by = 'predicted_doublet.rna.only', cols = c('black', 'yellow'), reduction = 'umap.rna') + DimPlot(panc, group.by = 'predicted_doublet.upd.rna.only', cols = c('black', 'yellow'), reduction = 'umap.rna')
p/p1
getwd()

atac.log <- data.frame(atac.log,  row.names = 1, check.rows = F, check.names = F) %>% 
  mutate(predicted_doublet.upd = case_when(doublet_score >= 0.32 ~ 'True',
                                           TRUE ~ 'False'))
colnames(atac.log) <- paste(colnames(atac.log) , 'atac.only', sep='-')
panc <- AddMetaData(panc, atac.log)
pa <- DimPlot(panc, group.by = 'predicted_doublet.atac.only', cols = c('grey', '#ffd92f'), reduction = 'umap.atac') + 
  DimPlot(panc, group.by = 'predicted_doublet.upd.atac.only', cols = c('grey', '#ffd92f'), reduction = 'umap.atac')
pr <- DimPlot(panc, group.by = 'predicted_doublet.rna.only', cols = c('grey', '#ffd92f'), reduction = 'umap.rna') + 
  DimPlot(panc, group.by = 'predicted_doublet.upd.rna.only', cols = c('grey', '#ffd92f'), reduction = 'umap.rna')
pa/pr
ggsave('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/HT235B1_dimplots_doublets_compare.pdf', width = 10, height = 7, useDingbats = F)

table(panc$predicted_doublet.rna.only, panc$predicted_doublet.atac.only)
table(panc$predicted_doublet.upd.rna.only, panc$predicted_doublet.upd.atac.only)

panc$doublets.rna.atac <- case_when(panc$predicted_doublet.upd.rna.only == 'True' & panc$predicted_doublet.upd.atac.only == 'True' ~ 'Both',
                                    panc$predicted_doublet.upd.rna.only == 'True' ~ 'Rna',
                                    panc$predicted_doublet.upd.atac.only == 'True' ~ 'Atac',
                                    TRUE ~ 'Singlet')
pra <- DimPlot(panc, reduction = 'wnn.umap', 
               cols.highlight = c('#66c2a5','#ffd92f', '#e78ac3', 'grey'), sizes.highlight = 0.4,
               cells.highlight = list(Both = rownames(subset(panc@meta.data, doublets.rna.atac=='Both')), 
                                      Atac = rownames(subset(panc@meta.data, doublets.rna.atac=='Atac')),
                                      Rna = rownames(subset(panc@meta.data, doublets.rna.atac=='Rna'))))
pra
ggsave('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/HT235B1_dimplots_doublets_intersection.pdf', width = 5, height = 4, useDingbats = F)

FeaturePlot(panc, features = c('doublet_score.atac.only', 'doublet_score.rna.only'), order = T, reduction = 'wnn.umap')
ggsave('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/HT235B1_featureplot_doublets_scores.pdf', width = 8, height = 4, useDingbats = F)

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

ggsave(paste0( "Dotplot_marker_gene_expression_", '235B1', "_RNA.pdf"),plot=p,height=50,width=50,useDingbats=FALSE,limitsize = FALSE)

DimPlot(panc, label = T, reduction = 'umap.rna')
panc$cell_type <- case_when(panc$seurat_clusters == 11 ~ 'Fibroblasts',
                            panc$seurat_clusters == 7 ~ 'Macrophages',
                            panc$seurat_clusters == 9 ~ 'T-cells',
                            panc$seurat_clusters == 12 ~ 'Plasma',
                            TRUE ~ 'Tumor')
VlnPlot(panc, assay = 'RNA', features = 'nFeature_RNA', split.by = 'predicted_doublet.rna.only',group.by  = 'cell_type') +
VlnPlot(panc, assay = 'RNA', features = 'nFeature_RNA', split.by = 'predicted_doublet.upd.rna.only',group.by  = 'cell_type')


VlnPlot(panc, assay = 'RNA', features = 'atac_fragments', split.by = 'predicted_doublet',group.by  = 'cell_type') +
  VlnPlot(panc, assay = 'RNA', features = 'atac_fragments', split.by = 'predicted_doublet.upd',group.by  = 'cell_type')

panc@meta.data


##### take care of regular ATAC scrublet results
setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/regATAC_15p_rate/')
df.res <- list.files(path = './', pattern = 'output*', recursive = T, full.names = F)
df.meta.total <- lapply(df.res, FUN = function (x) {
  sample <- str_split_fixed(x, '_scrub', 2)[,1]
  df.metadata <- fread(x)
  df.metadata$Sample <- sample
  df.metadata$merged_barcode <- paste(sample, df.metadata$Barcodes, sep ='_')
  df.metadata
  return (df.metadata)
})
df.meta.total <- rbindlist(df.meta.total)
fwrite(data.table(Sample = df.meta.total$Sample %>% unique %>% sort, Cutoff = 0.3), 'Sample_cutoffs.tsv', sep = '\t')

adj.cut <- fread('Sample_cutoffs_adj.txt', data.table = F)
df.meta.total <- merge(df.meta.total, adj.cut, by = 'Sample')
df.meta.total <- df.meta.total %>% mutate(predicted.doublet.upd = ifelse(doublet_score >= Cutoff, 'True', 'False'))

gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE') %>%
  dplyr::filter(`Cellranger version` == 'v2.0')%>%
  dplyr::select(Sample, Piece_ID)
df.meta.total <- merge(df.meta.total, samples, by = 'Sample')

fwrite(df.meta.total, 'Scrublet_regATAC_15p_rate_scores_adjusted.cutoff.tsv', sep = '\t')

### now take care of combo data scrublet calls by overlapping doublet calls on atac and rna level
setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_RNAonly_withlog_mean_var/')
df.res <- list.files(path = './', pattern = 'output*', recursive = T, full.names = F)

df.meta.total <- lapply(df.res, FUN = function (x) {
  sample <- str_split_fixed(x, '_scrub', 2)[,1]
  df.metadata <- fread(x)
  df.metadata$Sample <- sample
  df.metadata$merged_barcode <- paste(sample, df.metadata$Barcodes, sep ='_')
  
  df.metadata
  return (df.metadata)
})
df.meta.total <- rbindlist(df.meta.total)

fwrite(data.table(df.meta.total$Sample %>% unique %>% sort), 'Sample_cutoffs.tsv', sep = '\t')
adj.cut <- fread('Sample_cutoffs_adj.txt', data.table = F)
df.meta.total <- merge(df.meta.total, adj.cut, by = 'Sample')
df.meta.total <- df.meta.total %>% mutate(predicted.doublet.upd = ifelse(doublet_score >= Cutoff, 'True', 'False'))
df.meta.total <- merge(df.meta.total, samples, by = 'Sample')
fwrite(df.meta.total, 'Scrublet_regATAC_scores_adjusted.cutoff.tsv', sep = '\t')


setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboATAC_withlog/')
df.res <- list.files(path = './', pattern = 'output*', recursive = T, full.names = F)

df.meta.total <- lapply(df.res, FUN = function (x) {
  sample <- str_split_fixed(x, '_scrub', 2)[,1]
  df.metadata <- fread(x)
  df.metadata$Sample <- sample
  df.metadata$merged_barcode <- paste(sample, df.metadata$Barcodes, sep ='_')
  
  df.metadata
  return (df.metadata)
})
df.meta.total <- rbindlist(df.meta.total)

fwrite(data.table(df.meta.total$Sample %>% unique %>% sort), 'Sample_cutoffs.tsv', sep = '\t')
adj.cut <- fread('Sample_cutoffs_adj.txt', data.table = F)
df.meta.total <- merge(df.meta.total, adj.cut, by = 'Sample')
df.meta.total <- df.meta.total %>% mutate(predicted.doublet.upd = ifelse(doublet_score >= Cutoff, 'True', 'False'))
df.meta.total <- merge(df.meta.total, samples, by = 'Sample')
fwrite(df.meta.total, 'Scrublet_comboATAC_scores_adjusted.cutoff.tsv', sep = '\t')

rna.scrub <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboRNA_RNAonly_withlog_mean_var/Scrublet_comboRNA_scores_adjusted.cutoff.tsv')
atac.scrub <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/comboATAC_withlog/Scrublet_comboATAC_scores_adjusted.cutoff.tsv')

total.scrub <- merge(rna.scrub, atac.scrub[,c('doublet_score', 'merged_barcode', 'Cutoff', 'predicted_doublet', 'predicted.doublet.upd')], by = 'merged_barcode', 
                     suffixes = c('.rna', '.atac'))
total.scrub
total.scrub$predicted.doublet.upd.rna.atac <- case_when(total.scrub$predicted.doublet.upd.rna & total.scrub$predicted.doublet.upd.atac ~ 'TRUE',
                                    TRUE ~ 'FALSE')
fwrite(total.scrub, '/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_comboATAC_comboRNA_scores_adjusted.cutoff.tsv', sep = '\t')

# merge regATAC scrublet and combo scrublet
reg.scrub <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/regATAC_15p_rate/Scrublet_regATAC_15p_rate_scores_adjusted.cutoff.tsv')
total.scrub <- fread('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_comboATAC_comboRNA_scores_adjusted.cutoff.tsv')
reg.scrub <- reg.scrub %>% dplyr::rename (predicted.doublet.upd.atac = predicted.doublet.upd,
                             predicted_doublet.atac = predicted_doublet,
                             doublet_score.atac = doublet_score,
                             Cutoff.atac = Cutoff)

total.scrub[,final.doublet:= predicted.doublet.upd.rna.atac]
reg.scrub[,final.doublet:=predicted.doublet.upd.atac]

reg.scrub[,`:=`(doublet_score.rna=NA, predicted_doublet.rna=NA, Cutoff.rna=NA, predicted.doublet.upd.rna=NA, predicted.doublet.upd.rna.atac=NA)]


all.scrub <- rbindlist(list(reg.scrub[,colnames(total.scrub), with = F], total.scrub))

fwrite(all.scrub, '/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_allATAC__15p_rate_adjusted.cutoff.tsv', sep= '\t')

### make HTAN BRCA combo scrublet calls
total.scrub.brca <- total.scrub %>% filter(grepl('BRCA', Sample)) %>% mutate(Piece_ID=gsub('-','_' ,Piece_ID)) %>%
  mutate(brca.barcode = paste('Multiome', Piece_ID, Barcodes, sep = '_'))
fwrite(total.scrub.brca, '/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/scrublet/HTAN_BRCA_Scrublet_comboATAC_comboRNA_scores_adjusted.cutoff.tsv', sep = '\t')




