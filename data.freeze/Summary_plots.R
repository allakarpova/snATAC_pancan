library(googlesheets4)
library(stringr)
library(ggplot2)
library(dplyr)
library(cowplot)

gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
#samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::select(`Disease Type`, Sample,`Sample Type`, `Data Type`, `Data folder`, Chemistry)
samples$`Disease Type`[samples$`Disease Type`=='PKD'] <- 'ccRCC'
samples

dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v4.0_data_freeze/')
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v4.0_data_freeze/')
ggplot(data = samples, aes (x = `Disease Type`)) +
  geom_bar(aes(fill =  `Data Type`)) +
  theme_cowplot() +
  scale_fill_brewer(palette = 'Paired') +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, colour = "black")
ggsave('Summary_plot_data.type.pdf', width = 10, height = 5, useDingbats = F)

ggplot(data = samples, aes (x = `Disease Type`)) +
  geom_bar(aes(fill =  `Sample Type`)) +
  theme_cowplot() +
  scale_fill_brewer(palette = 'Set2') +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, colour = "black")
ggsave('Summary_plot_sample_type.pdf', width = 10, height = 5, useDingbats = F)


library(Signac)
library(Seurat)
library(data.table)
library(RColorBrewer)

add_filename <- 'HNSCC_5_snATAC_Merged_06283021'
panc <- readRDS(paste0(add_filename,'.rds'))
meta <- fread('HNSCC_5_samples_metadata_data_freeze_v2.0.tsv', data.table = F) %>% data.frame(row.names = 1, check.rows = F, check.names = F)
panc <- AddMetaData(panc, meta)

n <- length(unique(panc$cell_type.harmonized))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p1 <- DimPlot(panc, group.by = 'data.type', pt.size = 0.1, cols = 'Paired') + 
  ggplot2::ggtitle("Combined snATAC data.type")

p2 <- DimPlot(panc, pt.size = 0.1,label=T, group.by = 'cell_type.harmonized') + 
  ggplot2::ggtitle("Combined snATAC cell_type.harmonized")
pdf(paste0(add_filename,'_celltype.pdf'),height=6,width=16, useDingbats = F)
p1+p2
dev.off()

p <- DimPlot(panc, group.by = 'Sample_type', pt.size = 0.1, cols = 'Dark2') + 
  ggplot2::ggtitle("Combined snATAC Sample_type")
pdf(paste0(add_filename,'_tumor_normal.pdf'),height=6,width=8, useDingbats = F)
p
dev.off()


# plot QC plots
meta <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/data_freeze/v4.0_data_freeze/All_149_samples_metadata_data_freeze_v4.1.tsv')
head(meta)
ggplot(meta, aes_string(y = 'TSS.enrichment', x = 'Cancer', fill = 'Cancer')) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA,# inherit.aes = FALSE,
              #mapping = aes_string(x = 'ident', y = 'TSS.enrichment', fill = 'ident'),
               alpha = 0.5,
               color = "black") +
  scale_fill_brewer(palette = 'Spectral') +
  theme_cowplot()
ggsave('TSS.enrochment.plot.pdf', width = 10, height = 5)


ggplot(meta, aes(y = log10(passed_filters), x = Cancer, fill = Cancer)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA,# inherit.aes = FALSE,
               #mapping = aes_string(x = 'ident', y = 'TSS.enrichment', fill = 'ident'),
               alpha = 0.5,
               color = "black") +
  scale_fill_brewer(palette = 'Spectral') +
  theme_cowplot()
ggsave('log10_passed_filters.plot.pdf', width = 10, height = 5)



ggplot(meta, aes_string(y = 'TSS.enrichment', x = 'Sample_type', fill = 'Sample_type')) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA,# inherit.aes = FALSE,
               #mapping = aes_string(x = 'ident', y = 'TSS.enrichment', fill = 'ident'),
               alpha = 0.5,
               color = "black") +
  facet_wrap(~Cancer, scales = 'free_x') +
  scale_fill_brewer(palette = 'Spectral') +
  theme_cowplot()
ggsave('TSS.enrochment.plot.pdf', width = 10, height = 5)

ggplot((meta %>% filter(cell_type.harmonized.cancer %in% c('T-cells', 'Tumor', 'Macrophages', 'B-cells', 'Plasma', 'Fibroblasts', 'Endothelial'))), 
       aes(y = log10(passed_filters), x = cell_type.harmonized.cancer, fill = cell_type.harmonized.cancer)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA,# inherit.aes = FALSE,
               #mapping = aes_string(x = 'ident', y = 'TSS.enrichment', fill = 'ident'),
               alpha = 0.5,
               color = "black") +
  #facet_wrap(~Cancer, scales = 'free_x') +
  scale_fill_brewer(palette = 'Paired') +
  theme_cowplot()
ggsave('log10_passed_filters.plot_by_cell_type.pdf', width = 10, height = 5)





