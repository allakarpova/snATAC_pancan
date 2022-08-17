library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)


dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/clinical_links')
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/clinical_links')

meta <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/v4.0_script/immune_cells_integrated_cluster_peaks_chemistry_metadata.tsv') 
colnames(meta)
meta <- meta %>% dplyr::select (V1, Sample, Case_ID,Piece_ID,  Sample_type, data.type, Chemistry,Cancer)
ct <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/v4.0_script/mark_doublets/Celltypedata_v3.0_doublets_immune_cells_integrated_cluster_peaks_chemistry.tsv', header = T)

meta <- left_join(meta, ct, by = 'V1')
meta$cell_lin <- case_when(grepl('NK|CD|Treg', meta$Cell_type_combo_reg_doublets) ~ 'T-cells_NK',
                           grepl('DC|TAM|Micro|Mast', meta$Cell_type_combo_reg_doublets) ~ 'Myeloid',
                           grepl('B|Plasma', meta$Cell_type_combo_reg_doublets) ~ 'B-cell_Plasma',
                           TRUE ~ meta$Cell_type_combo_reg_doublets)

conditions <- c('T-cells_NK', 'Myeloid', 'B-cell_Plasma') 


ggplot(meta, aes(x=Cell_type_combo_reg_doublets, fill = Cancer)) +
  geom_bar(position = 'fill') +
  scale_fill_brewer(palette = 'Spectral') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Not_refined_Celltypedata_v3.0_doublets_immune_cells.pdf')

ggplot(meta, aes(x=Cell_type_combo_reg_doublets, fill = Cancer)) +
  geom_bar(position = 'stack') +
  scale_fill_brewer(palette = 'Spectral') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Not_refined_Celltypedata_v3.0_doublets_immune_cells_cancer_type_stack.pdf')


ggplot(subset(meta, cell_lin =='Myeloid'), aes(x=Piece_ID, fill = Cell_type_combo_reg_doublets)) +
  geom_bar(position = 'fill') +
  scale_fill_brewer(palette = 'Paired') +
  facet_grid(.~Cancer + Sample_type, scales = 'free',space = 'free') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Not_refined_Celltypedata_v3.0_doublets_Myeloid_cells_in_cancers_sample_type.pdf', width = 35, height = 6)

toplot <- meta %>% 
  filter(Cell_type_combo_reg_doublets != 'Doublet') %>%
  group_by(Cancer,Piece_ID, Cell_type_combo_reg_doublets) %>% 
  summarise(cnt = n()) %>%
  mutate(freq = round(cnt / sum(cnt), 3))

ggplot(toplot, aes(x=Cancer, y = freq, fill = Cancer)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, shape = 21, stroke = 0.1) +
  facet_wrap(~Cell_type_combo_reg_doublets, scales = 'free') +
  scale_fill_brewer(palette = 'Spectral') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Not_refined_Celltypedata_v3.0_doublets_immune_cells_boxplot_pct.pdf', width = 15, height = 10, useDingbats = F)







