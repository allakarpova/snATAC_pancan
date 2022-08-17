library(data.table)
library(dplyr)
library(stringr)
library(patchwork)
library(tidyverse)

dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/')
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/PDAC_HT-018_cellTyped.meta.data')
x <- x %>% mutate(Sample_ID = case_when(dataset=='HT-018-P-Slice2' ~ 'TWCE-HT-018-P-Slice2_ATAC-lib1',
                                   dataset=='HT-018-P-Slice1' ~ 'TWCE-HT-018-P-Slice1_ATAC-lib1',
                                   dataset=='HT-018-P-Slice2_V2' ~ 'TWCE-HT-018-P-Slice2_ATAC-lib1_V2'),
             Piece_ID = dataset,
             Case_ID = 'HT-018-P') %>%
  mutate(cancer.type = 'PDAC') %>%
  mutate(V1 = case_when (dataset=='HT-018-P-Slice2_V2' ~ str_split_fixed(V1, '_V2_', 2)[,2],
                          TRUE ~ str_split_fixed(V1, '_', 2)[,2]))

x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })
fwrite(x, '~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_HT-018_cellTyped.meta.data', sep = '\t')

########
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/PDAC_HT-022_cellTyped.meta.data')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID = case_when(dataset=='PDAC_HT-022-P-PunchA1' ~ 'TWCE-HT-022-P-PunchA1_ATAC-lib1',
                                        dataset=='PDAC_HT-022-P-PunchB2' ~ 'TWCE-HT-022-P-PunchB2_ATAC-lib1',
                                        dataset=='PDAC_HT-022-P-Slice1' ~ 'TWCE-HT-022-P-Slice1_ATAC-lib1',
                                        dataset=='PDAC_HT-022-P-Slice2' ~ 'TWCE-HT-022-P-Slice2_ATAC-lib1'),
                  Piece_ID = str_split_fixed(dataset, '_', 2)[,2],
                  Case_ID = 'HT-022-P') %>%
  mutate(cancer.type = 'PDAC') %>%
  mutate(V1 = str_split_fixed(V1, '_', 3)[,3])
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })
fwrite(x, '~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_HT-022_cellTyped.meta.data', sep = '\t')

#########
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/PDAC_TWCE-HT-020-P-PunchA2_ATAC-lib1_cellTyped.meta.data')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID = 'TWCE-HT-020-P-PunchA2_ATAC-lib1',
                  Piece_ID = 'HT-020-P-PunchA2',
                  Case_ID = 'HT-020-P') %>%
  mutate(cancer.type = 'PDAC')
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })

###########
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/PDAC_TWHG-CPT0092630013-XBa1_cellTyped.meta.data.txt')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID = 'TWHG-CPT0092630013-XBa1',
                  Piece_ID = 'CPT0092630013',
                  Case_ID = 'CPT0092630013') %>%
  mutate(cancer.type = 'PDAC')
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })
#########
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/PDAC_TWHG-CPT0123650013-XBa1_cellTyped.meta.data')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID = 'TWHG-CPT0123650013-XBa1',
                  Piece_ID = 'CPT0123650013',
                  Case_ID = 'CPT0123650013') %>%
  mutate(cancer.type = 'PDAC')
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })
###############
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/PDAC_TWHG-CPT0123750013-XBa1_cellTyped.meta.data')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID = 'TWHG-CPT0123750013-XBa1',
                  Piece_ID = 'CPT0123750013',
                  Case_ID = 'CPT0123750013') %>%
  mutate(cancer.type = 'PDAC')
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/PDAC/cell_typed/new_columns/PDAC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })

######## CRC
dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/CRC/cell_typed/new_columns/')
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/CRC/cell_typed/CRC_CM556C1-T1Y2_2N1Ba1_1_cellTyped.meta.data')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID = 'CM556C1-T1Y2_2N1Ba1_1',
                  Piece_ID = 'CM556C1-T1Y2',
                  Case_ID = 'CM556C1') %>%
  mutate(cancer.type = 'CRC')
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/CRC/cell_typed/new_columns/CRC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })

######
x <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/CRC/cell_typed/CRC_CM618_CellTyped.meta.data')
x$dataset %>% unique()
x <- x %>% mutate(Sample_ID =case_when(dataset=='CRC_CM618C1-T1Y2' ~ 'CM618C1-T1Y2_2N1Ba1_1',
                                       dataset=='CRC_CM618C2-S1Y2' ~ 'CM618C2-S1Y2_2N1Ba1_1'),
                  Piece_ID = str_split_fixed(Sample_ID, '_', 2)[,1],
                  Case_ID = str_split_fixed(Sample_ID, '-', 2)[,1]) %>%
  mutate(cancer.type = 'CRC') %>%
  mutate(V1 = str_split_fixed(V1, '_', 3)[,3])
x
x$Sample_ID %>% unique() %>%
  walk(function(Sample) {
    x.sub <- subset(x, Sample_ID == Sample)
    fwrite(x.sub, glue::glue('~/lab_Ding/work/single_cell/snATAC_pancan/CRC/cell_typed/new_columns/CRC_{Sample}_cellTyped.meta.data'), sep = '\t')
  })
fwrite(x, '~/lab_Ding/work/single_cell/snATAC_pancan/CRC/cell_typed/new_columns/CRC_CM618_CellTyped.meta.data', sep = '\t')
