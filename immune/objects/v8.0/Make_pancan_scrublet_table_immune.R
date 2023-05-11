#Alla Karpova
## make pancan doublet table
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))



###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = "Scrublet location", trim_ws = T)
samples

samples$Keep <- samples$`Include in immune` %>% unlist()
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::select(`Disease Type`, Sample,Piece_ID, `Data Type`, `Scrublet doublet Table`)

data.path <- samples$`Scrublet doublet Table`
names(data.path) <- samples$Sample

cancer.piece <- paste(samples$`Disease Type`, samples$Piece_ID, sep = '_')
names(cancer.piece) <- samples$Sample

regular.sample <- samples %>% subset (`Data Type` == 'snATAC' | `Data Type` == 'scATAC') %>% pull (Sample)
combo.sample <- samples %>% subset (`Data Type` == '10x_SC_Multi_ATAC_SEQ') %>% pull (Sample)
regular.sample.rna <- samples %>% subset (`Data Type` == 'snRNA' | `Data Type` =='10x_SC_Multi_GEX') %>% pull (Sample)

print('combo samples:')
print(combo.sample)
print('regular samples:')
print(regular.sample)


regular.scrublet <- regular.sample %>% map(function(x){
  #print(x)
  toret <- fread(data.path[x], header = TRUE, data.table = F) 
  #print((toret))
  toret$Barcodes <- paste(cancer.piece[x], toret$Barcodes, sep ='_')
  toret <- toret %>% 
    
    rename(doublet_score_atac=doublet_score,
                   predicted_doublet_atac=predicted_doublet) %>%
    mutate(doublet_score_rna = NA,
           predicted_doublet_rna = NA,
           Doublet_final = predicted_doublet_atac)
  return(toret)
}) %>% rbindlist() %>% data.frame()
head(regular.scrublet)

regular.scrublet.rna <- regular.sample.rna %>% map(function(x){
  print(x)
  toret <- fread(data.path[x], header = TRUE, data.table = F) 
  #print((toret))
  toret$Barcodes <- paste(cancer.piece[x], toret$Barcodes, sep ='_')
  toret <- toret %>% 
    
    rename(doublet_score_rna=doublet_score,
           predicted_doublet_rna=predicted_doublet) %>%
    mutate(doublet_score_atac = NA,
           predicted_doublet_atac = NA,
           Doublet_final = predicted_doublet_rna)
  return(toret)
}) %>% rbindlist() %>% data.frame()
head(regular.scrublet.rna)

combo.scrublet <- combo.sample %>% map(function(x){
  #print(x)
  toret <- fread(data.path[x], header = TRUE, data.table = F) 
  #print((toret))
  toret$Barcodes <- paste(cancer.piece[x], toret$Barcodes, sep ='_')
  toret <- toret %>% 
    mutate(Doublet_final = case_when(predicted_doublet_atac=='True' &  predicted_doublet_rna=='True' ~ 'True',
                                     TRUE ~ 'False'))
  return(toret)
}) %>% rbindlist()  %>% data.frame()

total.scrublet.atac <- rbind(regular.scrublet[colnames(combo.scrublet)], combo.scrublet)

fwrite(total.scrublet.atac, 
       paste0('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/snATAC/scrublet/All_',(length(regular.sample) + length(combo.sample)), '_samples_metadata_immune_data_freeze_v8.0.tsv'), 
       sep ='\t', row.names = FALSE)

total.scrublet.rna <- rbind(regular.scrublet.rna[colnames(combo.scrublet)], combo.scrublet)

fwrite(total.scrublet.rna, 
       paste0('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/snATAC/scrublet/All_',(length(regular.sample.rna) + length(combo.sample)), '_samples_metadata_immune_data_freeze_v8.0.tsv'), 
       sep ='\t', row.names = FALSE)


