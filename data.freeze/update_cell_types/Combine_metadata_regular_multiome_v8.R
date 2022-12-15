
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))


option_list = list(
  make_option(c("-c", "--cancer"),
              type="character",
              default="BRCA", 
              help="cancer type",
              metavar="character"),
  make_option(c("-i", "--input.folder"),
              type="character",
              default="/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v3.1", 
              help="output folder path",
              metavar="character"),
  make_option(c("-v", "--version"),
              type="character",
              default="v2.0", 
              help="data freeze version",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cancer.type <- opt$cancer
input.folder <- opt$input.folder
v <- opt$version

# this script will take metadata tables from regular ATAC and combo objects and create a single metadats file
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
#samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
#samples$Cancer_piece = paste(samples$`Disease Type`, samples$Piece_ID, sep = '_')
samples$`Disease Type`[samples$`Disease Type`=='PKD'] <- 'ccRCC'
samples <- samples %>% dplyr::filter(`Disease Type` == cancer.type)

samples.id <- samples$Sample %>% as.character()
length(samples.id)


setwd(input.folder)
meta.files <- list.files(path = input.folder, pattern = '*data', full.names = T)
if (cancer.type=='ccRCC') {
  meta.files <- meta.files[grepl('ccRCC|PKD', meta.files)]
} else {
  meta.files <- meta.files[grepl(cancer.type, meta.files)]
}
print(meta.files)

regular.sample <- samples %>% subset (`Data Type` == 'snATAC' | `Data Type` == 'scATAC') %>% pull (Sample)
combo.sample <- samples %>% subset (`Data Type` == '10x_SC_Multi_ATAC_SEQ') %>% pull (Sample)
print('combo samples:')
print(combo.sample)
print('regular samples:')
print(regular.sample)

#do regular ATAC metadata first
if (length(regular.sample) > 0) {
  cat('making regular meta\n')
  regular.metadata <- lapply(regular.sample, function(x) {
    caseid <- samples$`Case ID`[samples$Sample==x]
    pieceid <- samples$Piece_ID[samples$Sample==x]
    type <- samples$`Sample Type`[samples$Sample==x]
    datatype <- samples$`Data Type`[samples$Sample==x]
    chemistry <- samples$`Chemistry`[samples$Sample==x]
    toret <- fread(meta.files[grepl(x,meta.files)], data.table = F) %>%
      select(-starts_with('prediction.score')) %>% # this removes columns with prediction scores for different cell types because each file can have different # of those
      mutate(Sample = x,
             Case_ID = caseid,
             Piece_ID = pieceid,
             Sample_type = type,
             data.type = datatype,
             Chemistry = chemistry,
             Cancer = cancer.type)
    return(toret)
  })
  regular.metadata <- do.call('rbind', regular.metadata)
  rownames(regular.metadata) <- paste(regular.metadata$Cancer,regular.metadata$Piece_ID, regular.metadata$V1, sep = '_')
} else {
  regular.metadata <- NULL
}

load.combo <- function(x) {
  caseid <- samples$`Case ID`[samples$Sample==x]
  pieceid <- samples$Piece_ID[samples$Sample==x]
  type <- samples$`Sample Type`[samples$Sample==x]
  datatype <- samples$`Data Type`[samples$Sample==x]
  chemistry <- samples$`Chemistry`[samples$Sample==x]
  toret <- fread(meta.files[grepl(x,meta.files)], data.table = F) %>%
    select(-starts_with('gex_')) %>% # this removes all gex data quality columns
    mutate(Sample = x,
           Case_ID = caseid,
           Piece_ID = pieceid,
           Sample_type = type,
           data.type = datatype,
           Chemistry = chemistry,
           Cancer = cancer.type)
  return(toret)
}
#then if there is combo metadata
if (length(combo.sample) > 0) {
  cat('making combo meta\n')
  if (length(combo.sample) == 1) {
    combo.metadata <- load.combo(combo.sample)
    if ('pct_read_in_peaks_500MACS2' %in% colnames(combo.metadata)) {combo.metadata <- combo.metadata %>% dplyr::rename(frip_500MACS2=pct_read_in_peaks_500MACS2)}
    if ('peak_region_fragments_500MACS2' %in% colnames(combo.metadata)) {combo.metadata <- combo.metadata %>% dplyr::rename(peak_RF_500MACS2=peak_region_fragments_500MACS2)}
  } else {
    combo.metadata <- lapply(combo.sample, load.combo)
    combo.metadata <- do.call('rbind', combo.metadata)
  }
  
  rownames(combo.metadata) <- paste(combo.metadata$Cancer,combo.metadata$Piece_ID, combo.metadata$V1, sep = '_')
  
  # make colnames uniform
  combo.metadata <- combo.metadata %>% dplyr::rename(is__cell_barcode = is_cell,
                                                     chimeric = atac_chimeric_reads,
                                                     unmapped = atac_unmapped_reads,
                                                     lowmapq = atac_lowmapq,
                                                     duplicate = atac_dup_reads,
                                                     total = atac_raw_reads,
                                                     mitochondrial = atac_mitochondrial_reads,
                                                     TSS_fragments = atac_TSS_fragments,
                                                     passed_filters = atac_fragments,
                                                     peak_region_fragments = atac_peak_region_fragments,
                                                     peak_region_cutsites = atac_peak_region_cutsites,
                                                     #peak_RF_500MACS2 = peak_region_fragments_500MACS2,
                                                     #pct_reads_in_peaks = pct_read_in_peaks_500MACS2
  )
  
  rownames(combo.metadata) <- paste(combo.metadata$Cancer,combo.metadata$Piece_ID, combo.metadata$V1, sep = '_')
  
  if (length(regular.sample) > 0) {
    cat('rbind combo and regular\n')
    common.columns <- intersect(colnames(regular.metadata), colnames(combo.metadata))
    total.metadata <- rbind(regular.metadata[common.columns], 
                            combo.metadata[common.columns])
  } else {
    total.metadata <- combo.metadata
  }
  
} else {
  total.metadata <- regular.metadata
}

#print(head(total.metadata))
# #now harmonize cell types
cat('harmonize cell types\n')
#print(head(total.metadata))
total.metadata$cell_type %>% unique

if(cancer.type == 'BRCA') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial|LSEC', cell_type) ~ 'Endothelial',
                                                                                grepl('CAF|Fibroblast|Stell', cell_type) ~ 'Fibroblasts',
                                                                                grepl('CD4_T|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                cell_type =='B' | grepl('B-cells', cell_type)  ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('Macro', cell_type) ~ 'Macrophages',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('Basal prog', cell_type) ~ 'Basal progenitors',
                                                                                grepl('Luminal prog', cell_type) ~ 'Luminal progenitors',
                                                                                grepl('Luminal mature', cell_type) ~ 'Luminal mature',
                                                                                grepl('Breast|Basal', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('Adipo', cell_type) ~ 'Adipocytes',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Tumor', cell_type) ~ 'Tumor',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('Hepatocyte', cell_type) ~ 'Hepatocytes'))
} else if (cancer.type == 'CRC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial|LSEC', cell_type) ~ 'Endothelial',
                                                                                grepl('tellate|ibroblast', cell_type) ~ 'Fibroblasts',
                                                                                cell_type =='T' | grepl('CD4_T|CD8|T-cells|2', cell_type) ~ 'T-cells',
                                                                                cell_type =='B' | grepl('B-cells', cell_type) ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('Macro|Kupffer', cell_type) ~ 'Macrophages',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('Hepatocyte', cell_type) ~ 'Hepatocytes',
                                                                                grepl('Cholangiocytes', cell_type) ~ 'Cholangiocytes',
                                                                                grepl('Goblet', cell_type) ~ 'Goblet',
                                                                                grepl('Alveolar', cell_type) ~ 'Alveolar',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Tumor', cell_type) ~ 'Tumor',
                                                                                grepl('Unknown', cell_type) ~ 'Unknown',
                                                                                grepl('Immune', cell_type) ~ 'Immune cells',
                                                                                grepl('Lymph', cell_type) ~ 'Lymphocytes',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('Oligodendrocyte', cell_type) ~ 'Oligodendrocytes',
                                                                                grepl('Neuron', cell_type) ~ 'Neurons',
                                                                                grepl('Enterocyte', cell_type) ~ 'Enterocytes',
                                                                                grepl('Epithelial', cell_type) ~ 'Epithelial',
                                                                                TRUE ~ cell_type))
} else if (cancer.type == 'ccRCC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial|Endo', cell_type) ~ 'Endothelial',          
                                                                                grepl('ibroblast', cell_type) ~ 'Fibroblasts',
                                                                                grepl('CD4_T|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('CD|DCT|LOH|Podo|PT|Mesa', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('B-cells', cell_type) ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('Macro', cell_type) ~ 'Macrophages',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('umor', cell_type) ~ 'Tumor',
                                                                                grepl('Unknown', cell_type) ~ 'Unknown',
                                                                                grepl('Immune', cell_type) ~ 'Immune cells',
                                                                                TRUE ~ cell_type))
} else if (cancer.type == 'GBM') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                cell_type =='T cell' | grepl('CD4_T|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('B-cell', cell_type) ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('TAM', cell_type) ~ 'Macrophages',
                                                                                grepl('umor', cell_type) ~ 'Tumor',
                                                                                grepl('Unknown', cell_type) ~ 'Unknown',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('Oligodendrocyte', cell_type) ~ 'Oligodendrocytes',
                                                                                grepl('Neuron', cell_type) ~ 'Neurons',
                                                                                TRUE ~ cell_type))
} else if (cancer.type == 'CESC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('B-cell', cell_type) ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('ibroblast', cell_type) ~ 'Fibroblasts',
                                                                                grepl('Macrop|Mono', cell_type) ~ 'Macrophages',
                                                                                grepl('Malignant|Tumor', cell_type) ~ 'Tumor',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('Epithelial', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('Erythrocyte', cell_type) ~ 'Erythrocytes',
                                                                                TRUE ~ cell_type))
} else if (cancer.type == 'HNSCC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                cell_type =='B' | grepl('B-cell', cell_type) ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('ibroblast', cell_type) ~ 'Fibroblasts',
                                                                                grepl('Macrop|Mono', cell_type) ~ 'Macrophages',
                                                                                grepl('Malignant|Tumor', cell_type) ~ 'Tumor',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('Epithelial', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('Erythrocyte', cell_type) ~ 'Erythrocytes',
                                                                                TRUE ~ cell_type))
} else if (cancer.type == 'MM') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                cell_type =='B' | grepl('preB', cell_type) ~ 'B-cells',
                                                                                cell_type =='Plasma' & Sample_type == 'Normal' ~ 'Plasma',
                                                                                cell_type =='Plasma' & Sample_type != 'Normal' ~ 'Tumor',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('ibroblast', cell_type) ~ 'Fibroblasts',
                                                                                grepl('Macrop|Mono', cell_type) ~ 'Macrophages',
                                                                                grepl('Neutrophil', cell_type) ~ 'Neutrophils',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('Epithelial', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('Erythrocyte', cell_type) ~ 'Erythrocytes',
                                                                                TRUE ~ cell_type))
} else if (cancer.type=='OV') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                grepl('B-cell', cell_type) ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('umor', cell_type) ~ 'Tumor',
                                                                                grepl('ibroblast', cell_type) ~ 'Fibroblasts',
                                                                                grepl('Macrop|Mono', cell_type) ~ 'Macrophages',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Unknown', cell_type) ~ 'Unknown',
                                                                                TRUE ~ cell_type))
} else if (cancer.type =='PDAC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial|LSEC', cell_type) ~ 'Endothelial',
                                                                                grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                grepl('Lymph', cell_type) ~ 'Lymphocytes',
                                                                                cell_type =='B-cells' | cell_type =='B' ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('umor', cell_type) ~ 'Tumor',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('Hepatocyte', cell_type) ~ 'Hepatocytes',
                                                                                grepl('Cholangiocytes', cell_type) ~ 'Cholangiocytes',
                                                                                grepl('CAF|Fibroblast|tellate', cell_type) ~ 'Fibroblasts',
                                      
                                                                                grepl('Macrop|Mono|Myeloid|Kupffer', cell_type) ~ 'Macrophages',
                                                                                grepl('Islet', cell_type) ~ 'Islets',
                                                                                grepl('Acinar|ADM', cell_type) ~ 'Acinar',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Normal_duct|Ductal', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('Immune', cell_type) ~ 'Immune cells',
                                                                                TRUE ~ cell_type))
} else if (cancer.type =='PKD') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endo', cell_type) ~ 'Endothelial',
                                                                                cell_type =='T' ~ 'T-cells',
                                                                                cell_type =='B-cells' | cell_type =='B' ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('cyst', cell_type) ~ 'Cyst',
                                                                                grepl('Macro|Mono', cell_type) ~ 'Macrophages',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('unknown', cell_type) ~ 'Unknown',
                                                                                grepl('CD|DCT|LOH|Podo|PT|Mesa', cell_type) ~ 'Normal epithelial cells',
                                                                                TRUE ~ cell_type))
} else if (cancer.type =='UCEC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                grepl('Epithelial', cell_type) ~ 'Normal epithelial cells',
                                                                                cell_type =='T' | grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                cell_type =='B-cells'  ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('umor', cell_type) ~ 'Tumor',
                                                                                grepl('Fibroblast', cell_type) ~ 'Fibroblasts',
                                                                                grepl('Macrop|Mono|Myeloid', cell_type) ~ 'Macrophages',
                                                                                grepl('Immune', cell_type) ~ 'Immune cells',
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                TRUE ~ cell_type))
} else if (cancer.type =='PBMC') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (cell_type =='T' | grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                grepl('B ', cell_type) ~ 'B-cells',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('Mono', cell_type) ~ 'Monocytes',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                TRUE ~ cell_type))
} else if (cancer.type =='SKCM') {
  total.metadata <- total.metadata %>% mutate(cell_type.harmonized = case_when (grepl('Endothelial', cell_type) ~ 'Endothelial',
                                                                                grepl('CD4|CD8|T-cells', cell_type) ~ 'T-cells',
                                                                                grepl('Treg', cell_type) ~ 'Tregs',
                                                                                
                                                                                cell_type =='B-cells' | cell_type =='B' ~ 'B-cells',
                                                                                grepl('Plasma', cell_type) ~ 'Plasma',
                                                                                grepl('DC', cell_type) ~ 'DC',
                                                                                grepl('umor', cell_type) ~ 'Tumor',
                                                                                grepl('Pericyte', cell_type) ~ 'Pericytes',
                                                                                grepl('Fibroblast|tellate', cell_type) ~ 'Fibroblasts',
                                                                                grepl('Keratino', cell_type) ~ 'Keratinocytes',
                                                                                
                                                                                grepl('Macrop|Mono|Myeloid|Kupffer', cell_type) ~ 'Macrophages',
                                          
                                                                                grepl('Mast', cell_type) ~ 'Mast',
                                                                                grepl('NK', cell_type) ~ 'NK',
                                                                                grepl('Normal_duct', cell_type) ~ 'Normal epithelial cells',
                                                                                grepl('Immune', cell_type) ~ 'Immune cells',
                                                                                TRUE ~ cell_type))
print(paste('Cells not cell typed:', sum(is.na(total.metadata$cell_type.harmonized))))

#now label transfer did not work well for some samples in regular ATAC, I updated some cell types
fwrite(total.metadata, paste0(cancer.type,'_',length(samples.id),'_samples_metadata_v8_data_freeze_',v,'.tsv'), sep = '\t', row.names = T)


