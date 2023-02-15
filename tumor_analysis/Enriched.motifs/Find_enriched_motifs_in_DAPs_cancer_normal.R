suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))
suppressMessages(library(JASPAR2020))
suppressMessages(library(TFBSTools))
suppressMessages(library(ChIPseeker))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
library("BSgenome.Hsapiens.UCSC.hg38")
library(motifmatchr)

plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))

### FUNCTIONS ###

findEnrichedMotifs <- function(cancer, open.peaks) {
  peaks.to.test <- dap.relaxed %>% filter(p_val_adj < 0.05 & Disease==cancer & avg_log2FC > 0) %>% pull(peak)
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(obj, assay = "pancan", slot = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[peaks.to.test, ],
    n = 5000
  )
  
  enriched.motifs <- FindMotifs(
    #assay = 'pancan',
    object = obj,
    features = peaks.to.test,
    background=peaks.matched
  )
  return(enriched.motifs)
}

findDepletedMotifs <- function(cancer, open.peaks) {
  peaks.to.test <- dap.relaxed %>% filter(p_val_adj < 0.05 & Disease==cancer & avg_log2FC < 0) %>% pull(peak)
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(obj, assay = "pancan", slot = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[peaks.to.test, ],
    n = 5000
  )
  
  enriched.motifs <- FindMotifs(
    #assay = 'pancan',
    object = obj,
    features = peaks.to.test,
    background=peaks.matched
  )
  return(enriched.motifs)
}
###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="Input multiome object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--dap"),
              type="character",
              default=NULL, 
              help="path to DAP table",
              metavar="character"),
  make_option(c("--cancer"),
              type="character",
              default=NULL, 
              help="cancer type",
              metavar="character"),
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output
dap.path <- opt$dap
cancer.type <- opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)

obj <- readRDS(input_path)
DefaultAssay(obj) <- 'pancan'

obj <- RegionStats(obj, assay = 'pancan', genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(obj) <- 'pancan'
motif <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/Motif.object/Motif.object.pancan.peak.set.rds')
Motifs(obj) <- motif

dap.relaxed <- fread(dap.path, header = TRUE)

atac.meta <- fread('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/All_225_samples_metadata_data_freeze_v7.0.tsv') %>%
  data.frame(row.names = 1)
obj <- AddMetaData(obj, atac.meta)

if (cancer.type %in% c('ccRCC', 'PDAC', 'CRC')) {
  if (cancer.type %in% c('ccRCC')) {
    normal.meta <- fread(glue::glue('/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/v7.0/{cancer.type}/ATAC/cell_typing/{cancer.type}_ATAC_normal_cell_type_regular.tsv'), header = T) %>%
      data.frame(row.names = 1)
    
  } else {
    normal.meta <- fread(glue::glue('/diskmnt/Projects/snATAC_analysis/tumor_Alla/Normal_cell_typing_detailed/v7.0/{cancer.type}/ATAC/cell_typing/{cancer.type}_ATAC_normal_cell_type_combo_regular.tsv'), header = T) %>%
      data.frame(row.names = 1)
    
  }
  
  obj <- AddMetaData(obj, normal.meta)
  obj$cell_type.normal <- case_when(is.na(obj$cell_type.normal) ~ obj$cell_type.harmonized.cancer,
                                    TRUE ~ obj$cell_type.normal)
  
  table(obj$cell_type.normal)
  Idents(obj) <- 'cell_type.normal'
} else {
  Idents(obj) <- 'cell_type.harmonized.cancer'
}

if (cancer.type=='BRCA') {
  obj$Cancer <- case_when(obj$Piece_ID %in% c("HT268B1-Th1H3", "HT029B1-S1PC", "HT035B1-S1PA",
                                              "HT1408-06","HT141B1-S1H1", "HT206B1-S1H4", "HT271B1-S1H3",
                                              "HT378B1-S1H1", "HT378B1-S1H2", "HT384B1-S1H1", "HT517B1-S1H1") ~ 'BRCA_Basal',
                          TRUE ~ 'BRCA')
  obj$Cancer_cell_type <- case_when(obj$cell_type.harmonized.cancer=='Tumor', paste(Cancer, cell_type.harmonized.cancer, sep='__'),
                                    TRUE ~ obj$cell_type.harmonized.cancer)
  Idents(obj) <- 'Cancer_cell_type'
}

table(Idents(obj))

cancer.normal.pairs <- data.frame(Cancer = c('PDAC', 'ccRCC', "CRC", 'UCEC', 'CESC', 'GBM', 'MM', 'OV', 'HNSCC', 'SKCM'),
                                  Normal = c('Ductal-like2', 'Proximal Tubule', 'Distal Stem Cells', 'Secretory Endometrial epithelial cells', 
                                             'Normal squamous cells', 'OPC', 'B-cells', 'Secretory Endometrial epithelial cells', 'Normal squamous cells', 'Melanocytes'))
#find peaks accessible in tumor cells and prximal tubules
if (cancer.type=='BRCA') {
   walk(c('BRCA', 'BRCA_Basal'), c('Luminal mature', 'Luminal progenitor'), function(c, n) {
    open.peaks <- AccessiblePeaks(obj, idents = c(glue::glue("{c}__Tumor"), n))
    enriched.m <- findEnrichedMotifs(cancer=c, open.peaks)
    fwrite(enriched.m, glue::glue('Motifs_enriched_in_{c}_cancer_vs_{n}.tsv'), row.names = F, sep='\t')
    depleted.m <- findDepletedMotifs(cancer=c, open.peaks)
    fwrite(depleted.m, glue::glue('Motifs_depleted_in_{c}_cancer_vs_{n}.tsv'), row.names = F, sep='\t')
    
  })
  
} else {
    normal <- cancer.normal.pairs %>% filter(Cancer==cancer.type) %>% pull(Normal)
    open.peaks <- AccessiblePeaks(obj, idents = c("Tumor", normal))
    enriched.m <- findEnrichedMotifs(cancer=cancer.type, open.peaks)
    fwrite(enriched.m, glue::glue('Motifs_enriched_in_{cancer.type}_cancer_vs_{normal}.tsv'), row.names = F, sep='\t')
    depleted.m <- findDepletedMotifs(cancer=cancer.type, open.peaks)
    fwrite(depleted.m, glue::glue('Motifs_depleted_in_{cancer.type}_cancer_vs_{normal}.tsv'), row.names = F, sep='\t')
    
}



