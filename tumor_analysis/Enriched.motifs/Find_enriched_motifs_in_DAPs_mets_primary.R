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
    n = length(peaks.to.test) + 500
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
  print(length(peaks.to.test))
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(obj, assay = "pancan", slot = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[peaks.to.test, ],
    n = length(peaks.to.test) + 500
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
              metavar="character")
  
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


obj@meta.data <- obj@meta.data %>% mutate(Cancer_stage_cell_type=case_when(cell_type.harmonized.cancer=='Tumor' ~ paste(Cancer,Sample_type, 'Tumor', sep='__'),
                                                                       TRUE ~ cell_type.harmonized.cancer))


Idents(obj) <- 'cell_type.harmonized.cancer'
table(Idents(obj))

open.peaks <- AccessiblePeaks(obj, idents = 'Tumor')
tryCatch( {
  enriched.m <- findEnrichedMotifs(cancer=cancer.type, open.peaks)
  fwrite(enriched.m, glue::glue('Motifs_enriched_in_{cancer.type}_mets_vs_primary.tsv'), row.names = F, sep='\t')
},
error = function(e){
  message(e)
})
depleted.m <- findDepletedMotifs(cancer=cancer.type, open.peaks)
fwrite(depleted.m, glue::glue('Motifs_depleted_in_{cancer.type}_mets_vs_primary.tsv'), row.names = F, sep='\t')




