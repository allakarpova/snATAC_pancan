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
  make_option(c("--links"),
              type="character",
              default='/diskmnt/Projects/snATAC_analysis/tumor_Alla/linkingGenes/analysis_results/tumor_only_multiome/annotate.links/Pancan_all.links.filtered_diffused_links_annotated_by_GeneHancer_scEnhancer_with_CNV_no_gained_genes_DAP_DEG.txt', 
              help="path to metadata")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output
links_path <- opt$links
#cancer <- opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)


all.links.no.gained.tumor.normal <- fread(links_path)

obj <- readRDS(input_path)
DefaultAssay(obj) <- 'pancan'
Idents(obj) <- 'Cancer'

obj <- RegionStats(obj, assay = 'pancan', genome = BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

c('BRCA','HNSCC', 'CRC', 'CESC', 'PDAC', 'OV', 'UCEC') %>% walk(function(cancer) {
  peaks.to.test <- all.links.no.gained.tumor.normal %>% filter(Cancer==cancer) %>% pull(peak) %>% unique
  Idents(obj) <- 'cell_type.harmonized.cancer.rna'
  #find peaks accessible in tumor cells
  open.peaks <- AccessiblePeaks(obj, idents = 'Tumor')
  
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(obj, assay = "pancan", slot = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[peaks.to.test, ],
    n = 50000
  )
  
  enriched.motifs <- FindMotifs(
    #assay = 'pancan',
    object = obj,
    features = peaks.to.test,
    background=peaks.matched
  )
  
  
  fwrite(enriched.motifs, glue::glue('Enriched.motifs.{cancer}.bckgr.All.cancers.open.peaks.GC.tsv'), sep='\t', row.names = F, col.names = T)
  
  
})


