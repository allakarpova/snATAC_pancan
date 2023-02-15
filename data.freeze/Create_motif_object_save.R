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
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output


dir.create(out_path, showWarnings = F)
setwd(out_path)


obj <- readRDS(input_path)
DefaultAssay(obj) <- 'pancan'

obj <- RegionStats(obj, assay = 'pancan', genome = BSgenome.Hsapiens.UCSC.hg38)


pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(obj),
  pwm = pfm,
  genome = 'hg38'
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

saveRDS(motif, 'Motif.object.pancan.peak.set.rds')
