library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

library(optparse)
library(data.table)

option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to rds object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="my.data", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-a", "--assay"),
              type="character",
              default="X500peaksMACS2", 
              help="which assay should be used to merge objects? X500peaksMACS2, peaks",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay


dir.create(out_path, showWarnings = F)
setwd(out_path)

atac=readRDS(input.path)
DefaultAssay(atac)=assay.towork

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(atac),
  pwm = pfm,
  genome = 'BSgenome.Hsapiens.UCSC.hg38',
  use.counts = FALSE
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
atac <- SetAssayData(
  object = atac,
  assay = assay.towork,
  slot = 'motifs',
  new.data = motif
)

atac <- RegionStats(object = atac, genome = BSgenome.Hsapiens.UCSC.hg38)

atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(atac, paste0(add_filename,".chromvar.rds"))



