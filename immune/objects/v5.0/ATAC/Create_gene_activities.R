#Create gene activities on all protein coding genes and save
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))



###options###
######################
option_list = list(
  make_option(c("-a", "--input.atac.object"),
              type="character",
              default=NULL, 
              help="path to integrated ATAC object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character")
 
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

# read in initial arguments

input.path.atac <- opt$input.atac.object
out_path <- opt$output
add_filename <- opt$extra


dir.create(out_path, showWarnings = F)
setwd(out_path)


panc.atac <- readRDS(input.path.atac)

plan("multicore", workers = 20)
options(future.globals.maxSize = 250 * 1024^3) # for 250 Gb RAM

gene.activities <- GeneActivity(panc.atac)
  
saveRDS(gene.activities, paste0('ATAC_GENEACTIVITY_', add_filename, '.rds'))






