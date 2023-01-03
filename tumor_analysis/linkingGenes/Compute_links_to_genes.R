# Alla Karpova create pancan combo object with RNA and ATAC

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))


###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="input object",
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

dir.create('out')


cat('opening object \n')
obj <- readRDS(input_path)
DefaultAssay(obj) <- "pancan"
cat('done \n')

# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "pancan",
  distance = 5e+05,
  n_sample = 1000,
  expression.assay = "SCT"
)

out.obj <- str_replace(input_path, pattern = 'rds', replacement = 'linked.rds')
out.obj <- str_split(out.obj, '[/]')[[1]][str_count(out.obj, '/')+1]

print(out.obj)
saveRDS(obj, paste0('out/', out.obj))
