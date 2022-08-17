# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))

################################

###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object",
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
              metavar="character"),
  make_option(c("--res"),
              type="double",
              default=NULL, 
              help="resolution used for clustering, if nothing provided then the current seurat clusters will be used")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
resol <- opt$res

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
panc <- readRDS(input.path)

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

if (!is.null(resol)) {
  panc <- FindClusters(panc, resolution = resol)
}

resol.add.name <- ifelse(is.null(resol), 'current', as.character(resol))

Idents(panc) <- 'seurat_clusters'
DefaultAssay(panc) <- 'SCT'
unique(Idents(panc))

deg <- FindAllMarkers(
  object = panc,
  only.pos = F,
  test.use = 'LR',
  latent.vars = 'nCount_SCT'
)

fwrite(deg, paste0('DEG_findAllMarkers_by_clusters_resolution_',resol.add.name,'_',add_filename,'.txt'), sep = '\t', row.names = T)












