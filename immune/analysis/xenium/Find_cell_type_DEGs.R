# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################


suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(presto))


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
  make_option(c("--metadata"),
              type="character",
              default=NULL, 
              help="path to metadata")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
panc <- readRDS(input.path)

library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024 ^ 3)

if(!is.null(meta.path)) {
  
  meta <- fread(meta.path, header = TRUE) %>%
    column_to_rownames(var = 'V1') 
  
  meta <- meta %>%    
    select(all_of(cell_column))
  
  print (head(meta))
  panc <- AddMetaData(panc, meta)
  
}

Idents(panc) <- cell_column
DefaultAssay(panc) <- 'SCT'
unique(Idents(panc))

deg <- FindAllMarkers(
  object = panc,
  only.pos = F
)

fwrite(deg, paste0('DEG_findAllMarkers_by_',cell_column,'_SCT_',add_filename,'.txt'), sep = '\t', row.names = T)

