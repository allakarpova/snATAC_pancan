# Fubd DEGs by cell types
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
  make_option(c("--metadata"),
              type="character",
              default=NULL, 
              help="path to metadata"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
panc <- readRDS(input.path)

library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024 ^ 3)

if(!is.null(meta.path)) {
  atac.meta <- fread(meta.path) %>% data.frame(row.names = 1) %>%
    dplyr::select(all_of(cell_column))
  panc <- AddMetaData(panc, atac.meta)
  
}

print(head(panc@meta.data))
Idents(panc) <- cell_column
DefaultAssay(panc) <- 'SCT'
unique(Idents(panc))

deg <- FindAllMarkers(
  object = panc,
  only.pos = F,
  test.use = 'LR',
  latent.vars = 'nCount_SCT'
)

fwrite(deg, paste0('DEG_findAllMarkers_by_',cell_column,'_SCT_',add_filename,'.txt'), sep = '\t', row.names = T)


Idents(panc) <- cell_column
DefaultAssay(panc) <- 'RNA'
unique(Idents(panc))

deg <- FindAllMarkers(
  object = panc,
  only.pos = F,
  test.use = 'LR',
  latent.vars = 'nCount_RNA'
)

fwrite(deg, paste0('DEG_findAllMarkers_by_',cell_column,'_RNA_',add_filename,'.txt'), sep = '\t', row.names = T)










