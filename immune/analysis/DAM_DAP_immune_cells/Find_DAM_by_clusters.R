# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(org.Hs.eg.db))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

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

if(!is.null(meta.path)) {
  atac.meta <- fread(meta.path) %>% data.frame(row.names = 1) %>%
    dplyr::select(seurat_clusters)
  panc <- AddMetaData(panc, atac.meta)
  
}

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

Idents(panc) <- 'seurat_clusters'
DefaultAssay(panc) <- 'chromvar'
unique(Idents(panc))

da_motif <- FindAllMarkers(
  object = panc,
  only.pos = F,
  mean.fxn=rowMeans,
  fc.name="avg_diff",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

fwrite(da_motif, paste0('DAM_findAllMarkers_by_clusters_',add_filename,'.txt'), sep = '\t', row.names = T)












