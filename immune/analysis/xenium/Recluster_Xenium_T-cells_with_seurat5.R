
library(future)

#plan("multicore", workers = 20)
#options(future.globals.maxSize = 100 * 1024 ^ 3)

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))


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
  make_option(c("-m","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

renormalize <- function(obj.sub, num_pcs=30) {
  
  obj.sub <- SCTransform(obj.sub, assay = "Xenium", return.only.var.genes = T, vst.flavor="v2")
  obj.sub <- RunPCA(obj.sub, npcs = num_pcs)
  obj.sub <- RunUMAP(obj.sub, dims = 1:num_pcs, reduction.name = paste("umap.",num_pcs,"PC", sep = ""), reduction.key = paste("UMAP",num_pcs,"PC_",sep=""))
  obj.sub <- FindNeighbors(obj.sub, reduction = "pca", dims = 1:num_pcs)
  obj.sub <- FindClusters(obj.sub, resolution = 0.4)
  return(obj.sub)
}

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################
# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
  column_to_rownames('V1') 

options(Seurat.object.assay.version = "v3")

panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)

ct='T-cell'
panc.sub <- subset(x = panc.my, 
                  cells = rownames(dplyr::filter(panc.my@meta.data, 
                                                 grepl(ct, .data[[cell_column]])
                  )
                  )
)

panc.sub <- renormalize(panc.sub)

saveRDS(panc.sub,  paste0(add_filename,"_",make.names(ct), ".rds"))





