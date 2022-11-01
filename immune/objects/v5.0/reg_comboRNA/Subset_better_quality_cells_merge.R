# refilter RNA cells and remove shitty cells
###libraries
##################
library(future)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))


suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(stringr))



################################

#####################################
####### FUNCTIONS ##################
####################################


runAllNormalization <- function(obj, dims) {
  #### run normalization to get initial clusters ###
  ########
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 500) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 ) %>% 
    FindNeighbors(
      reduction = 'lsi',
      dims = 2:dims ) %>% 
    FindClusters(
      algorithm = 3,
      resolution = 1,
      verbose = FALSE
    ) %>% 
    RunUMAP(dims = 2:dims,
            reduction = 'lsi')
  return(obj)
}


############################################

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
  make_option(c("--nfeature_min"),
              type="integer",
              default=200,
              help="nFeature_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--ncount_min"),
              type="integer",
              default=1000,
              help="nCount_RNA min value for filtering",
              metavar="integer")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
obj <- readRDS(input.path)

obj <- subset(obj, subset = nFeature_RNA > opt$nfeature_min &
                nCount_RNA > opt$ncount_min)

obj <- runAllNormalization(obj, dims = 40)

batch <- obj$Batches %>% unique

saveRDS(obj, glue::glue('PanImmune_{batch}_merged_ATAC_{add_filename}.rds'))
p2 <- DimPlot(obj, group.by = cell_type_column, label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_{cell_type_column}.pdf'), plot = p2, width = 12, height = 6)

p2 <- DimPlot(obj, group.by = "Cancer", cols = colors$Cancer, label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_Cancer.pdf'), plot = p2, width = 5.5, height = 4.5)

p2 <- DimPlot(obj, group.by = "Sample_type", cols = 'Paired', label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_Sample_type.pdf'), plot = p2, width = 5.5, height = 4.5)

p2 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_seurat_clusters.pdf'), plot = p2, width = 6, height = 4.5)

p2 <- DimPlot(obj, group.by = "data.type", cols = 'Paired', label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_data.type.pdf'), plot = p2, width = 5.5, height = 4.5)

p2 <- DimPlot(obj, group.by = "Piece_ID", label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_Piece_ID.pdf'), plot = p2, width = 12, height = 7)

fwrite(cbind(Embeddings(obj, reduction='umap'), obj@meta.data), glue::glue('Panimmune_{batch}_merged_{add_filename}_metadata.tsv'),
       sep='\t', row.names = T)


################









