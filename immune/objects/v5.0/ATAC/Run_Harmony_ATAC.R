# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################
library(future)

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



################################

#####################################
####### FUNCTIONS ##################
####################################


runHarmonyNormalization <- function(obj, dims) {
  obj$Data.source <- ifelse(obj$Cancer == 'PBMC', '10x', 'DingLab')
  if(conditions=='B-cell_Plasma') {
    obj$Batches <- case_when(obj$Cancer %in% c('PBMC') ~ 'PBMC',
                                   obj$Cancer %in% c('MM') ~ obj$Cancer,
                                   TRUE ~ obj$Chemistry)
  }  else {
    obj$Batches <- case_when(obj$Cancer %in% c('PBMC') ~ paste(obj$Cancer, obj$data.type, sep = '__'),
                             obj$Cancer %in% c('MM') ~ obj$Cancer,
                                   TRUE ~ obj$Chemistry)
  }
  #### run normalization to get initial clusters ###
  ########
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 500) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 ) %>% 
    RunHarmony( group.by.vars = c("Batches"), 
                reduction = 'lsi', 
                project.dim = F, 
                max.iter.harmony = 20)
    FindNeighbors(
      reduction = 'harmony',
      dims = 2:dims ) %>% 
    FindClusters(
      algorithm = 4,
      resolution = 1,
      verbose = FALSE
    ) %>% 
    RunUMAP(dims = 2:dims,
            reduction = 'harmony')
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
  make_option(c("-t","--metadata.file"),
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

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v5.0/Colors_panatac_v2.0.rds')


dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
obj <- readRDS(input.path)
DefaultAssay(obj) <- 'ATAC_immune'

obj <- runHarmonyNormalization(obj, dims = 40)

ct='Batches'
saveRDS(obj, glue::glue('PanImmune_merged_ATAC_{add_filename}_harmony_{ct}.rds'))

p2 <- DimPlot(obj, group.by = "Cancer", cols = colors$Cancer, label = TRUE)
ggsave(glue::glue('Dimplot_{add_filename}_Cancer.pdf'), plot = p2, width = 5.5, height = 4.5)

p2 <- DimPlot(obj, group.by = cell_column, label = TRUE)
ggsave(glue::glue('Dimplot_{add_filename}_{cell_column}.pdf'), plot = p2, width = 12, height = 7)

p2 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE)
ggsave(glue::glue('Dimplot_{add_filename}_seurat_clusters.pdf'), plot = p2, width = 6, height = 4.5)

p2 <- DimPlot(obj, group.by = "data.type", cols = 'Paired', label = TRUE)
ggsave(glue::glue('Dimplot_{add_filename}_data.type.pdf'), plot = p2, width = 5.5, height = 4.5)

p2 <- DimPlot(obj, group.by = "Piece_ID", label = TRUE)
ggsave(glue::glue('Dimplot_{add_filename}_Piece_ID.pdf'), plot = p2, width = 12, height = 7)

p2 <- DimPlot(obj, group.by =ct, label = TRUE)
ggsave(glue::glue('Dimplot_{add_filename}_{ct}.pdf'), plot = p2, width = 12, height = 7)

fwrite(cbind(Embeddings(obj, reduction='umap'), obj@meta.data), glue::glue('Panimmune_merged_ATAC_{add_filename}_harmony_{ct}_metadata.tsv'),
       sep='\t', row.names = T)










