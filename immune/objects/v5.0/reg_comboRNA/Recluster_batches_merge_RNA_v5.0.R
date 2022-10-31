# create a separate object for each batch RNA
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


runAllNormalization <- function(obj, dims=30) {
  DefaultAssay(obj) <- 'RNA'
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T,
      verbose = FALSE
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims, assay = 'SCT')
  
  obj <- NormalizeData(obj, assay = 'RNA')
  obj <- FindNeighbors(object = obj,  dims = 1:dims)
  obj <- FindClusters(object = obj,resolution = 1, verbose = FALSE)
  return(obj)
  
}

filter <- dplyr::filter
select <- dplyr::select


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
              help = "path to RNA metadata file with cell types, make cell barcodes in the 1st column",
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
cell_type_column <- opt$cell_type_column
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v5.0/Colors_panatac_v2.0.rds')


dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
int <- readRDS(input.path)
my.metadata <- fread(meta.path, data.table = F, header=TRUE) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)
int <- AddMetaData(int, my.metadata)


int$Batches <- case_when(int$Cancer %in% c('PBMC') ~ paste(int$Cancer, int$data.type, sep = '__'),
                         int$Cancer %in% c('MM') ~ int$Cancer,
                         TRUE ~ int$Chemistry)

print(table(int$Batches))

rna.split <- SplitObject(int, split.by = 'Batches')

walk(rna.split, function(obj) {
  if (!file.exists(glue::glue('PanImmune_{batch}_merged_RNA_{add_filename}.rds'))) {
    obj <- runAllNormalization(obj, dims = 40)
    batch <- obj$Batches %>% unique
    print(batch)
    saveRDS(obj, glue::glue('PanImmune_{batch}_merged_RNA_{add_filename}.rds'))
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
  }
  
  
})

################









