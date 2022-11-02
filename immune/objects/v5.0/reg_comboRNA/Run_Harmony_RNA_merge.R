# Run Harmony over RNA object using specfic column
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


runHarmonyNormalization <- function(obj, dims=30, column = 'Piece_ID') {
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
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    FindNeighbors(reduction = "harmony", dims = 2:dims) %>%
    FindClusters(verbose = FALSE, resolution = 1) %>%
    RunUMAP(reduction = "harmony", dims = 2:dims)
  
  obj <- NormalizeData(obj, assay = 'RNA')
  
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
  make_option(c("--columnHarmony"),
              type="character",
              default="Piece_ID",
              help="Column in metadata to do Harmony over",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
ct <- opt$columnHarmony
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v5.0/Colors_panatac_v2.0.rds')

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
obj <- readRDS(input.path)


obj <- runAllNormalization(obj, dims = 40, column=ct)

batch <- obj$Batches %>% unique

saveRDS(obj, glue::glue('PanImmune_{batch}_merged_RNA_{add_filename}_harmony_{ct}.rds'))

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

p2 <- DimPlot(obj, group.by =ct, label = TRUE)
ggsave(glue::glue('Dimplot_{batch}_{add_filename}_{ct}.pdf'), plot = p2, width = 12, height = 7)

fwrite(cbind(Embeddings(obj, reduction='umap'), obj@meta.data), glue::glue('Panimmune_{batch}_merged_RNA_{add_filename}_harmony_{ct}_metadata.tsv'),
       sep='\t', row.names = T)


################



