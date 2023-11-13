#Run harmony between multiome object and any other RNA object to figure out cell types on the external object


suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
library(harmony)

normalize_rna_harmony <- function(obj, dims=30, column = 'Piece_ID_RNA') {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      vars.to.regress =  c("percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T,
      verbose = FALSE)
  
  obj <- obj %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    RunUMAP(reduction = "harmony", dims = 1:dims, reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  
  obj <- NormalizeData(obj, assay = 'RNA')
  
  return(obj)
  
}


###options###
######################
option_list = list(
  make_option(c("-r", "--input.reference"),
              type="character",
              default=NULL, 
              help="path to integrated multiome object",
              metavar="character"),
  make_option(c("-q", "--input.query.object"),
              type="character",
              default=NULL, 
              help="path to reference rna object",
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
  make_option(c("-t", "--meta"),
              type="character",
              default="./", 
              help="metadata path",
              metavar="character"),
  make_option(c("-c", "--cell.column"),
              type="character",
              default="cell_type_v8.5_multi", 
              help="column with cell type to transfer",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
query.path <- opt$input.query.object
ref.path <- opt$input.reference
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$meta
cell_column <- opt$cell.column

dir.create(out_path, showWarnings = F)
setwd(out_path)


colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/Colors_panatac_v3.0.rds')
set.seed(666)
query.obj <- readRDS(query.path)
ref.obj <- readRDS(ref.path)

meta <- fread(meta.path) %>%
  column_to_rownames(var = 'V1') %>%    
  select(all_of(cell_column))
ref.obj <- AddMetaData(ref.obj, meta)


DefaultAssay(ref.obj) <- 'RNA'
  ref.obj <- DietSeurat(ref.obj, assays = 'RNA', counts = TRUE, data = FALSE)
ref.obj@meta.data$object.type <- 'Internal'
DefaultAssay(query.obj) <- 'RNA'
query.obj <- DietSeurat(query.obj, assays = 'RNA', counts = TRUE, data = FALSE)
query.obj@meta.data$object.type <- 'External'
query.obj@meta.data$Batches <- 'External'

new.obj <- merge(ref.obj, query.obj)

new.obj <- normalize_rna_harmony(new.obj, dims = 50, column = 'Batches')


p1 = DimPlot(new.obj, reduction = "rna.umap", group.by = "Batches", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(new.obj, reduction = "rna.umap", group.by = cell_column, label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p1+p2
ggsave(glue::glue('Dimplot_batches_cell_type_{add_filename}.pdf'), width = 12, height = 5)

saveRDS(new.obj, glue::glue('Merged_harmony_object_{add_filename}.rds'))
new.obj@meta.data %>% fwrite( glue::glue('Metadata_Merged_harmony_object_{add_filename}.tsv'), sep = '\t', row.names = T)

# 
# ref.obj <- DietSeurat(ref.obj, counts = FALSE, dimreducs = "pca")
# query.obj <- DietSeurat(query.obj, counts = FALSE, dimreducs = "ref.pca")
# 
# 
# ref.obj$id <- 'reference'
# query.obj$id <- 'query'
# refquery <- merge(ref.obj, query.obj)
# refquery[["pca"]] <- merge(ref.obj[["pca"]], query.obj[["ref.pca"]])
# refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50, reduction.name = "ref.query.umap", reduction.key = "refqueryUMAP_")
# refquery@meta.data[[cell_column]][is.na(refquery@meta.data[[cell_column]])] <- refquery$predicted.celltype.l1[refquery@meta.data[[cell_column]]]
# saveRDS(refquery, glue::glue('Reference_query_new_umap_object_{add_filename}.rds'))
# 
# 
# p1 <- DimPlot(refquery,reduction = "ref.query.umap", group.by = 'id', shuffle = TRUE)
# p2 <- DimPlot(refquery,reduction = cell_column, group.by = 'id', shuffle = TRUE)
# p1+p2
# ggsave(glue::glue('Dimplot_reference_query_merged_predicted.celltype.l1_{add_filename}.pdf'), width = 12, height = 5)




