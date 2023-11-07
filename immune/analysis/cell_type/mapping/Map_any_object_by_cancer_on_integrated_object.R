#Map reg RNA object on multiome object

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
library("glmGamPoi")



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

query.obj.list <- SplitObject(query.obj, split.by = 'Cancer')

query.obj.list <- lapply(X = query.obj.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- CellCycleScoring(x, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  x <- x %>% SCTransform(
    assay = 'RNA',
    method = "glmGamPoi",
    vars.to.regress =  c( "percent.mt", "S.Score", "G2M.Score"),
    conserve.memory = T,
    verbose = F,
    return.only.var.genes = F
  )
  return(x)
})


DefaultAssay(ref.obj) <- 'integrated'
ref.obj <- RunUMAP(ref.obj, nn.name = "weighted.nn", 
                   reduction.name = "wnn.umap", 
                   reduction.key = "wnnUMAP_", return.model = TRUE)
# 
# 
# 
# ref.obj <- FindMultiModalNeighbors(ref.obj, 
#                         reduction.list = list("pca", "lsi"), 
#                         dims.list = list(1:50, 2:50),
#                         knn.graph.name = "wsnn", 
#                         cache.index = TRUE,
#                         return.neighbor = TRUE,
#                         l2.norm = TRUE)

DefaultAssay(query.obj) <- 'SCT'

anchors <- list()
for (i in 1:length(query.obj.list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = ref.obj,
    query = query.obj.list[[i]],
    normalization.method = "SCT",
   # k.filter = NA,
    reference.reduction = "pca", 
    #reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

query.obj <- MapQuery(
  anchorset = anchors,
  query = query.obj,
  reference = ref.obj,
  refdata = list(
    celltype.l1 = cell_column
  ),
  reference.reduction = "pca", 
  reduction.model = "wnn.umap"
)

for (i in 1:length(query.obj)) {
  query.obj[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query.obj[[i]],
    reference = ref.obj, 
    refdata = list(
      celltype.l1 = cell_column
    ),
    reference.reduction = "pca", 
    reduction.model = "wnn.umap"
  )
}

names(query.obj.list) %>% walk(function(c) {
  p1 = DimPlot(query.obj.list[[c]], reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
  p2 = DimPlot(ref.obj, reduction = "wnn.umap", group.by = cell_column, label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
  p1+p2
  ggsave(glue::glue('Dimplot_{c}_predicted.celltype.l1_{add_filename}.pdf'), width = 12, height = 5)
  query.obj.list[[c]]@meta.data %>% fwrite( glue::glue('Metadata_{c}_mapped_object_{add_filename}.tsv'), sep = '\t', row.names = T)
  
})



ref.obj <- DietSeurat(ref.obj, counts = FALSE, dimreducs = "pca")
query.obj <- DietSeurat(query.obj, counts = FALSE, dimreducs = "ref.pca")


ref.obj$id <- 'reference'
query.obj$id <- 'query'
refquery <- merge(ref.obj, query.obj)
refquery[["pca"]] <- merge(ref.obj[["pca"]], query.obj[["ref.pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50, reduction.name = "ref.query.umap", reduction.key = "refqueryUMAP_")
refquery@meta.data[[cell_column]][is.na(refquery@meta.data[[cell_column]])] <- refquery$predicted.celltype.l1[refquery@meta.data[[cell_column]]]
saveRDS(refquery, glue::glue('Reference_query_new_umap_object_{add_filename}.rds'))


p1 <- DimPlot(refquery,reduction = "ref.query.umap", group.by = 'id', shuffle = TRUE)
p2 <- DimPlot(refquery,reduction = cell_column, group.by = 'id', shuffle = TRUE)
p1+p2
ggsave(glue::glue('Dimplot_reference_query_merged_predicted.celltype.l1_{add_filename}.pdf'), width = 12, height = 5)




