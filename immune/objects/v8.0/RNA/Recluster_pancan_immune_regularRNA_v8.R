# take integrated RNA object and keep only regular RNA cells - remove multiome RNA cells
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
library(harmony)



normalize_rna <- function(obj, dims=50) {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      vars.to.regress =  c( "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T,
      verbose = FALSE) %>% 
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims,reduction = 'pca', reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  return(obj)
}



###options###
######################
option_list = list(
  make_option(c("-r", "--input.rna.object"),
              type="character",
              default=NULL, 
              help="path to integrated RNA object",
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
  make_option(c("--metadata.rna"),
              type="character",
              default=NULL, 
              help="path to rna metadata"),
  make_option(c("--metadata.multi"),
              type="character",
              default=NULL, 
              help="path to multiomic metadata"),
  make_option(c("-c", "--cell.column"),
              type="character",
              default="cell_type_v8.5_multi", 
              help="column with cell type to transfer",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)



# read in initial arguments
input.path.rna <- opt$input.rna.object

out_path <- opt$output
add_filename <- opt$extra
meta.rna.path <- opt$metadata.rna
meta.multi.path <- opt$metadata.multi
cell_column <- opt$cell.column

r.obj <- readRDS(input.path.rna)
DefaultAssay(r.obj) <- 'RNA'
r.obj <- DietSeurat(r.obj, counts = TRUE,assays = 'RNA')

meta <- fread(meta.path) %>%
  column_to_rownames(var = 'V1') %>%    
  select(all_of(cell_column))
r.obj <- AddMetaData(r.obj, meta)


multi.cells <- fread(meta.multi.path) %>% pull(V1)
r.obj$Multiome <- rownames(r.obj@meta.data) %in% multi.cells

#remove multiome cells
r.obj <- subset(r.obj, Multiome, invert = TRUE)

r.obj <- r.obj %>% normalize_rna(r.obj)
saveRDS(r.obj, glue::glue('PanImmune_merged_regularRNA_{add_filename}.rds'))

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/Colors_panatac_v3.0.rds')


DimPlot(r.obj, group.by = "Cancer", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$Cancer)
ggsave(glue::glue("Dimplot_{add_filename}_Cancer.pdf"), height=8,width=11)

DimPlot(r.obj, group.by = "Piece_ID_RNA", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(glue::glue("Dimplot_{add_filename}_Piece_ID.pdf"), height=8,width=20)

DimPlot(r.obj, group.by = cell_column, label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(glue::glue("Dimplot_{add_filename}_{cell_column}.pdf"), height=8,width=15)





