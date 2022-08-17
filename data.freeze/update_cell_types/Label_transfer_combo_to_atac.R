suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))



option_list = list(
  make_option(c( "--combo.obj"),
              type="character",
              default=NULL, 
              help="path to combo RDS object",
              metavar="character"),
  make_option(c("-m", "--combo.meta"),
              type="character",
              default=NULL, 
              help="path to combo metadata file if needed",
              metavar="character"),
  make_option(c("-a", "--atac.obj"),
              type="character",
              default=NULL, 
              help="path to atac RDS object",
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
  make_option(c("-c", "--cell.type.column.name"),
              type="character",
              default="cell_type", 
              help="the name of the column containing cell type info in the RNA object",
              metavar="character"),
  make_option(c("-s", "--assay.atac"),
              type="character",
              default="X500peaksMACS2", 
              help="the name of ATAC assay in combo object",
              metavar="character")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
sample_combo <- opt$combo.obj
sample_atac <- opt$atac.obj
output <- opt$output
add_filename <- opt$extra
ct.column <- opt$cell.type.column.name
combo.meta <- opt$combo.meta
ass <-  opt$assay.atac

dir.create(output, showWarnings = F)
setwd(output)

combo.meta <- fread(combo.meta) %>% 
  data.frame(row.names = 1, check.rows=F, check.names=F)

combo.obj <- readRDS(sample_combo)
combo.obj <- AddMetaData(combo.obj, combo.meta[,c('orig.ident', ct.column)])
DefaultAssay(combo.obj) <- ass
combo.obj@assays$ATAC@key <- "atac_"

atac.obj <- readRDS(sample_atac)
DefaultAssay(atac.obj) <- ass
atac.obj@assays$ATAC@key <- "atac_"
atac.obj@assays$peaks@key <- 'peaks_'

plan("multicore", workers = 20)
options(future.globals.maxSize = 500 * 1024^3) # for 500 Gb RAM


combo.obj <- RunUMAP(combo.obj, reduction = "lsi", dims = 2:50, return.model = TRUE)

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = combo.obj,
  query = atac.obj,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:50
)

# map query onto the reference dataset
atac.obj <- MapQuery(
  anchorset = transfer.anchors,
  reference = combo.obj,
  query = atac.obj,
  refdata = as.character(combo.obj@meta.data[,ct.column]),
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)


p1 <- DimPlot(combo.obj, reduction = "umap", group.by = ct.column, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p2 <- DimPlot(atac.obj, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query")
p3 <- DimPlot(atac.obj, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query original UMAP")

p1 | p2 | p3
ggsave(paste0(add_filename, '_cell_type_predicted.pdf'), width = 18, height = 4.5)

fwrite(atac.obj@meta.data, paste0( add_filename, '_processed_multiomic_cellTyped.meta.data'), row.names = T, sep = '\t')




