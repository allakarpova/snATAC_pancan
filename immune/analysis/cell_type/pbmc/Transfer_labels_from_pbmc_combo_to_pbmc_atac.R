suppressMessages(library(SeuratDisk))
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(GenomicRanges))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(future))
suppressMessages(library(optparse))

###################################
option_list = list(
  make_option(c("-i", "--input.obj"),
              type="character",
              default=NULL, 
              help="path to PBMC rds object",
              metavar="character"),
  make_option(c("-r", "--reference"),
              type="character",
              default=NULL, 
              help="path to annotated combo rds object",
              metavar="character"),
  make_option(c("-m", "--metadata"),
              type="character",
              default=NULL, 
              help="path to metadata of reference object",
              metavar="character"),
  make_option(c("-c", "--cell.type.column"),
              type="character",
              default='cell_type', 
              help="cell.type.column",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="foo", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-a", "--assay"),
              type="character",
              default="SCT", 
              help="which assay should be used to transfer labels? X500peaksMACS2, assay.towork",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.obj
ref.path <- opt$reference
meta.path <- opt$metadata
ct <- opt$cell.type.column
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay

out_path <- paste0(out_path, '/', add_filename)
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)

pbmc <- readRDS(input.path)

# load PBMC reference
reference <- readRDS(ref.path)

metadata <- fread(meta.path) %>% data.frame(row.names = 1) %>% select(matches(ct))
reference <- AddMetaData(reference, metadata)
print(head(reference@meta.data))

DefaultAssay(pbmc) <- assay.towork
DefaultAssay(reference) <- assay.towork

reference <- RunUMAP(reference, reduction = "lsi", dims = 2:50, return.model = TRUE)

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc, 
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:50
)

# map query onto the reference dataset
pbmc <- MapQuery(
  anchorset = transfer.anchors,
  reference = reference,
  query = pbmc,
  refdata = reference@meta.data[,ct],
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)


# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"
pbmc$cell_type <- pbmc$predicted.id

p1 <- DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "atac.umap")

pdf(paste0(add_filename, '_processed_atac_cell_types.pdf'), width = 8, height = 5)
p1
dev.off()



#saveRDS(pbmc, 'PBMC_granulocyte_sorted_10k_processed_multiomic_cell_typed.rds')
fwrite(pbmc@meta.data, paste0(add_filename,'_processed_atac_cellTyped.meta.data'), sep = '\t', row.names = T)

