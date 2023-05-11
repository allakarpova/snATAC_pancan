#library(SeuratDisk)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(tidyverse)
library(data.table)
library(doParallel)
library(future)
library(optparse)

###################################
option_list = list(
  make_option(c("-i", "--input.obj"),
              type="character",
              default=NULL, 
              help="path to PBMC rds object",
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
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay

out_path <- paste0(out_path, '/', add_filename)
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)

pbmc <- readRDS(input.path)

# load PBMC reference
reference <- LoadH5Seurat("/diskmnt/primary/published_data/seurat_vignettes_data/pbmc_multimodal.h5seurat")
reference@assays
DefaultAssay(pbmc) <- assay.towork


# transfer cell type labels from reference to query

transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = 'SCT',
  reference.reduction = 'spca',
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"
pbmc$cell_type <- pbmc$predicted.id
p1 <- DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "atac.umap")
p2 <- DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "wnn.umap")

pdf(paste0(add_filename, '_processed_', data.type, '_cell_types.pdf'), width = 13, height = 5)
p1 + p2
dev.off()



#saveRDS(pbmc, 'PBMC_granulocyte_sorted_10k_processed_multiomic_cell_typed.rds')
fwrite(pbmc@meta.data, paste0(add_filename,'_processed_',data.type,'_cellTyped.meta.data'), sep = '\t', row.names = T)

