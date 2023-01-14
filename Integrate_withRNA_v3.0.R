#Author:Nadezhda Terekhanova
####2020-08-27
####This script performs the cell type annotation, using annotated RNA rds-object
####It will transfer the labels from the "predicted_cell_type" or "cell_type" (change accordingly) meta-data field of RNA RDS-object to the ATAC object, and will save it into the out/
####Upd: we use manually annotated RNA-objects, so the cell type annotation is in the "cell_type" meta.data field
### Adapted by Alla. Added parameters and the metadata will be saved. 

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
 
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
#library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)

library(optparse)
library(data.table)

option_list = list(
  make_option(c("-r", "--rna.obj"),
              type="character",
              default=NULL, 
              help="path to rna RDS object",
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
  make_option(c("-m", "--metadata"),
              type="character",
              default=NULL, 
              help="path to metadata for RNA oject if needed",
              metavar="character"),
  make_option(c("-c", "--cell.type.column.name"),
              type="character",
              default="peaks", 
              help="the name of the column containing cell type info in the RNA object",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
sample_rna <- opt$rna.obj
sample_atac <- opt$atac.obj
out_path <- opt$output
add_filename <- opt$extra
ct.column <- opt$cell.type.column.name
meta.path <- opt$metadata

dir.create(out_path, showWarnings = F)
setwd(out_path)


rna=readRDS(sample_rna)
rna@assays$ATAC@key <- 'atac_'
atac=readRDS(sample_atac)
atac@assays$ATAC@key <- 'atac_'

DefaultAssay(rna) <- 'RNA'
rna <- NormalizeData(rna,assay = 'RNA')

DefaultAssay(atac) <- 'ATACGeneActivity'


if(!is.null(meta.path)) {
  meta <- fread(meta.path, data.table = F, header = T) %>%
    data.frame(row.names = 1)
  rna <- AddMetaData(rna,meta)
}
atac <- NormalizeData(
  object = atac,
  assay = 'ATACGeneActivity',
  normalization.method = 'LogNormalize'#,
  #scale.factor = median(atac$nCount_ATACGeneActivity)
)

cat('doing FindTransferAnchors\n')
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  features = VariableFeatures(object = rna), 
  reference.assay = "RNA", 
  query.assay = "ATACGeneActivity", 
  reduction = 'cca'
)

cat('doing TransferData\n')
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  reference = rna,
  refdata = ct.column,
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

print(head(predicted.labels))

cat('Adding metadata\n')
atac <- AddMetaData(object = atac, metadata = predicted.labels)


cat('saving shit\n')
saveRDS(atac,paste0(add_filename, "_cellTyped.rds"))
fwrite(atac@meta.data, paste0(add_filename, '_snATAC_cellTyped.meta.data'), row.names = T, sep = '\t')

###Making plots:
cat('plotting plot1\n')
plot1 <- DimPlot(
  object = rna,
  group.by = ct.column,
  label = TRUE,
  repel = TRUE) + ggtitle(paste(add_filename,'snRNA-seq'))

cat('plotting plot2\n')
plot2 <- DimPlot(
  object = atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(add_filename,'snATAC-seq'))

cat('combining plots\n')
pdf(paste0(add_filename,"_snRNA_snATAC_integrated.pdf"),height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()


