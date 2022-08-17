#Author:Nadezhda Terekhanova
####2020-08-27
####This script performs the cell type annotation, using annotated RNA rds-object
####It will transfer the labels from the "predicted_cell_type" or "cell_type" (change accordingly) meta-data field of RNA RDS-object to the ATAC object, and will save it into the out/
####Upd: we use manually annotated RNA-objects, so the cell type annotation is in the "cell_type" meta.data field

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
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
args = commandArgs(trailingOnly=TRUE)

sample_atac=args[1]
sample_rna=args[2]

rna=readRDS(paste('../snRNA_processed/Inputs_v.2.0/',sample_rna,'_cell_types_annotated.rds',sep=""))

sample=sample_atac
print(sample)

atac=readRDS(paste("../2.Create_rds/out/",sample,'/',sample,'_processed_atac.rds',sep=""))

DefaultAssay(atac) <- 'RNA'

dir.create(paste("out/",sample,sep=""))

atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)

transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
#  refdata = rna$predicted_cell_type,
  refdata = rna$cell_type,
  weight.reduction = atac[['lsi']]
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)

###Making plots:

plot1 <- DimPlot(
  object = rna,
  group.by = 'cell_type',
#  group.by = 'predicted_cell_type',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(sample_rna,' snRNA-seq',sep=""))

plot2 <- DimPlot(
  object = atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(sample,' snATAC-seq',sep=""))

pdf(paste("out/",sample,"/",sample,"_snRNA_snATAC_integrated.pdf",sep=""),height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()

saveRDS(atac,paste("out/",sample,"/",sample,"_cellTyped.rds",sep=""))
