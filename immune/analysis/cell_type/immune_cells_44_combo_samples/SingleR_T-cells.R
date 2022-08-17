library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(SingleR)
library(celldex)
library(tidyverse)
library(RColorBrewer)

###Load in seurat object
obj<-readRDS("~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/comboRNA/snRNA_combo_Merged_is_Tcell_lin_44_samples.rds")
DefaultAssay(obj) <- 'RNA'
add_filename <- 'is_Tcell_lin' 
dir.create('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells/singleR')
setwd('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/comboRNA_immune_cells/singleR')

VlnPlot(obj, features = 'nCount_RNA')
#Convert seurat object to single-cell format for singleR 
#https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
Seurat_Object_Diet <- DietSeurat(obj, graphs = "umap", assays = 'RNA',counts = T, scale.data = F) #https://github.com/satijalab/seurat/issues/4633
obj.sce <- as.SingleCellExperiment(Seurat_Object_Diet)

#Reference SingleR
#https://bioconductor.org/packages/3.13/data/experiment/vignettes/celldex/inst/doc/userguide.html

ref <- MonacoImmuneData()
#singleR usage
#http://bioconductor.org/books/devel/SingleRBook/using-the-classic-mode.html#annotating-the-test-dataset
pred <- SingleR(test = obj.sce, ref = ref, 
                labels = ref$label.fine, assay.type.test=1)
#colnames(pred)
ref_pred<-as.data.frame(pred)

#Extract columns of interest to add to seurat object
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "MonacoImmuneData"

##Add predictions to seurat object
#Add metadata to seurat object
obj <- AddMetaData(
  object = obj,
  metadata = ref_predictions)

col_vecor <- c(brewer.pal(n = 12, name = 'Paired'), brewer.pal(n = 12, name = 'Set3'), brewer.pal(n = 9, name = 'Dark2'), brewer.pal(n = 9, name = 'Set1'),)
DimPlot(obj, group.by = 'MonacoImmuneData', label = F, cols = col_vecor)
ggsave(paste0( "Dimplot_MonacoImmuneData_", add_filename, ".pdf"),height=7,width=11,useDingbats=FALSE)

table(obj$MonacoImmuneData) %>% sort

##### TRY another reference
ref <- BlueprintEncodeData()
pred <- SingleR(test = obj.sce, ref = ref, 
                labels = ref$label.fine, assay.type.test=1)
#colnames(pred)
ref_pred<-as.data.frame(pred)

#Extract columns of interest to add to seurat object
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "BlueprintEncodeData"

##Add predictions to seurat object
#Add metadata to seurat object
obj <- AddMetaData(
  object = obj,
  metadata = ref_predictions)
DimPlot(obj, group.by = 'BlueprintEncodeData', label = F, cols = col_vecor)
ggsave(paste0( "Dimplot_BlueprintEncodeData_", add_filename, ".pdf"),height=7,width=11,useDingbats=FALSE)


### TRY another reference
ref <- DatabaseImmuneCellExpressionData()
pred <- SingleR(test = obj.sce, ref = ref, 
                labels = ref$label.fine, assay.type.test=1)
#colnames(pred)
ref_pred<-as.data.frame(pred)

#Extract columns of interest to add to seurat object
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "DatabaseImmuneCellExpressionData"

##Add predictions to seurat object
#Add metadata to seurat object
obj <- AddMetaData(
  object = obj,
  metadata = ref_predictions)
DimPlot(obj, group.by = 'DatabaseImmuneCellExpressionData', label = F, cols = col_vecor)
ggsave(paste0( "Dimplot_DatabaseImmuneCellExpressionData_", add_filename, ".pdf"),height=7,width=11,useDingbats=FALSE)


DimPlot(obj, group.by = 'Sample_type', label = F, cols = 'Dark2')
ggsave(paste0( "Dimplot_Sample_type_", add_filename, ".pdf"),height=7,width=8,useDingbats=FALSE)
DimPlot(obj, group.by = 'Cancer', label = F, cols = 'Dark2')
ggsave(paste0( "Dimplot_Cancer", add_filename, ".pdf"),height=7,width=8,useDingbats=FALSE)
DimPlot(obj, group.by = 'cell_type_upd', label = F, cols = 'Dark2')
ggsave(paste0( "Dimplot_cell_type_upd_", add_filename, ".pdf"),height=7,width=8,useDingbats=FALSE)


tb <- table(obj$seurat_clusters, obj$DatabaseImmuneCellExpressionData)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
obj$DatabaseImmuneCellExpressionData_majority <- cluster.match.celltype[as.character(obj$seurat_clusters)]


DimPlot(obj, group.by = 'DatabaseImmuneCellExpressionData_majority', label = F, cols = col_vecor)
ggsave(paste0( "Dimplot_DatabaseImmuneCellExpressionData_majority_", add_filename, ".pdf"),height=7,width=8,useDingbats=FALSE)


