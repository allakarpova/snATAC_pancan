## Alla Karpova
### recluster T-cells and NK cells by sub lineage and integrate RNA with seurat integration

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

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))



runAllNormalization <- function(obj, dims=30) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims, assay = 'SCT')
  
  obj <- NormalizeData(obj, assay = 'RNA')
  obj <- FindNeighbors(object = obj,  dims = 1:dims)
  obj <- FindClusters(object = obj,resolution = 2, verbose = FALSE)
  return(obj)
  
}

filter <- dplyr::filter
select <- dplyr::select


option_list = list(
  make_option(c("-i", "--input.obj"),
              type="character",
              default=NULL,
              help="path to input object",
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
  make_option(c("-t","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type.harmonized',
              help = "cell_type_column",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.obj
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column


dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


plan("multicore", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM


all.rna <- readRDS(input.path)
# add meta data if provided
if (!is.null(meta.path)) {
  my.metadata <- fread(meta.path, header = TRUE, data.table = F) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F) %>% select(all_of(cell_column))
  all.rna <- AddMetaData(all.rna, metadata = my.metadata)
}

DefaultAssay(all.rna) <- 'RNA'
all.rna <- all.rna %>%  DietSeurat(assay = 'RNA', counts = TRUE, data = TRUE)

all.rna$is_NK <- grepl('trNK|NK', as.character(unlist(all.rna[[cell_column]])))
all.rna$is_CD8 <- grepl('CD8|MAIT', as.character(unlist(all.rna[[cell_column]]))) & !grepl('REMOVE', as.character(unlist(all.rna[[cell_column]]))) 
all.rna$is_CD4 <- grepl('CD4', as.character(unlist(all.rna[[cell_column]]))) & !grepl('REMOVE', as.character(unlist(all.rna[[cell_column]]))) 

conditions <- c('is_NK', 'is_CD8', 'is_CD4') 

conditions %>% walk (function(column) {
  all.rna.sub <- subset(x = all.rna, cells = rownames(dplyr::filter(all.rna@meta.data, .data[[column]])))
  
  all.rna.sub$Data.source <- ifelse(all.rna.sub$Cancer == 'PBMC', '10x', 'DingLab')
  all.rna.sub$Batches <- case_when(all.rna.sub$Cancer %in% c('PBMC') ~ paste(all.rna.sub$Cancer, all.rna.sub$Chemistry, sep = '__'),
                                   all.rna.sub$Cancer %in% c('MM') ~ all.rna.sub$Cancer,
                                   TRUE ~ all.rna.sub$Chemistry)
  
  cat ('Integrate regular RNA and combo RNA by batches \n')
  all.rna.list <- SplitObject(all.rna.sub, split.by = 'Batches')
  
  all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA'
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- CellCycleScoring(x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
    x <- x %>% SCTransform(
      assay = 'RNA',
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    )
    return(x)
  })
  
  features <- SelectIntegrationFeatures(object.list = all.rna.list, nfeatures = 3000)
  all.rna.list <- PrepSCTIntegration(object.list = all.rna.list, anchor.features = features)
  all.rna.list <- lapply(X = all.rna.list, FUN = RunPCA, features = features)
  
  rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
  int <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:30)
  int <- RunPCA(int, verbose = FALSE)
  int <- RunUMAP(int, reduction = "pca", dims = 1:40)
  int <- FindNeighbors(int, reduction = "pca", dims = 1:40)
  int <- FindClusters(int, resolution = 2)
  
  cat('saving the object...\n')
  saveRDS(int,   paste0("PanImmune_", column,"_integrated_object_regularRNA_multiome_",add_filename,".rds"))
  
  dimplot=DimPlot(object = int, label = TRUE) + NoLegend()
  pdf(paste0("Dimplot_", column,"_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int, label = F, group.by = "Piece_ID")
  pdf(paste0("Dimplot_", column,"_Piece_ID_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int, label = F, group.by = "data.type")
  pdf(paste0("Dimplot_", column,"_data.type_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=14, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized.cancer")
  pdf(paste0("Dimplot_", column,"_cell_type.harmonized.cancer_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int, label = TRUE, group.by = cell_column)
  pdf(paste0("Dimplot_", column,"_",cell_column,"_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int, label = TRUE, group.by = "Cancer")
  pdf(paste0("Dimplot_", column,"_Cancer_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
})