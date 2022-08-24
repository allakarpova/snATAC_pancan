## Alla Karpova
### merge regular RNA data for immune cells and integrate with seurat

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

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
  my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F) %>% select(all_of(cell_column))
  all.rna <- AddMetaData(all.rna, metadata = my.metadata)
}
all.rna$to_remove <- grepl('oublet', as.character(unlist(all.rna[[cell_column]])))
DefaultAssay(all.rna) <- 'RNA'
all.rna <- all.rna %>%  DietSeurat(assay = 'RNA', counts = TRUE, data = TRUE)
all.rna <- subset(x = all.rna, subset = to_remove, invert = TRUE)


all.rna$Data.source <- ifelse(all.rna$Cancer == 'PBMC', '10x', 'DingLab')
all.rna$Batches <- case_when(all.rna$Cancer %in% c('PBMC') ~ paste(all.rna$Cancer, all.rna$Chemistry, sep = '__'),
                             all.rna$Cancer %in% c('MM') ~ all.rna$Cancer,
                             TRUE ~ all.rna$Chemistry)

cat ('Integrate regular RNA and combo RNA by batches \n')
all.rna.list <- SplitObject(all.rna, split.by = 'Batches')

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
int <- RunUMAP(int, reduction = "pca", dims = 1:30)
int <- FindNeighbors(int, reduction = "pca", dims = 1:30)
int <- FindClusters(int, resolution = 2)

tb <- table(int$seurat_clusters, int$cell_type.harmonized.cancer)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
int$cell_type.immune <- cluster.match.celltype[as.character(int$seurat_clusters)]

cat('saving the object...\n')
saveRDS(int,   paste0("PanImmune_integrated_object_regularRNA_multiome_",add_filename,".rds"))

dimplot=DimPlot(object = int, label = TRUE) + NoLegend()
pdf(paste0("Dimplot_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "Piece_ID")
pdf(paste0("Dimplot_Piece_ID_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "data.type")
pdf(paste0("Dimplot_data.type_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=14, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized.cancer")
pdf(paste0("Dimplot_cell_type.harmonized.cancer_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = cell_column)
pdf(paste0("Dimplot_",cell_column,"_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "Cancer")
pdf(paste0("Dimplot_Cancer_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

