#!/usr/bin/env Rscript --vanilla
# script to process smart-seq data into seurat object
# load required libraries
library(optparse)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(data.table)

# create user options
option_list = list(
  make_option(c("-i", "--input.tpm"),
              type="character",
              default=NULL,
              help="path to TPM matrix",
              metavar="character"),
  make_option(c("-i", "--input.raw"),
              type="character",
              default=NULL,
              help="path to raw count matrix",
              metavar="character"),
  make_option(c( "--metadata"),
              type="character",
              default=NULL,
              help="path to metadata file",
              metavar="character"),
  make_option(c("--nfeature_min"),
              type="integer",
              default=200,
              help="nFeature_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--nfeature_max"),
              type="integer",
              default=10000,
              help="nFeature_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--ncount_min"),
              type="integer",
              default=1000,
              help="nCount_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--ncount_max"),
              type="integer",
              default=80000,
              help="nCount_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--mito_max"),
              type="double",
              default=10,
              help="maximum allowed mitochondrial percentage",
              metavar="double"),
  make_option(c("-o", "--output"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-s","--sample_id"),
              type="character",
              default="single_cell_study",
              help="Name of your sample",
              metavar="character"),
  make_option(c("--pc_num"),
              type="integer",
              default=30,
              help="number of principal components to use",
              metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}

# read in initial arguments
sample_id <- opt$sample_id
out_path <- opt$output
matrix.tpm = opt$input.tpm
matrix.raw = opt$input.raw
meta = opt$metadata

# make output dir if it doesn't exist
out_path=paste(out_path, '/', sample_id,"/",sep="")
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


# read in matrix
cat('read in matrix\n')

  input.tpm <- fread(matrix.tpm, data.table = F) %>% 
    column_to_rownames(var='V1')
}

if(!is.null(matrix.raw)) {
  input.raw <- fread(matrix.raw, data.table = F) %>% 
    column_to_rownames(var='V1')
}

if(!is.null(meta) & !is.null(matrix.raw)) {
  metadata <- fread(meta, data.table = F) %>% 
    data.frame(row.names = 1) 
  panc = CreateSeuratObject(counts = input.raw, data = input.tpm, meta.data = metadata)
  
} else if (is.null(meta) & !is.null(matrix.raw)) {
  panc = CreateSeuratObject(counts = input.raw, data = input.tpm)
} else if(!is.null(meta) & is.null(matrix.raw)) {
  panc = CreateSeuratObject( data = input.tpm, meta.data = metadata)
} else {
  panc = CreateSeuratObject( data = input.tpm)
}

### QC
# get percent mitochondrial content
cat('get percent mitochondrial content\n')
panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^MT-")


# plot pre-filter metadata
#panc$percent.mito<-percent.mito
pdf(paste("QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# plot metadata associations
pdf(paste0("FeatureScatter_in_sample_",sample_id,".pdf",sep=""),width=12,height=7)
plot1 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter step
panc<-subset(x = panc, subset = nFeature_RNA > opt$nfeature_min & 
               nFeature_RNA < opt$nfeature_max & 
               nCount_RNA > opt$ncount_min & 
               nCount_RNA < opt$ncount_max & 
               percent.mt<opt$mito_max)
panc %>% dim

# plot post-filter metadata
pdf(paste("After_QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
dev.off()

# Run the standard workflow for visualization and clustering

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
panc <- NormalizeData(panc, assay = 'RNA')
panc <- CellCycleScoring(panc, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
print(head(panc@meta.data))
panc <- SCTransform(panc, 
                    vars.to.regress = c("percent.mt","S.Score", "G2M.Score"),return.only.var.genes = T)
panc <- RunPCA(panc, npcs = opt$pc_num, verbose = FALSE)

# t-SNE and Clustering
panc <- RunUMAP(panc, reduction = "pca", dims = 1:opt$pc_num)
panc <- FindNeighbors(panc, reduction = "pca", dims = 1:opt$pc_num)
panc <- FindClusters(panc,resolution = 1,algorithm = 4,
                     method='igraph',  verbose = FALSE)

# plot the clusters
pdf(paste0("DimPlot_",sample_id,".pdf"),useDingbats=FALSE)
DimPlot(object = panc, reduction = "umap",label=TRUE,label.size=6)
dev.off()


# save object so far
saveRDS(panc,file = paste(sample_id, "_processed.rds", sep=""))
