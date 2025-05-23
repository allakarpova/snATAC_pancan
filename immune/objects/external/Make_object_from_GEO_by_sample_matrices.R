#!/usr/bin/env Rscript --vanilla
#
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
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. cellranger output's raw matrices folder)",
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
matrices.folder = opt$input
#meta = opt$metadata

# make output dir if it doesn't exist
out_path=paste(out_path, '/', sample_id,"/",sep="")
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)

if(!file.exists(paste0(sample_id, "_raw.rds"))) {
  # read in matrix
  list.of.objects <- list.files(path = matrices.folder, pattern = '.gz', full.names = T) %>%
    map(function(x) {
    
      file <- basename(x)
      print(file)
      sample.name <- str_split_fixed(str_sub(file, 12, -1 ), '[.]', 2)[,1]
      input <- fread(glue::glue("{matrices.folder}/{file}")) %>% 
        data.frame(row.names = 1, check.names = F, check.rows = F)
      obj = CreateSeuratObject(counts = input, project = sample.name, min.cells = 5)
      return(obj)
    })
  list.of.samples <- map(list.of.objects, function(x) {x@meta.data$orig.ident[1]})
  list.of.samples <- unlist(list.of.samples)
  panc = merge(list.of.objects[[1]], list.of.objects[-1])
  
  saveRDS(panc, file = paste0(sample_id, "_raw.rds"))
} else {
  panc <- readRDS(paste0(sample_id, "_raw.rds"))
}



### QC
# get percent mitochondrial content
cat('get percent mitochondrial content\n')
panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^MT-")

print(head(panc@meta.data))
# plot pre-filter metadata

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
if (max(panc[["percent.mt"]])>0) {
panc <- SCTransform(panc, 
                    vars.to.regress = c("percent.mt","S.Score", "G2M.Score"),return.only.var.genes = T)
} else{
  panc <- SCTransform(panc, 
                      vars.to.regress = c("S.Score", "G2M.Score"),return.only.var.genes = T)
  
}
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
