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

runAllNormalization <- function(obj, dims=30) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress =  c( "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    ) 
  remove.genes <- rownames(obj)[grepl("TRBV|TRAV",rownames(obj))]
  VariableFeatures(obj) <- setdiff(VariableFeatures(obj), remove.genes)
  obj@assays$SCT@scale.data[setdiff(VariableFeatures(obj), remove.genes),]
  
  obj <- obj %>% RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims, assay = 'SCT')
  
  obj <- NormalizeData(obj, assay = 'RNA')
  obj <- FindNeighbors(object = obj,  dims = 1:dims)
  obj <- FindClusters(object = obj,resolution = 0.5,algorithm = 4,
                      method='igraph',  verbose = FALSE)
  return(obj)
  
}


# create user options
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
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.obj
out_path <- opt$output
add_filename <- opt$extra

# make output dir if it doesn't exist
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


# read in matrix
cat('read in objects \n')

list.samples <- list.dirs(input.path,full.names = FALSE)[-1]

list.obj <- list.samples %>% map(function(sample) {
  obj <- readRDS(glue::glue('{input.path}/{sample}/{sample}_processed.rds'))
  obj@meta.data$sample <- sample
  vdj.path <- glue::glue('{input.path}/{sample}/{sample}_VDJ_filtered_contig_annotations.csv.gz')
  
  if(file.exists(vdj.path)) {
    vdj <- fread(vdj.path) %>%
      distinct(barcode, raw_clonotype_id) %>%
      column_to_rownames('barcode')
    
    obj <- AddMetaData(obj,vdj)
  }
  
  return(obj)
})

obj <- merge(list.obj[[1]], list.obj[-1], add.cell.ids = list.samples)
obj <- runAllNormalization(obj)


# plot the clusters
pdf(paste0("DimPlot_clusters.pdf"),useDingbats=FALSE)
DimPlot(object = obj, reduction = "umap",label=TRUE,label.size=6)
dev.off()

pdf(paste0("DimPlot_sample.pdf"),useDingbats=FALSE)
DimPlot(object = obj,group.by = sample, reduction = "umap",label=TRUE,label.size=6)
dev.off()

# save object so far
saveRDS(panc,file = paste0(add_filename, "_merged.rds"))
