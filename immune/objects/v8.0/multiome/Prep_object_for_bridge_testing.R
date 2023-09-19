# Alla Karpova split existing multiome object into 60 by 40% for testing bridge integration

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
library("glmGamPoi")
library(harmony)



assign_batches <- function(obj) {
  if(opt$int_batch=='weird_Brca_Ov') {
    wierd.brca <- c('HT206B1-S1H4', 'HT378B1-S1H2')
    wierd.ov <- c('VF031V1-Tm1Y1', 'VF027V1-S1Y1', 'VF034V1-T1Y1')
    obj@meta.data$Batches <- case_when(obj$Piece_ID_RNA %in% wierd.brca ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Piece_ID_RNA %in% wierd.ov ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       
                                       TRUE ~ 'All_other')
  }  else if(opt$int_batch=='weird_Brca_Ov_PBMC') {
    wierd.brca <- c('HT206B1-S1H4', 'HT378B1-S1H2')
    wierd.ov <- c('VF031V1-Tm1Y1', 'VF027V1-S1Y1', 'VF034V1-T1Y1')
    obj@meta.data$Batches <- case_when(obj$Piece_ID_RNA %in% wierd.brca ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Piece_ID_RNA %in% wierd.ov ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Cancer == 'PBMC' ~  obj$Cancer,
                                       TRUE ~ 'All_other')
  }  else if(opt$int_batch=='weird_Brca_Ov_crc_PBMC') {
    wierd.brca <- c('HT206B1-S1H4', 'HT378B1-S1H2')
    wierd.ov <- c('VF031V1-Tm1Y1', 'VF027V1-S1Y1', 'VF034V1-T1Y1')
    wierd.crc <- c('HT307C1-Th1K1', 'HT270P1-S1H2')
    obj@meta.data$Batches <- case_when(obj$Piece_ID_RNA %in% wierd.brca ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Piece_ID_RNA %in% wierd.ov ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Piece_ID_RNA %in% wierd.crc ~  paste('CRC_PDAC', 'weird', sep = '__'),
                                       obj$Cancer == 'PBMC' ~  obj$Cancer,
                                       TRUE ~ 'All_other')
  } else if(opt$int_batch=='sample') {
    obj@meta.data$Batches <- obj@meta.data$Piece_ID_RNA
    
  } else if(opt$int_batch=='cancer') {
    obj@meta.data$Batches <- case_when(obj$Cancer=='GBM' ~ 'ccRCC',
                                       TRUE ~ obj$Cancer)
  } else if (opt$int_batch=='ov') {
    wierd.brca <- c('HT206B1-S1H4', 'HT378B1-S1H2')
    wierd.ov <- c('VF031V1-Tm1Y1', 'VF027V1-S1Y1', 'VF034V1-T1Y1')
    obj@meta.data$Batches <- case_when(obj$Piece_ID_RNA %in% wierd.brca ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Cancer == 'OV' ~  obj$Cancer,
                                       TRUE ~ 'All_other')
  } 
  return(obj)
}


integrate_rna <- function(obj) {
  
  obj <- assign_batches(obj)
  print(table(obj$Batches))
  
  cat ('Run SCT on batches\n')
  all.rna.list <- SplitObject(obj, split.by = 'Batches')
  
  batches <- names(all.rna.list)
  reference.batches.n <- which(grepl('All_other', batches))
  
  all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA'
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- CellCycleScoring(x, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
    regress.me <- c("percent.mt", "S.Score", "G2M.Score")
    
    x <- x %>% SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      vars.to.regress =  c( "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      verbose = F,
      return.only.var.genes = F
    )
    
    return(x)
  })
  
  message('Selecting integration features')
  features <- SelectIntegrationFeatures(object.list = all.rna.list, nfeatures = 4500)
  
  
  
  all.rna.list <- PrepSCTIntegration(object.list = all.rna.list, anchor.features = features)
  message('Run PCA on integration features')
  all.rna.list <- lapply(X = all.rna.list, FUN = RunPCA, features = features)
  message('Run FindIntegrationAnchors')
  
  rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, normalization.method = "SCT",
                                          anchor.features = features, dims = 1:50, reduction = "rpca")
    
  
  message('Run IntegrateData')
  int <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:50)
  
  
  int <- RunPCA(int, verbose = FALSE)
  int <- RunUMAP(int, reduction = "pca", dims = 1:50, reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  int <- PrepSCTFindMarkers(int)
  int <- NormalizeData(int, assay = 'RNA')
  return(int)
  
}

normalize_rna <- function(obj, dims=50) {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      vars.to.regress =  c( "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T,
      verbose = FALSE) %>% 
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims,reduction = 'pca', reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  return(obj)
}


normalize_atac <- function(obj, dims=50) {
  DefaultAssay(obj) <- "ATAC_immune"
  
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 500) %>%
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400) %>%
    RunUMAP(dims = 2:dims,reduction = 'lsi', reduction.name = "atac.umap", reduction.key = "atacUMAP_")
  return(obj)
}

normalize_multiome_with_integration <- function(obj,dims = 50) {
  obj <- integrate_rna(obj)
  obj <- normalize_atac(obj)
  obj <- FindMultiModalNeighbors(obj, 
                                 reduction.list = list("pca", "lsi"), 
                                 dims.list = list(1:dims, 2:dims))
  obj <- RunUMAP(obj, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 4,  resolution=1.2, verbose = T)
  
  return(obj)
}


###options###
######################
option_list = list(
  make_option(c("-m", "--input.multi.object"),
              type="character",
              default=NULL, 
              help="path to integrated multiome object",
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
  make_option(c("--do.integration"),
              type="character",
              default=FALSE,
              help = "Do RNA integartion or not",
              metavar="logical"),
  make_option(c("--int_batch"),
              type="character",
              default="chemistry",
              help = "options include 'weird_Brca_Ov','sample', 'cancer','ov'",
              metavar="character")
  
  
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.multi.object
out_path <- opt$output
add_filename <- opt$extra


dir.create(out_path, showWarnings = F)
setwd(out_path)


colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/Colors_panatac_v3.0.rds')

#if (!file.exists( glue::glue("PanImmune_merged_RNA_ATAC_{add_filename}.rds"))) {

set.seed(666)
r.obj <- readRDS(input.path)

for.multi <- r.obj@meta.data %>% group_by(seurat_clusters) %>% sample_frac(size = 0.6, replace = FALSE) %>% rownames()
for.test <- rownames(r.obj@meta.data[-for.multi,])

print(length(for.multi))
print(length(for.test))

r.obj <- DietSeurat(r.obj, assays = c('RNA', 'ATAC_immune'))
r.obj.multi <- subset(r.obj, cells = for.multi)

r.obj.test <- subset(r.obj, cells = for.test)
a.obj.test <- r.obj.test
a.obj.test[['RNA']] <- NULL
r.obj.test[['ATAC_immune']] <- NULL

rm(r.obj)
gc()
r.obj.multi <- normalize_multiome_with_integration(r.obj.multi)
saveRDS(r.obj.multi, glue::glue("PanImmune_int_RNA_ATAC_{add_filename}.rds"))

rm(r.obj.multi)
gc()
r.obj.test <- normalize_rna(r.obj.test)
saveRDS(r.obj.test, glue::glue("PanImmune_merged_RNA_{add_filename}.rds"))

rm(r.obj.test)
gc()
a.obj.test <- normalize_atac(a.obj.test)
saveRDS(a.obj.test, glue::glue("PanImmune_merged_ATAC_{add_filename}.rds"))








