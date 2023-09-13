# Alla Karpova create pancan combo object with RNA and ATAC

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


######################
### FUNCTIONS #####

integrate_rna <- function(obj) {
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
      wierd.crc <- c('HT307C1-Th1K1')
      obj@meta.data$Batches <- case_when(obj$Piece_ID_RNA %in% wierd.brca ~  paste(obj$Cancer, 'weird', sep = '__'),
                                         obj$Piece_ID_RNA %in% wierd.ov ~  paste(obj$Cancer, 'weird', sep = '__'),
                                         obj$Cancer == 'PBMC' ~  obj$Cancer,
                                         TRUE ~ 'All_other')
  } else if(opt$int_batch=='sample') {
    obj@meta.data$Batches <- obj@meta.data$Piece_ID_RNA
    
  } else if(opt$int_batch=='cancer') {
    obj@meta.data$Batches <- obj$Cancer
  } else if (opt$int_batch=='ov') {
    wierd.brca <- c('HT206B1-S1H4', 'HT378B1-S1H2')
    wierd.ov <- c('VF031V1-Tm1Y1', 'VF027V1-S1Y1', 'VF034V1-T1Y1')
    obj@meta.data$Batches <- case_when(obj$Piece_ID_RNA %in% wierd.brca ~  paste(obj$Cancer, 'weird', sep = '__'),
                                       obj$Cancer == 'OV' ~  obj$Cancer,
                                       TRUE ~ 'All_other')
  } 
  
  print(table(obj$Batches))
  
  cat ('Run SCT on batches\n')
  all.rna.list <- SplitObject(obj, split.by = 'Batches')
  
  batches <- names(all.rna.list)
  reference.batches.n <- which(grepl('All_other', batches))
  
  all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA'
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- CellCycleScoring(x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
    if(opt$regress.cc.interferon.genes) {
      x <- InterferonScoring(x)
      regress.me <- c("percent.mt", "S.Score", "G2M.Score", "Interferon.response")
    } else {
      regress.me <- c("percent.mt", "S.Score", "G2M.Score")
    }
    x <- x %>% SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      verbose = F,
      return.only.var.genes = T
    )
    return(x)
  })
  
  message('Selecting integration features')
  features <- SelectIntegrationFeatures(object.list = all.rna.list, nfeatures = 3500)
  if(opt$regress.cc.interferon.genes) {
    interferon.genes <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/v8.0/allRNA2/by_lineage/no_doublets/Int_chemistry_cancer_reference_multiome2/no_doublets2/Interferon.reponse.genes.tsv', data.table = F, header = T) %>%
      pull(V1)
    features <- unique(c(features, cc.genes$s.genes, cc.genes$g2m.genes, interferon.genes))
    features <- setdiff(features, c('MLF1IP', 'FAM64A', 'HN1'))
  }
  print(length(features))
  
  all.rna.list <- PrepSCTIntegration(object.list = all.rna.list, anchor.features = features)
  message('Run PCA on integration features')
  all.rna.list <- lapply(X = all.rna.list, FUN = RunPCA, features = features)
  message('Run FindIntegrationAnchors')
  if(opt$do.reference) {
    rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, reference = reference.batches.n,
                                          normalization.method = "SCT",
                                          anchor.features = features, dims = 1:50, reduction = "rpca")
    
  } else {
    rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, normalization.method = "SCT",
                                          anchor.features = features, dims = 1:50, reduction = "rpca")
    
  }
  message('Run IntegrateData')
  int <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:50)
  
  if(opt$regress.cc.interferon.genes) { 
    message('Run ScaleData')
    int <- ScaleData(int, 
                     vars.to.regress = c("S.Score", "G2M.Score", "Interferon.response"),
                     do.scale = FALSE,
                     do.center = TRUE)
  }
  
  int <- RunPCA(int, verbose = FALSE)
  int <- RunUMAP(int, reduction = "pca", dims = 1:50, reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  
  return(int)
  
}

normalize_rna <- function(obj, dims=50) {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
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

InterferonScoring <- function(obj) {
  interferon.genes <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/v8.0/allRNA2/by_lineage/no_doublets/Int_chemistry_cancer_reference_multiome2/no_doublets2/Interferon.reponse.genes.tsv', data.table = F, header = T) %>%
    pull(V1)
  print(interferon.genes)
  features <- list('Interferon.response' = interferon.genes)
  ctrl <- length(features[[1]])
  
  obj <- AddModuleScore(
    object = obj,
    features = features,
    name = 'Interferon.response',
    ctrl = ctrl)
  
  obj@meta.data <- obj@meta.data %>% 
    rename(Interferon.response=Interferon.response1)
  
  return(obj)
}

normalize_multiome <- function(obj,dims = 50) {
  obj <- normalize_rna(obj)
  obj <- normalize_atac(obj)
  obj <- FindMultiModalNeighbors(obj, 
                                 reduction.list = list("pca", "lsi"), 
                                 dims.list = list(1:dims, 2:dims))
  obj <- RunUMAP(obj, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 4, resolution=1.2, verbose = T)
  
  return(obj)
}

normalize_multiome_with_integration <- function(obj,dims = 50) {
  obj <- integrate_rna(obj)
  obj <- normalize_atac(obj)
  obj <- FindMultiModalNeighbors(obj, 
                                 reduction.list = list("pca", "lsi"), 
                                 dims.list = list(1:30, 2:dims))
  obj <- RunUMAP(obj, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 4,  resolution=1, verbose = T)
  
  return(obj)
}


###options###
######################
option_list = list(
  make_option(c("-r", "--input.rna.object"),
              type="character",
              default=NULL, 
              help="path to integrated RNA object",
              metavar="character"),
  make_option(c("-a", "--input.atac.object"),
              type="character",
              default=NULL, 
              help="path to integrated ATAC object",
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
  make_option(c("--metadata.rna"),
              type="character",
              default=NULL, 
              help="path to rna metadata"),
  make_option(c("--metadata.atac"),
              type="character",
              default=NULL, 
              help="path to atac metadata"),
  make_option(c("--do.integration"),
              type="character",
              default=FALSE,
              help = "Do RNA integartion or not",
              metavar="logical") ,
  make_option(c("--int_batch"),
              type="character",
              default="chemistry",
              help = "options include 'weird_Brca_Ov','sample', 'cancer','ov'",
              metavar="character"),
  make_option(c("--do.reference"),
              type="logical",
              default=FALSE,
              help = "Do reference integartion against all_other samples or not (int_batch must be weird.brca.ov or ov",
              metavar="logical"),
  make_option(c("--regress.cc.interferon.genes"),
              type="logical",
              default=FALSE,
              help = "Regress interferon and cell cycle genes from integrated assay",
              metavar="logical")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path.rna <- opt$input.rna.object
input.path.atac <- opt$input.atac.object
out_path <- opt$output
add_filename <- opt$extra
meta.rna.path <- opt$metadata.rna
meta.atac.path <- opt$metadata.atac


dir.create(out_path, showWarnings = F)
setwd(out_path)


colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/Colors_panatac_v2.0.rds')

#if (!file.exists( glue::glue("PanImmune_merged_RNA_ATAC_{add_filename}.rds"))) {
  

r.obj <- readRDS(input.path.rna)
a.obj <- readRDS(input.path.atac)

likely.combo.cells <- intersect(colnames(r.obj), colnames(a.obj)) 
print(length(likely.combo.cells))
atac.combo.cells <- rownames(filter(a.obj@meta.data, data.type=='10x_SC_Multi_ATAC_SEQ'))
combo.cells <- intersect(likely.combo.cells, atac.combo.cells)
print(length(combo.cells))

a.obj <- subset(a.obj, cells = combo.cells)
r.obj <- subset(r.obj, cells = combo.cells)


r.obj[['ATAC_immune']] <- a.obj[['ATAC_immune']]

if(opt$do.integration) {
  cat('normalize_multiome_with_integration \n')
  r.obj <- normalize_multiome_with_integration(r.obj)
} else {
  cat('normalizing all \n')
  r.obj <- normalize_multiome(r.obj)
  cat('done \n')
}


saveRDS(r.obj, glue::glue("PanImmune_merged_RNA_ATAC_{add_filename}.rds"))
#} else {
#  r.obj <- readRDS(glue::glue("PanImmune_merged_RNA_ATAC_{add_filename}.rds"))
#}
p1 <- DimPlot(r.obj, reduction = "rna.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(r.obj, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(r.obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(glue::glue("Dimplot_{add_filename}_seurat_clusters.pdf"),height=8,width=24)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()


DimPlot(r.obj, reduction = "wnn.umap", group.by = "Cancer", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$Cancer)
ggsave(glue::glue("Dimplot_{add_filename}_Cancer.pdf"), height=8,width=11)


DimPlot(r.obj, reduction = "wnn.umap", group.by = "Piece_ID", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(glue::glue("Dimplot_{add_filename}_Piece_ID.pdf"), height=8,width=20)



