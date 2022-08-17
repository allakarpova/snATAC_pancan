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



######################
### FUNCTIONS #####

normalize_atac<- function(obj) {
  
  DefaultAssay(obj) <- "pancan"
  
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures( min.cutoff = 20) %>%
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400) %>%
    RunUMAP(dims = 2:30,reduction = 'lsi', reduction.name = "atac.umap", reduction.key = "atacUMAP_") %>%
    FindNeighbors( reduction = 'lsi', dims = 2:30) %>%
    FindClusters(resolution = 0.5, algorithm = 3, verbose = T)
  return(obj)
}


###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="input object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--metadata.atac"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/All_159_samples_metadata_data_freeze_v5.0.tsv', 
              help="path to atac metadata")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output
meta.atac.path <- opt$metadata.atac

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')
dir.create('plots')

annotation <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/Colors_panatac_v1.0.rds')

atac.meta <- fread(meta.atac.path) %>% data.frame(row.names = 1)

  cat('opening objects \n')
  a.obj <- readRDS(input_path)
  cat('done \n')
  
  Annotation(a.obj) <-annotation 
  
  a.obj <- AddMetaData(a.obj, atac.meta)
  
  a.obj$Case_ID_new <- str_replace(a.obj$Case_ID, pattern = 'C[1-2]$', 'C')
  
  a.obj$Case_ID_new %>% unique
  
  cases.to.use <- a.obj@meta.data %>% 
    dplyr::select(Case_ID_new, Piece_ID) %>% distinct() %>%
    group_by(Case_ID_new) %>% 
    tally() %>% 
    filter(n > 1) %>% 
    pull(Case_ID_new) %>% unique()
  
  print(cases.to.use)
  
  a.obj <- subset(a.obj, Case_ID_new %in% cases.to.use)
  a.obj <- subset(a.obj, cell_type.harmonized.cancer %in% c('Tumor', 'Normal epithelial cells', 'Secretory Endometrial epithelial cells', 'Ciliated Endometrial epithelial cells'))
  
  a.obj.split <- SplitObject(a.obj, split.by = 'Case_ID_new')
  
  
  cat('normalizing all \n')
  a.obj.split %>% walk ( function(obj) {
    case <- obj$Case_ID_new[1]
    print(case)
    obj <- normalize_atac(obj)
    saveRDS(obj, paste0("out/",case,"_PrimaryMet_TumorNormal.", format(Sys.Date(), format="%Y%m%d"),".rds"))
    
    DimPlot(obj, reduction = "atac.umap", group.by = "cell_type.harmonized.cancer", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal.s_cell_type.harmonized.cancer.pdf"), height=8,width=11)
    
    DimPlot(obj, reduction = "atac.umap", group.by = "Piece_ID", label = TRUE, label.size = 2.5, repel = TRUE)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal._Piece_ID.pdf"), height=8,width=11)
    
    DimPlot(obj, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal._seurat_clusters.pdf"), height=8,width=11)
    
    obj <- FindClusters(obj, algorithm = 3, resolution = 0.1)
    DimPlot(obj, reduction = "atac.umap", group.by = 'pancan_snn_res.0.1', label = TRUE, label.size = 2.5, repel = TRUE)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal._pancan_snn_res.0.1.pdf"), height=8,width=11)
    
    obj <- FindClusters(obj, algorithm = 3, resolution = 0.2)
    DimPlot(obj, reduction = "atac.umap", group.by = 'pancan_snn_res.0.2', label = TRUE, label.size = 2.5, repel = TRUE)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal._pancan_snn_res.0.2.pdf"), height=8,width=11)
    
    obj <- FindClusters(obj, algorithm = 3, resolution = 0.3)
    DimPlot(obj, reduction = "atac.umap", group.by = 'pancan_snn_res.0.3', label = TRUE, label.size = 2.5, repel = TRUE)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal._pancan_snn_res.0.3.pdf"), height=8,width=11)
    
    obj <- FindClusters(obj, algorithm = 3, resolution = 0.4)
    DimPlot(obj, reduction = "atac.umap", group.by = 'pancan_snn_res.0.4', label = TRUE, label.size = 2.5, repel = TRUE)
    ggsave(paste0("plots/",case,"_PrimaryMet_TumorNormal._pancan_snn_res.0.4.pdf"), height=8,width=11)
    
    obj <- FindClusters(obj, algorithm = 3, resolution = 0.5)
  })
  cat('done \n')










