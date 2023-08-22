# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################
library(future)

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



################################

#####################################
####### FUNCTIONS ##################
####################################


runAllNormalization <- function(obj, dims) {
  #### run normalization to get initial clusters ###
  ########
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 500) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 ) %>% 
    FindNeighbors(
      reduction = 'lsi',
      dims = 2:dims ) %>% 
    FindClusters(
      algorithm = 3,
      resolution = 1,
      verbose = FALSE
    ) %>% 
    RunUMAP(dims = 2:dims,
            reduction = 'lsi')
  return(obj)
}

doIntegration <- function (int.sub.f,  k.w = 100, k.filter = 200) {
  int.sub.f$Data.source <- ifelse(int.sub.f$Cancer == 'PBMC', '10x', 'DingLab')
  if(conditions=='B-cell_Plasma') {
    int.sub.f@meta.data$Batches <- case_when(int.sub.f$Cancer %in% c('PBMC') ~ 'PBMC',
                                   int.sub.f$Cancer %in% c('MM') ~ int.sub.f$Cancer,
                                   TRUE ~ int.sub.f$Chemistry)
  }  else {
    int.sub.f@meta.data$Batches <- case_when(int.sub.f$Cancer %in% c('PBMC') ~ paste(int.sub.f$Cancer, int.sub.f$data.type, sep = '__'),
                                   int.sub.f$Cancer %in% c('MM') ~ int.sub.f$Cancer,
                                   TRUE ~ int.sub.f$Chemistry)
  }
  
  
  print(table(int.sub.f$Batches))
  
  atac.split <- SplitObject(int.sub.f, split.by = 'Batches')
  
  atac.split <- map(atac.split, function(obj) {
    obj <- FindTopFeatures(obj, min.cutoff = 500) %>%
      RunTFIDF() %>%
      RunSVD(reduction.key = 'LSI_',
             reduction.name = 'lsi',
             irlba.work = 400)
    return(obj)
  })
  
  #######integration############
  plan("multiprocess", workers = 10)
  options(future.globals.maxSize = 100 * 1024^3)
  
  cat('FindIntegrationAnchors\n')
  integration.anchors <- FindIntegrationAnchors(
    object.list = atac.split,
    anchor.features = rownames(int.sub.f),
    reduction = "rlsi",
    k.filter=k.filter,
    dims = 2:50
  )
  
  # integrate LSI embeddings
  cat('IntegrateEmbeddings\n')
  integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = int.sub.f[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:50, 
    k.weight = k.w
  )
  
  # create a new UMAP using the integrated embeddings
  integrated <- RunUMAP(integrated, 
                        reduction = "integrated_lsi", 
                        dims = 2:50)
  
  integrated  <-  integrated %>% 
    FindNeighbors(
      reduction = 'integrated_lsi',
      dims = 2:40
    ) %>% 
    FindClusters( 
      algorithm = 4,
      method='igraph',
      resolution = 1,
      verbose = FALSE
    )
  
  
  #Annotation(integrated) <- annotations.f
  return(integrated)
}

############################################

###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object",
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
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')


dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
int <- readRDS(input.path)

int$cell_lin <- case_when(grepl('NK|CD|Treg|Tfh|MAIT|ILC|dnT|gdT|T-cell', int$cell_type_v8_atac) ~ 'T-cells_NK',
                          grepl('DC|Macro|Micro|Mast|Mono', int$cell_type_v8_atac) ~ 'Myeloid',
                          grepl('B-cell|Plasma', int$cell_type_v8_atac) ~ 'B-cell_Plasma',
                          TRUE ~ int$cell_type_v8_atac)
#print(head(int@meta.data, n = 5))

conditions <- c( 'T-cells_NK', 'Myeloid','B-cell_Plasma')

conditions %>% walk (function(column) {
  print(column)
  if(!file.exists(paste0(column, '_object_same_peaks_normalized_MERGED_', add_filename, '.rds'))) {
    int.sub <- subset(x = int, subset = cell_lin == column)
    
    #normalize original object
    cat('normalize original object\n')
    int.sub <- int.sub %>% 
      FindTopFeatures(min.cutoff = 500) %>%
      RunTFIDF() %>%
      RunSVD() %>%
      RunUMAP(reduction = "lsi", dims = 2:50)
    saveRDS(int.sub,  paste0(column, '_object_same_peaks_normalized_MERGED_', add_filename, '.rds'))
  } else {
    int.sub <- readRDS(paste0(column, '_object_same_peaks_normalized_MERGED_', add_filename, '.rds'))
  }
  
  
  if (!file.exists(paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '.rds'))) {
    integrated <- doIntegration(int.sub, 
                                
                                k.w = ifelse(conditions=='B-cell_Plasma', 28, 100),
                                k.filter = ifelse(conditions=='B-cell_Plasma', 50, 200))
    
    saveRDS(integrated, paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '.rds'))
  } else {
    integrated <- readRDS(paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '.rds'))
  }
  
  p2 <- DimPlot(integrated, group.by = "Chemistry", label = TRUE)
  p1 <- DimPlot(int.sub, group.by = "Chemistry", label = TRUE)
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Chemistry.pdf'), width = 13, height = 4.5)
  #saveRDS(integrated, paste0(add_filename, '_Chemistry.rds'))
  
  p2 <- DimPlot(integrated, group.by = 'cell_type_v5_atac', label = TRUE)
  p1 <- DimPlot(int.sub, group.by = 'cell_type_v5_atac', label = TRUE)
  (p1 + ggtitle("Merged") + NoLegend()) | (p2 + ggtitle("Integrated"))
  ggsave(glue::glue('{column}_cell_type_v5_atac.pdf'), width = 15, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "Cancer", cols = 'Spectral', label = TRUE)
  p1 <- DimPlot(int.sub, group.by = "Cancer", cols = 'Spectral', label = TRUE)
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Cancer.pdf'), width = 12, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "Sample_type", cols = 'Paired', label = TRUE)
  p1 <- DimPlot(int.sub, group.by = "Sample_type", cols = 'Paired', label = TRUE)
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Sample_type.pdf'), width = 12, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "seurat_clusters", label = TRUE)
  p1 <- DimPlot(int.sub, group.by = "seurat_clusters", label = TRUE)
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_seurat_clusters.pdf'), width = 12, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "Batches", label = TRUE)
  p1 <- DimPlot(int.sub, group.by = "Batches", label = TRUE)
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Batches.pdf'), width = 12, height = 4.5)
  
  
  fwrite(integrated@meta.data, paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '_metadata.tsv'),
         sep='\t', row.names = T)
  
})

################









