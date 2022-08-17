# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################
library(future)

plan("multicore", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

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

runAllNormalization <- function(obj, dims) {
  #### run normalization to get initial clusters ###
  ########
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 20) %>% 
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

doIntegration <- function (int.sub.f, annotations.f) {
  ##########################
  ##### Integration #######
  ##########################
  
  atac.split <- SplitObject(int.sub.f, split.by = 'Chemistry')
  cat('normalize each split object\n')
  atac.split <- map(atac.split, function(obj) {
    obj <- FindTopFeatures(obj, min.cutoff = 10) %>%
      RunTFIDF() %>%
      RunSVD()
    return(obj)
  })
  
  #######integration############
  cat('integration\n')
  plan("multicore", workers = 10)
  options(future.globals.maxSize = 100 * 1024^3)
  
  integration.anchors <- FindIntegrationAnchors(
    object.list = atac.split,
    anchor.features = rownames(int.sub.f),
    reduction = "rlsi",
    dims = 2:50
  )
  
  # integrate LSI embeddings
  integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = int.sub.f[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:50
  )
  
  # create a new UMAP using the integrated embeddings
  
  integrated <- RunUMAP(integrated, 
                        reduction = "integrated_lsi", 
                        dims = 2:50)
  
  integrated  <-  integrated %>% 
    FindNeighbors(
      reduction = 'integrated_lsi',
      dims = 2:30
    ) %>% 
    FindClusters( 
      algorithm = 3,
      resolution = 2,
      verbose = FALSE
    )
  Annotation(integrated) <- annotations.f
  return(integrated)
}

############################################

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')


dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
int <- readRDS(input.path)



int$Cell_type_combo_reg_doublets <- case_when(int$Cancer == 'PBMC' ~ int$cell_type.harmonized.cancer,
                                              TRUE ~ int$Cell_type_combo_reg_doublets) 

tb <- table(int$seurat_clusters, int$Cell_type_combo_reg_doublets)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
int$Cell_type_combo_reg <- cluster.match.celltype[as.character(int$seurat_clusters)]

int$cell_lin <- case_when(grepl('NK|CD|Treg|MAIT|ILC|dnT|gdT', int$Cell_type_combo_reg) ~ 'T-cells_NK',
                           grepl('DC|TAM|Micro|Mast', int$Cell_type_combo_reg) ~ 'Myeloid',
                           grepl('B|Plasma', int$Cell_type_combo_reg) ~ 'B-cell_Plasma',
                           TRUE ~ int$Cell_type_combo_reg)
print(head(int@meta.data, n = 5))

conditions <- c('T-cells_NK', 'Myeloid', 'B-cell_Plasma')

conditions %>% walk (function(column) {
  print(column)
  if(!file.exists(paste0(column, '_object_same_peaks_normalized_MERGED_', add_filename, '.rds'))) {
    int.sub <- subset(x = int, subset = cell_lin == column)
    
    #normalize original object
    cat('normalize original object\n')
    int.sub <- int.sub %>% FindTopFeatures(min.cutoff = 10) %>%
      RunTFIDF() %>%
      RunSVD() %>%
      RunUMAP(reduction = "lsi", dims = 2:50)
    saveRDS(int.sub,  paste0(column, '_object_same_peaks_normalized_MERGED_', add_filename, '.rds'))
  } else {
    int.sub <- readRDS(paste0(column, '_object_same_peaks_normalized_MERGED_', add_filename, '.rds'))
  }
  
  if (!file.exists(paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '.rds'))) {
    integrated <- doIntegration(int.sub, annotations.f = annotations)
    saveRDS(integrated, paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '.rds'))
  } else {
    integrated <- readRDS(paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '.rds'))
  }
  
  p2 <- DimPlot(integrated, group.by = "Chemistry")
  p1 <- DimPlot(int.sub, group.by = "Chemistry")
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Chemistry.pdf'), width = 13, height = 4.5)
  #saveRDS(integrated, paste0(add_filename, '_Chemistry.rds'))
  
  p2 <- DimPlot(integrated, group.by = 'Cell_type_combo_reg')
  p1 <- DimPlot(int.sub, group.by = 'Cell_type_combo_reg')
  (p1 + ggtitle("Merged") + NoLegend()) | (p2 + ggtitle("Integrated"))
  ggsave(glue::glue('{column}_Cell_type_combo_reg.pdf'), width = 15, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "Cancer", cols = 'Spectral')
  p1 <- DimPlot(int.sub, group.by = "Cancer", cols = 'Spectral')
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Cancer.pdf'), width = 12, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "Sample_type", cols = 'Paired')
  p1 <- DimPlot(int.sub, group.by = "Sample_type", cols = 'Paired')
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_Sample_type.pdf'), width = 12, height = 4.5)
  
  p2 <- DimPlot(integrated, group.by = "seurat_clusters")
  p1 <- DimPlot(int.sub, group.by = "seurat_clusters")
  (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
  ggsave(paste0(column, '_seurat_clusters.pdf'), width = 12, height = 4.5)
  
  
  fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), paste0(column, '_object_same_peaks_normalized_INTEGRATED_', add_filename, '_metadata.tsv'),
         sep='\t', row.names = T)
  
})

################






