# subset cells with TSS more than 3 or 4
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
#suppressMessages(library(harmony))
library(SeuratWrappers)


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

doIntegration <- function (int.sub.f, k.w = 100) {
  int.sub.f$Data.source <- ifelse(int.sub.f$Cancer == 'PBMC', '10x', 'DingLab')
  int.sub.f$Batches <- case_when(int.sub.f$Cancer %in% c('PBMC') ~ paste(int.sub.f$Cancer, int.sub.f$data.type, sep = '__'),
                                 int.sub.f$Cancer %in% c('MM') ~ int.sub.f$Cancer,
                                 TRUE ~ int.sub.f$Chemistry)
  
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
      algorithm = 3,
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
              metavar="character"),
  make_option(c("--tss.cut"),
              type="numeric",
              default=3,
              help = "Threshold for removing cells based on TSS enrichemnt score",
              metavar="numeric")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
tss.cut <- opt$tss.cut

dir.create(out_path, showWarnings = F)
setwd(out_path)


panc.my <- readRDS(input.path)
print(dim(panc.my))
panc.my <- subset(x = panc.my, subset = TSS.enrichemnt >= tss.cut)
print(dim(panc.my))

#remove redundant fragment files
all.fragment.obj <- Fragments(panc.my)
all.fragment.obj.cell.count <- map_chr(all.fragment.obj, function(x) length(x@cells))
all.fragment.obj.upd <- all.fragment.obj[all.fragment.obj.cell.count > 0]
Fragments(panc.my) <- NULL
Fragments(panc.my) <- all.fragment.obj.upd

panc.my <- panc.my %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 500) %>% 
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400 )
print(panc.my@reductions)


####################################
##### Integration with seurat #######
####################################

integrated <- doIntegration(panc.my, 
                            #annotations.f = annotations, 
                            k.w = ifelse(grepl('B-cell', add_filename), 70, 100))

#remove redundant fragment files
all.fragment.obj <- Fragments(integrated)
all.fragment.obj.cell.count <- map_chr(all.fragment.obj, function(x) length(x@cells))
all.fragment.obj.upd <- all.fragment.obj[all.fragment.obj.cell.count > 0]
Fragments(integrated) <- NULL
Fragments(integrated) <- all.fragment.obj.upd

saveRDS(integrated, paste0('PanImmune_seurat_integrated_object_new_peaks_', add_filename, '_chemistry_data.source.rds'))


p2 <- DimPlot(integrated, group.by = "Chemistry")
p2 + ggtitle("Integrated")
ggsave(paste0(add_filename, '_Chemistry.pdf'), width = 6, height = 4.5)

p2 <- DimPlot(integrated, group.by = "cell_type.harmonized.cancer")
p2 + ggtitle("Integrated")
ggsave(paste0(add_filename, '_cell_type.harm.cancer.pdf'), width = 6, height = 4.5)

p2 <- DimPlot(integrated, group.by = "Cancer", cols = 'Spectral')
p2 + ggtitle("Integrated")
ggsave(paste0(add_filename, '_Cancer.pdf'), width = 6, height = 4.5)



fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_seurat_integrated_object_new_peaks_', add_filename, '_chemistry_data.source_matadata.tsv'),
       sep='\t', row.names = T)









