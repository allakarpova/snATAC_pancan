suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))

runNormalization <- function(obj, dims) {
  #### run normalization 
  ########
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 500) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 )
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
  make_option(c("-t","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

my.metadata <- fread(meta.path, data.table = F, header=TRUE) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)


panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)

panc.my$to_split <- is.na(unlist(panc.my[[cell_column]]))
panc.my$to_split[panc.my$Chemistry=='snATAC'] <- TRUE

panc.my.labeled <- subset(panc.my, to_split, invert=TRUE)
panc.my.unlabeled <- subset(panc.my, to_split)

cat('Integrate objects with RNA labels and without\n')
panc.my.labeled <- runNormalization(panc.my.labeled)
int.labeled <- doIntegration(panc.my.labeled, k.w=100)
saveRDS(int.labeled, paste0('PanImmune_comboATAC_with_comboRNA_labels_integrated_', add_filename, '.rds'))

DimPlot(int.labeled, reduction = "umap", group.by = cell_column, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference labeled")
ggsave(paste0('Dimplot_comboATAC_with_comboRNA_labels_integrated_', add_filename, '_',cell_column,'.pdf'), width = 6.5, height = 6)

panc.my.unlabeled <- runNormalization(panc.my.unlabeled)
int.unlabeled <- doIntegration(panc.my.unlabeled, k.w=100)
saveRDS(int.unlabeled, paste0('PanImmune_all_other_ATAC_unlabeled_integrated_', add_filename, '.rds'))


# compute UMAP and store the UMAP model
cat('Run UMAP on integrated labeled object and save the model\n')
int.labeled <- RunUMAP(int.labeled, reduction = "integrated_lsi", dims = 2:50, return.model = TRUE)

plan("multicore", workers = 20)
options(future.globals.maxSize = 500 * 1024^3) # for 500 Gb RAM

# find transfer anchors
cat('Run FindTransferAnchors\n')
transfer.anchors <- FindTransferAnchors(
  reference = int.labeled,
  query = int.unlabeled, 
  reference.reduction = "integrated_lsi",
  reduction = "lsiproject",
  dims = 2:50
)

cat('Run MapQuery\n')
# map query onto the reference dataset
int.unlabeled <- MapQuery(
  anchorset = transfer.anchors,
  reference = int.labeled,
  query = int.unlabeled,
  refdata =  as.charater(int.labeled@meta.data[,cell_column]),
  reference.reduction = "integrated_lsi",
  new.reduction.name = "ref.int.lsi",
  reduction.model = 'umap'
)


p1 <- DimPlot(int.labeled, reduction = "umap", group.by = cell_column, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p2 <- DimPlot(int.unlabeled, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query")

p1 | p2
ggsave(paste0('Dimplot_labeled_transferred_', add_filename, '_cell_type_predicted.pdf'), width = 12, height = 4.5)

p1 <- DimPlot(combo.obj, reduction = "umap", group.by = "Cancer",cols = colors$Cancer, label = TRUE, repel = TRUE, pt.size = 0.0005)  + ggtitle("Reference")
p2 <- DimPlot(int.obj, reduction = "ref.umap", group.by = "Cancer", cols = colors$Cancer, label = TRUE, repel = TRUE)  + ggtitle("Query")
p1 | p2
ggsave(paste0(add_filename, '_cancer.pdf'), width = 12, height = 4.5)



saveRDS(int.obj, paste0(add_filename, '_10142021.rds'))



