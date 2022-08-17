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


###################################
option_list = list(
  make_option(c("-i", "--input.folder"),
              type="character",
              default=NULL, 
              help="path to folder with rds objects",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="foo", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-a", "--assay"),
              type="character",
              default="X500peaksMACS2", 
              help="which assay should be used to merge objects? X500peaksMACS2, peaks",
              metavar="character")#,
  # make_option(c("-s", "--samples.file"),
  #             type="character",
  #             default=NULL, 
  #             help="path to file with a list of samples in one column names 'Sample' and the second column named 'Data Type' indicating if its combo (10x_SC_Multi_ATAC_SEQ) or regular ATAC sample (snATAC)",
  #             metavar="character")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay
#sample.path <- opt$samples.file


# input.path <- '/home/allakarpova/ATAC_immune/obj/v3.0/chromvar/Reclustered_immune_snATAC_Merged_137_samples_v3.0_cluster_peaks_harmony_MM_data.type.chromvar.rds'
# out_path <- '/home/allakarpova/ATAC_immune/obj/v3.0_integrated'
# add_filename <- '137_samples_v3.0_integrated'
# assay.towork <- 'ATAC_immune'

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


atac <- readRDS(input.path)
DefaultAssay(atac) <- assay.towork
atac$data.type[atac$Cancer=='MM'] <- 'scATAC'

atac.split <- SplitObject(atac, split.by = 'data.type')

atac.split <- map(atac.split, function(obj) {
  obj <- FindTopFeatures(obj, min.cutoff = 10) %>%
    RunTFIDF() %>%
    RunSVD()
  return(obj)
})

atac <- atac %>% FindTopFeatures(min.cutoff = 10) %>%
  RunTFIDF() %>%
  RunSVD() %>%
  RunUMAP(reduction = "lsi", dims = 2:50)

#######integration############
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100 * 1024^3)
integration.anchors <- FindIntegrationAnchors(
  object.list = atac.split,
  anchor.features = rownames(atac),
  reduction = "rlsi",
  dims = 2:50
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50)
p2 <- DimPlot(integrated, group.by = "data.type")
p1 <- DimPlot(atac, group.by = "data.type")
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_data_type.pdf'), width = 13, height = 4.5)
saveRDS(integrated, paste0(add_filename, '_data_type.rds'))

p2 <- DimPlot(integrated, group.by = "Cell_type_state")
p1 <- DimPlot(atac, group.by = "Cell_type_state")
(p1 + ggtitle("Merged") + NoLegend()) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_cell_type_state.pdf'), width = 15, height = 4.5)

p2 <- DimPlot(integrated, group.by = "Cancer", cols = 'Spectral')
p1 <- DimPlot(atac, group.by = "Cancer", cols = 'Spectral')
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_Cancer.pdf'), width = 12, height = 4.5)


integrated  <-  integrated %>% 
  FindNeighbors(
    reduction = 'integrated_lsi',
    dims = 2:50
  ) %>% 
  FindClusters( 
    algorithm = 3,
    resolution = 1.2,
    verbose = FALSE
    )
saveRDS(integrated, paste0(add_filename, '_data_type.rds'))

########## REFERENCE mapping




