# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################
library(future)

#plan("multicore", workers = 20)
#options(future.globals.maxSize = 100 * 1024 ^ 3)

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
library(rliger)

################################

#####################################
####### FUNCTIONS ##################
####################################


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
#cell_column <- opt$cell_type_column
#macs2_path <- opt$macs2_path
#chrom.size <- opt$chrom_size
#meta.path <- opt$metadata.file
#assay.towork <- opt$assay
#annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

# input.path <-'/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/137Samples_PanCan_merged_obj/137_snATAC_113K_peaks_diffPCs.motifsAdded.chromvar.20210826.rds.gz'
# out_path <- '/diskmnt/Projects/snATAC_analysis/immune/obj/v3.0'
# add_filename <- '137_samples_v3.0_cluster_peaks'
# cell_column <-'cell_type.harmonized.cancer'
# macs2_path <- '/diskmnt/Projects/Users/allakarpova/Tools/anaconda3/envs/signac/bin/macs2'
# chrom.size <- '/diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt'
# meta.path <- '/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/50PCs/All_137_samples_metadata_data_freeze_v3.1.tsv'
# assay.towork <- 'pancan'

#out_path <- '/diskmnt/Projects/snATAC_analysis/immune/obj/v1.0'
#add_filename <- '100_sample_obj.v1.0_old_macs2'
dir.create(out_path, showWarnings = F)
setwd(out_path)

panc.my <- readRDS(input.path)

####################################
##### Integration with LIGER #######
####################################


panc.my$Data.source <- ifelse(panc.my$Cancer == 'PBMC', '10x', 'DingLab')
panc.my$Batches <- case_when(panc.my$Cancer %in% c('PBMC') ~ paste(panc.my$Cancer, panc.my$data.type, sep = '__'),
                             panc.my$Cancer %in% c('MM') ~ panc.my$Cancer,
                             TRUE ~ panc.my$Chemistry)

#plan("multicore", workers = 10)
#options(future.globals.maxSize = 100 * 1024^3)

#integrated <- RunFastMNN(object.list = SplitObject(panc.my, split.by = "Batches"))
panc.my <- FindTopFeatures(panc.my, min.cutoff = 'q20')
panc.my <- ScaleData(panc.my, split.by = "Batches", do.center = FALSE, block.size = 100, 
                     features = VariableFeatures(panc.my))
integrated <- RunOptimizeALS(panc.my, k = 20, lambda = 5, split.by = "Batches")
integrated <- RunQuantileNorm(integrated, split.by = "Batches")

integrated <- integrated %>%
  FindNeighbors(reduction = "iNMF", dims = 2:50) %>%
  FindClusters(verbose = TRUE, resolution = 2) %>%
  RunUMAP(reduction = "iNMF", dims = 2:50)

print(integrated@reductions)
saveRDS(integrated, paste0('PanImmune_LIGER_integrated_object_new_peaks_', add_filename, '_chemistry_data.source.rds'))

fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_LIGER_integrated_object_new_peaks_', add_filename, '_chemistry_data.source_matadata.tsv'),
       sep='\t', row.names = T)


p2 <- DimPlot(integrated, group.by = "Chemistry")
p1 <- DimPlot(panc.my, group.by = "Chemistry")
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_Chemistry.pdf'), width = 13, height = 4.5)
#saveRDS(integrated, paste0(add_filename, '_Chemistry.rds'))

p2 <- DimPlot(integrated, group.by = "cell_type.harmonized.cancer")
p1 <- DimPlot(panc.my, group.by = "cell_type.harmonized.cancer")
(p1 + ggtitle("Merged") + NoLegend()) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_cell_type.harm.cancer.pdf'), width = 15, height = 4.5)

p2 <- DimPlot(integrated, group.by = "Cancer", cols = 'Spectral')
p1 <- DimPlot(panc.my, group.by = "Cancer", cols = 'Spectral')
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_Cancer.pdf'), width = 12, height = 4.5)


saveRDS(integrated, paste0('PanImmune_LIGER_integrated_object_new_peaks_', add_filename, '_chemistry_data.source.rds'))

fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_LIGER_integrated_object_new_peaks_', add_filename, '_chemistry_data.source_matadata.tsv'),
       sep='\t', row.names = T)









