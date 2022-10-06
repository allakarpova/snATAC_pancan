# Do merging nofmalication of integrated ATAC object
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

dir.create(out_path, showWarnings = F)
setwd(out_path)


panc.my <- readRDS(input.path)
print(dim(panc.my))


####################################
##### Simple normalization with seurat #######
####################################

panc.my <- runAllNormalization(panc.my, dims = 40)

# #remove redundant fragment files
# all.fragment.obj <- Fragments(integrated)
# all.fragment.obj.cell.count <- map_chr(all.fragment.obj, function(x) length(x@cells))
# all.fragment.obj.upd <- all.fragment.obj[all.fragment.obj.cell.count > 0]
# Fragments(integrated) <- NULL
# Fragments(integrated) <- all.fragment.obj.upd

saveRDS(panc.my, paste0('PanImmune_merged_object_new_peaks_', add_filename, '.rds'))


p2 <- DimPlot(panc.my, group.by = "Chemistry")
p2 + ggtitle("Merged")
ggsave(paste0(add_filename, '_Chemistry.pdf'), width = 6, height = 4.5)

p2 <- DimPlot(panc.my, group.by = "cell_type.harmonized.cancer")
p2 + ggtitle("Merged")
ggsave(paste0(add_filename, '_cell_type.harm.cancer.pdf'), width = 6, height = 4.5)

p2 <- DimPlot(panc.my, group.by = "Cancer", cols = 'Spectral')
p2 + ggtitle("Merged")
ggsave(paste0(add_filename, '_Cancer.pdf'), width = 6, height = 4.5)



fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_merged_object_new_peaks_', add_filename, '_matadata.tsv'),
       sep='\t', row.names = T)



