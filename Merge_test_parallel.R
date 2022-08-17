suppressMessages(library(doParallel))


suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
suppressMessages(library(reshape))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
library(future)
library(optparse)
library(data.table)

library(googlesheets4)
library(stringr)

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
              default="./", 
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

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)
dir.create('indiv_obj/')

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 50)
options(future.globals.maxSize = 450 * 1024^3) # for 400 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data folder`)
samples <- samples[15:20,]
samples.id <- samples$Sample %>% as.character()
samples.type <- samples$`Data Type` %>% as.character()

#########
cat (paste("Samples found:" ,length(samples.id), '\n'))
paths <- NULL
for (i in 1:length(samples.id)){
  print(samples.id[i])
  p <- list.files(path = input.path, full.names = T, pattern = paste0(str_split_fixed(samples.id[i], '_',2)[2],'.*rds'), all.files = T, recursive = T)
  print(length(p))
  paths <- c(paths, p)
}
#stop if not all samples have RDS object
print(length(samples.id))
print(length(paths))
stopifnot(length(samples.id)==length(paths))

# make the list of atac objects
registerDoParallel(cores=10)
cat ('Reading in objects\n')
#atac=vector(mode = "list", length = length(samples.id))
atac <- foreach (i=1:length(samples.id), p = paths, .combine=c) %dopar% {
  print(samples.id[i])
  obj=readRDS(p)
  DefaultAssay(obj) <- assay.towork
  if (!file.exists(Fragments(obj)[[1]]@path)) stop("Urgh, this sample object can't locate fragments file")
  obj<- DietSeurat(obj, assay = assay.towork)
  obj$dataset = samples.id[i]
  obj$Data.type = samples.type[i]
  return(obj)
}
stopImplicitCluster()

str(atac)

cat ('Reducing peaks\n')
combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
combined.peaks=combined.peaks
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
#peaks.use <- combined.peaks
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

registerDoParallel(cores=10)
cat ('creating matrix counts\n')
#matrix.counts=vector(mode = "list", length = length(samples.id))
matrix.counts <- foreach (obj = atac, .combine=c) %dopar% {
  FeatureMatrix(
    fragments = Fragments(obj@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(obj)
  ) 
}
stopImplicitCluster()

#str(atac)

registerDoParallel(cores=10)
cat ('creating peaksinters and removing useless assays\n')
atac <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
  obj[['peaksinters']] <- CreateChromatinAssay(counts = co,fragments=Fragments(obj@assays$X500peaksMACS2))
  #obj$dataset=samples.id[i]
  DefaultAssay(obj)<-'peaksinters'
  ###remove other assay
  obj[['X500peaksMACS2']]<-NULL
  return(obj)
}
stopImplicitCluster()

####Merging and reduction on old peaks
cat ('Merging\n')
combined <- merge(x = atac[[1]], y = atac[2:length(samples.id)], add.cell.ids = samples.id)
cat ('Done\n')
str(combined)
DefaultAssay(combined) <- "peaksinters"

#remove individual objects
rm(atac)
gc()

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_testing_5k_parallel_",add_filename,".rds"))



