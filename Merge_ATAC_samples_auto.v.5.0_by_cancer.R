## merge ATAC samples which can be either regular ATAC seq or combo ATAC seq
#v 4.1 Alla implemented checking that all fragments files exist before starting merging
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

require(magrittr)
require(readr)
suppressMessages(library(Matrix))
suppressMessages(library(tidyr))
set.seed(1234)

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
library(optparse)
suppressMessages(library(data.table))

library(googlesheets4)
library(stringr)
suppressMessages(library(doParallel))

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

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 50)
options(future.globals.maxSize = 400 * 1024^3) # for 400 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::select(`Disease Type`, Sample, `Data Type`, `Data folder`)
cancers <- unique(samples$`Disease Type`)


for (cancer in cancers) {
  #samples <- fread(sample.path, data.table = F, header = F)
  samples.cancer <- base::subset(samples, `Disease Type`==cancer)
  samples.id <- samples.cancer$Sample %>% as.character()
  samples.type <- samples.cancer$`Data Type` %>% as.character()
  if(length(samples.id)==1) {
    next()
  }
  if(file.exists(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_",add_filename,".rds"))) {
    next()
  }
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
  
  #####To obtain the best results - use ALL peaks!
  cat ('Combining peaks\n')
  
  combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  peaks.use=combined.peaks
  
  #peaks.use=sample(combined.peaks, size = 8000, replace = FALSE)
  
  #We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.
  cat ('creating matrix counts\n')
  matrix.counts=vector(mode = "list", length = length(samples.id))
  for (i in 1:length(samples.id)){
    matrix.counts[[i]] <- FeatureMatrix(
      fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
      features = peaks.use,
      sep = c("-","-"),
      cells = colnames(atac[[i]])
    )
  }
  
  registerDoParallel(cores=10)
  cat ('creating peaksinters and removing useless assays\n')
  atac <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
    obj[['peaksinters']] <- CreateChromatinAssay(counts = co,fragments=Fragments(obj@assays$X500peaksMACS2), min.cells = -1, min.features = -1)
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
  DefaultAssay(combined) <- "peaksinters"
  
  #remove individual objects
  rm(atac)
  gc()
  
  cat ('Normalizing and reduction\n')
  combined <- RunTFIDF(combined)
  combined <- FindTopFeatures(combined, min.cutoff = 20)
  combined <- RunSVD(
    combined,
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400
  )
  
  combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
  combined <- FindNeighbors(
    object = combined,
    reduction = 'lsi',
    dims = 2:30
  )
  combined <- FindClusters(
    object = combined,
    algorithm = 3,
    resolution = 1,
    verbose = FALSE
  )
  
  cat('saving the object...\n')
  saveRDS(combined, paste0(cancer,'_', length(samples.id),"_snATAC_Merged_",add_filename,".rds"))
  
  
  n <- length(unique(combined$dataset))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  cat('plotting...\n')
  p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1, cols = col_vector) + 
    ggplot2::ggtitle("Combined snATAC samples.id")
  p2 <- DimPlot(combined, pt.size = 0.1,label=T) + 
    ggplot2::ggtitle("Combined snATAC clusters")
  p3 <- DimPlot(combined, group.by = 'Data.type', pt.size = 0.1, cols = 'Dark2') + 
    ggplot2::ggtitle("Combined snATAC Data.type")
  pdf(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_", add_filename, ".pdf"),height=6,width=24, useDingbats = F)
  print(p1+p2 +p3)
  dev.off()
  write.table(combined@meta.data, paste0(cancer,'_', length(samples.id),"_snATAC_Merged_",add_filename,"_metaData.txt"),
              sep="\t",quote=FALSE, row.names = T)
  
}

