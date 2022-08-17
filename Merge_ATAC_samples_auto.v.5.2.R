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
dir.create('indiv_obj/')

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 100 * 1024^3) # for 100 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data folder`)

samples.id <- samples$Sample %>% as.character()
samples.type <- samples$`Data Type` %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))

#########################################################################
#if the object has been merged, just open it
if (file.exists(paste0(length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))) {
  cat('opening the object...\n')
  combined <- readRDS(paste0(length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))
  
} else {
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
  
  checking.n.cells <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
    return(ncol(obj)==ncol(co))
    #stopifnot(ncol(obj)==ncol(co))
  }
  names(checking.n.cells) <- samples.id
  #print(checking.n.cells)
  print(checking.n.cells[!checking.n.cells])
  stopImplicitCluster()
  #str(atac)
  
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
  
  cat('saving the object...\n')
  saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))
  
}
#########################################################################

#########################################################################
# if peaks were already taken care of, load this file
if(file.exists(paste0(length(samples.id),'_recentered_final.filtered.',add_filename,'.tsv'))) {
  recentered_final <- fread(paste0(length(samples.id),'_recentered_final.filtered.',add_filename,'.tsv'), data.table = F)
} else {
  # read in peaks
  cat('work on peaks...\n')
  all_peaks <- lapply(samples.id, FUN = function (sample) {
    peaks=fread(paste0(input.path,'/', sample,"/recentered_final.filtered",sample,".tsv"))
    peaks$Sample=sample
    peaks$new_peak = paste(peaks$seqnames, peaks$recentered_start, peaks$recentered_end, sep = '-')
    return(peaks)
  })
  all_peaks <- rbindlist(all_peaks) # peaks from all samples
  fwrite(all_peaks, paste0(length(samples.id),"_sample_MACS2_peaks_",add_filename,".tsv"),
         sep='\t',row.names=FALSE)

  recentered_p=StringToGRanges(all_peaks$new_peak, sep = c("-", "-"))
  cat('finding overlapping peaks \n')
  overlapping=as.data.table(findOverlaps(recentered_p,recentered_p)) # find which peaks overlap
  overlapping=overlapping[queryHits!=subjectHits,] # remove hits that show peaks that overlap themselves
  fwrite(overlapping, 'tmp_files/139_overlapping_MACS2_peak_number.tsv', sep = '\t')
  overlapping.peak.number <- unique(overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks table

  recentered_non_overlapping=all_peaks[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
  fwrite(recentered_non_overlapping,paste0(length(samples.id),'_recentered_nonOverlapping.filtered.',add_filename,'.tsv'),
         sep='\t',row.names=FALSE)
  rm(recentered_p)
  gc()
  
  #split overlapping peaks by chromosome and process 500k peaks at a time
  dir.create('tmp_files', showWarnings = F)
  tmp <- data.table(chr = all_peaks$seqnames[overlapping.peak.number], 
                    num = overlapping.peak.number)
  overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) # overlapping peak number split by chromosome
  
  scores <- all_peaks$score # save scores for each peak for faster access in the following code
  best_in_overlapping <- lapply (X = overlapping.peak.number.split, 
                                 FUN = function(overlapping.peak.number) { # apply function on each chromosome peaks, overlapping.peak.number contains positions of overlapping peaks in all_peaks table
                                   best_in_overlapping_chr=NULL
                                   chr = overlapping.peak.number$chr[1] # just extract which chromosome we are working with
                                   cat(paste0(chr, '\n'))
                                   cat(paste0(nrow(overlapping.peak.number), '\n'))
                                   
                                   # I'm going to iterate over 500k peaks at a time and do them in parallel. 
                                   registerDoParallel(cores=25)
                                   bin.num <- nrow(overlapping.peak.number)%/%500000 # this is how many times I will iterate over 500k peaks
                                   best_in_overlapping_num <- vector(mode = "list", length = bin.num) # pre-allocate output for unique and the best peaks 
                                   for (b in 0:bin.num) {
                                     cat(paste0('Parallelization round #', b, '\n'))
                                     vector.toiterate <- (b*500000+1):((b+1)*500000) # peaks numbers to iterate over
                                     vector.toiterate <- vector.toiterate[vector.toiterate<=nrow(overlapping.peak.number)] # remove peak numbers that are larger then the total number of peaks 
                                     best_in_overlapping_num_sub <- foreach (i=vector.toiterate, .combine = c) %dopar% { # do paralellization
                                       if (i%%100000 == 0) {
                                         print(i) #this is to track the process, print i every 100k peaks
                                       }
                                       n <- overlapping.peak.number$num[i] # position of a peak that has overlaps
                                       neighbor.peaks.num <- c(n,overlapping[queryHits==n, subjectHits]) #find positions of other peaks overlapping with the first one 
                                       max.peak.pos <- which.max(scores[neighbor.peaks.num]) #among them select the position of a peak with the highest score
                                       return(neighbor.peaks.num[max.peak.pos]) # return this position
                                     }
                                     cat('parallelization is over\n')
                                     best_in_overlapping_num[[(b+1)]] <- best_in_overlapping_num_sub #save the vector of positions of the best peaks
                                   }
                                   best_in_overlapping_num <- do.call('c', best_in_overlapping_num) #combine best peaks positions from all iterations
                                   stopImplicitCluster()
                                   ######
                                   best_in_overlapping_num <- unique(best_in_overlapping_num) #keep unique peak numbers, this will cut the length of the vector significantly
                                   best_in_overlapping_chr <- all_peaks[best_in_overlapping_num,] #get the actual peaks for a chromosome
                                   fwrite(best_in_overlapping_chr,paste0('tmp_files/','139','_recentered_Overlapping.filtered.',add_filename,'_by_chr_', chr,'.tsv'),sep='\t',
                                          row.names=FALSE) # and save them
                                   return(best_in_overlapping_chr)
                                 })
  best_in_overlapping <- rbindlist(best_in_overlapping)  #combine best peaks positions from all chromosomes
  fwrite(best_in_overlapping,paste0('139','_recentered_Overlapping.filtered.',add_filename,'all.tsv'),sep='\t',
         row.names=FALSE)
  
  
  recentered_final=rbindlist(list(recentered_non_overlapping,best_in_overlapping))
  fwrite(recentered_final,paste0(length(samples.id),'_recentered_final.filtered.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  rm(recentered_non_overlapping, best_in_overlapping)
  gc()
  ########
  
}
#########################################################################

recentered_p=StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
rm(recentered_final)
gc()

frag <- Fragments(combined@assays$peaksinters)
cat('Making a large count matrix...\n')
matrix.counts <- FeatureMatrix(
  fragments = frag,
  features = recentered_p,
  process_n = 200000,
  sep = c("-","-"),
  cells = colnames(combined)
)
# remove ATAC assay
combined[['peaksinters']] <- NULL

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

cat('Creating chromatin assay...\n')
combined[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                     annotation = annotations,
                                                     genome = 'hg38',
                                                     fragments = frag)
rm(matrix.counts, frag, recentered_p)
gc()

DefaultAssay(combined)<-'X500peaksMACS2'

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_new_peaks_not_normalized_",add_filename,".rds"))

#add some more QC stuff
peak.data <- GetAssayData(object = combined, assay = 'X500peaksMACS2', slot = "counts")
total_fragments_cell <- combined$atac_fragments
peak.counts <- colSums(x = peak.data)
frip <- peak.counts *100 / total_fragments_cell
combined <- AddMetaData(object = combined, metadata = frip, col.name = 'pct_read_in_peaks_500MACS2')
combined <- AddMetaData(object = combined, metadata = peak.counts, col.name = 'peak_region_fragments_500MACS2')

cat('saving the object with updated metadata...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_new_peaks_not_normalized_",add_filename,".rds"))


#run normalization
cat('Normalizing...\n')
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

options(future.globals.maxSize= 891289600)
combined <- FindClusters( 
  object = combined,
  algorithm = 3,
  resolution = 1,
  verbose = FALSE
)

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,".rds"))

n <- length(unique(combined$dataset))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cat('plotting...\n')
p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1, cols = col_vector) + 
  ggplot2::ggtitle("Combined snATAC samples.id")

p2 <- DimPlot(combined, pt.size = 0.1,label=T)
pdf(paste0(length(samples.id),"_snATAC_Merged_", add_filename, ".pdf"),height=6,width=16, useDingbats = F)
p1+p2
dev.off()

write.table(combined@meta.data, paste0(length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,"_metaData.txt"),sep="\t",quote=FALSE, row.names = T)
write.table(samples.id,paste0("Samples_snATAC_Merged_",add_filename,".txt"),sep="\t",quote=FALSE)
