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
plan("multiprocess", workers = 50)
options(future.globals.maxSize = 250 * 1024^3) # for 250 Gb RAM

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
# if peaks were already taken care of load this file
if(file.exists(paste0(length(samples.id),'_recentered_final.filtered.',add_filename,'.tsv'))) {
  recentered_final <- fread(paste0(length(samples.id),'_recentered_final.filtered.',add_filename,'.tsv'), data.table = F)
} else {
  # overlap peaks
  cat('work on peaks...\n')
  all_peaks=NULL
  for (sample in samples.id){
    peaks=read.table(paste0(input.path,'/', sample,"/recentered_final.filtered",sample,".tsv"),sep='\t',header=TRUE)
    peaks$Sample=sample
    peaks$new_peak = paste(peaks$seqnames, peaks$recentered_start, peaks$recentered_end, sep = '-')
    all_peaks=rbind(all_peaks,peaks)
  }
  
  fwrite(all_peaks, paste0(length(samples.id),"_sample_MACS2_peaks_",add_filename,".tsv"),
         sep='\t',row.names=FALSE)
  
  #recenter peaks
  recentered_p=StringToGRanges(all_peaks$new_peak, sep = c("-", "-"))
  cat('find overlaping peaks...\n')
  olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
  olap1=olap[olap$queryHits!=olap$subjectHits,]
  rm(recentered_p)
  gc()
  
  #select non overlapping peaks
  recentered_non_olap=all_peaks[-olap1$queryHits,]
  pairs_all=cbind(all_peaks[olap1$queryHits,c(1:3,7)],
                  olap1$queryHits,
                  all_peaks[olap1$subjectHits,c(1:3,7)],
                  olap1$subjectHits)
  colnames(pairs_all)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')
  
  pairs_all=pairs_all[pairs_all$score_1>=pairs_all$score_2,]
  pairs_all=pairs_all[order(-pairs_all$score_1),]
  
  #pairs_all <- split
  library(doParallel)
  registerDoParallel(cores=30)
  all_st=NULL
  all_st <- foreach(chr_n=c(1:22,"X","Y")) %dopar% {
    chr=paste("chr",chr_n,sep='')
    pairs=pairs_all[pairs_all$chr_1==chr,]
    pairs=pairs[,c(4,5,9,10)]
    all_st_chr=NULL
    for (i in 1:nrow(pairs)){
      if (nrow(pairs)>0){
        p_del=pairs[pairs$row_1==pairs[1,2],]
        all_st_chr=rbind(all_st_chr,all_peaks[rownames(all_peaks)==pairs[1,2],])
        pairs=pairs[!(pairs$row_1 %in% c(p_del$row_1[1],p_del$row_2)),]
      }
    }
    return(all_st_chr)
  }
  stopImplicitCluster()
  
  all_st_f=NULL
  for (i in 1:24){
    all_st_1=as.data.frame(all_st[[i]])
    all_st_1=all_st_1[!duplicated(all_st_1),]
    all_st_f=rbind(all_st_f,all_st_1)
  }
  
  cat('done...\n')
  recentered_final=rbind(recentered_non_olap,all_st_f, olap)
  
  fwrite(recentered_final,paste0(length(samples.id),'_recentered_final.filtered.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  fwrite(recentered_non_olap,paste0(length(samples.id),'_recentered_nonOverlapping.filtered.',add_filename,'.tsv'),
         sep='\t',row.names=FALSE)
  fwrite(all_st_f,paste0(length(samples.id),'_recentered_Overlapping.filtered.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  #remove big things
  rm(recentered_non_olap, all_st_f, all_st,all_st_1,olap, olap1, pairs_all)
  gc()
}
#########################################################################

recentered_p=StringToGRanges(recentered_final$new_peak, sep = c(":", "-"))
frag <- Fragments(combined@assays$peaksinters)
cat('Making large matrix counts...\n')
matrix.counts <- FeatureMatrix(
  fragments = frag,
  features = recentered_p,
  sep = c("-","-"),
  cells = colnames(combined)
)
# remove ATAC assay
combined[['peaksinters']] <- NULL

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
