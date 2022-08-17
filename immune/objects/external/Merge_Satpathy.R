## merge ATAC samples which can be either regular ATAC seq or combo ATAC seq
#v 4.1 Alla implemented checking that all fragments files exist before starting merging
#v5.5 removed filtering by reproducibility, only iterative removal and removal by chrY and NNNN regions remained
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

############## FUNCTIONS #####################
filter_N_peaks <- function(peak.dt) {
  cancer.type <- peak.dt$Cancer[1]
  gr <- StringToGRanges(peak.dt$new_peak, sep = c("-", "-")) #get GRanges object from peaks
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr) #extract fasta sequence
  names(seq) <- peak.dt$new_peak
  peaks.match.pattern <- vmatchPattern("N", seq) #match peak sequence with N in them
  peaks.withN <- names(peaks.match.pattern)[elementNROWS(peaks.match.pattern)>0] # these are peaks that contain N in their sequence
  toreturn <- peak.dt[! new_peak %in% peaks.withN,]
  # fwrite(toreturn,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.reproducible.filtered.',add_filename,'.tsv'),sep='\t',
  #        row.names=FALSE)
  return(toreturn)
}

filter_reproducible <- function(recentered_final.f, all_peaks.f) {
  cancer.type <- all_peaks.f$Cancer[1]
  print(cancer.type)
  current_peaks <- unique(recentered_final.f$new_peak)
  peaks.tokeep <- NULL
  all_peaks.current <- all_peaks.f[new_peak %in% current_peaks,]
  updated_peaks <- all_peaks.current[score.norm>=5, .N, by='new_peak'][N>=2,][['new_peak']] # this counts peak occurance across samples and selects peaks found in at least 2 samples
  toreturn <- all_peaks.f[new_peak %in% updated_peaks,] 
  # fwrite(toreturn,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.reproducible.',add_filename,'.tsv'),sep='\t',
  #        row.names=FALSE)
  return(toreturn)
}

filter_reproducible_by_overlap <- function(recentered_final.f, all_peaks.f) {
  cancer.type <- all_peaks.f$Cancer[1]
  print(cancer.type)
  
  current_peaks <- unique(recentered_final.f$new_peak)
  
  all_peaks.current <- all_peaks.f[new_peak %in% current_peaks,]
  #these are peaks that we going to keep because they are exactly reproducible in samples:
  peaks.exactly.reproduc <- all_peaks.current[score.norm>=5, .N, by='new_peak'][N>=2,][['new_peak']] # this counts peak occurance across samples and selects peaks found in at least 2 samples
  
  # from these peaks we are trying to rescue some peaks by overlap
  peaks.in.question <- setdiff(current_peaks, peaks.exactly.reproduc)
  
  all_peaks.highscore <- all_peaks.f[!(new_peak %in% peaks.in.question) & score.norm>=5,] # remove peaks in question from the total table so that we dont overlap them with themselves
  peaks.in.question <- setdiff(current_peaks, peaks.exactly.reproduc)
  peaks.in.question.gr <- StringToGRanges(peaks.in.question)
  all_peaks.highscore.gr <- StringToGRanges(all_peaks.highscore$new_peak)
  #overlap peaks in question with filtered all peaks table
  overlapping <- as.data.table(findOverlaps(peaks.in.question.gr, all_peaks.highscore.gr, minoverlap = 251))
  #summarize how many hits each peak in question has, select peaks with >=2 hits
  peaks.rescued.num <- overlapping[,.N, by='queryHits'][N>=2,][['queryHits']]
  peaks.rescued <- peaks.in.question[peaks.rescued.num]
  updated_peaks <- c(peaks.exactly.reproduc, peaks.rescued)
  toreturn <- all_peaks.f[new_peak %in% updated_peaks,] 
  # fwrite(toreturn,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.reproducible.',add_filename,'.tsv'),sep='\t',
  #        row.names=FALSE)
  return(toreturn)
}

iterative_removal <- function(all_peaks.f, cancer.type) {
  print(cancer.type)
  
  #all_peaks.f <- all_peaks.f[order(score.norm, decreasing = T), ]
  #just load existing peaks if any
  if (file.exists(paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'))) {
    recentered_final.f <- fread(paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'))
  } else {
    
    recentered_p=StringToGRanges(regions = all_peaks.f$new_peak, sep = c("-", "-"))
    
    cat(paste0('finding overlapping peaks in ',cancer.type,'\n'))
    overlapping=as.data.table(x = findOverlaps(query = recentered_p, 
                                               subject = recentered_p)) # find which peaks overlap
    overlapping=overlapping[queryHits!=subjectHits,]
    overlapping.peak.number <- unique(x = overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks.f table
    recentered_non_overlapping=all_peaks.f[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
    fwrite(recentered_non_overlapping,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_nonOverlapping.',add_filename,'.tsv'),
           sep='\t',row.names=FALSE)
    if (length(overlapping.peak.number)>0) {
      tmp <- data.table(chr = all_peaks.f$seqnames[overlapping.peak.number], 
                        num = overlapping.peak.number)
      overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) #split peaks by chromosome 
      registerDoParallel(cores=25)
      #this is where iterative removal of peaks is done
      best_in_overlapping_num <- foreach(peak.numbers=overlapping.peak.number.split) %dopar% {
        cat('removing overlapping peaks in each chromosome\n')
        iterative_removal_core (peak.numbers = peak.numbers, overlapping.f = overlapping)
      }
      stopImplicitCluster()
      best_in_overlapping_num <- do.call('c', best_in_overlapping_num) #combine best peak numbers from all chromosomes
      best_in_overlapping_cancer <- all_peaks.f[best_in_overlapping_num,] #extract peaks themselves
      fwrite(best_in_overlapping_cancer,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_Overlapping.',add_filename,'.tsv'),
             sep='\t',row.names=FALSE)
      recentered_final.f=rbindlist(list(recentered_non_overlapping,best_in_overlapping_cancer))
    } else {
      recentered_final.f=recentered_non_overlapping
    }
    final.overlaps <-  recentered_final.f$new_peak %>% 
      unique %>% 
      StringToGRanges %>% 
      countOverlaps
    if (sum(final.overlaps>1)>0) {
      stop("Execution stopped. Overlapping peaks remained")
    }
    # fwrite(recentered_final.f,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'),sep='\t',
    #        row.names=FALSE)
  }
  return(recentered_final.f)
}

# this works like a charm
iterative_removal_core <- function(peak.numbers, overlapping.f) {
  chr = peak.numbers$chr[1]
  running.vector <- peak.numbers$num
  peaks.to.trash <- NULL
  peaks.to.keep <- NULL
  while (length(running.vector) != 0) {
    n <- running.vector[1] # this is the first and the best peak since peaks are sorted by scores
    neighbor.peaks.num.discard <- overlapping.f[queryHits==n, subjectHits] #find positions of other peaks overlapping with the first one 
    running.vector <- setdiff(running.vector, neighbor.peaks.num.discard) # remove them from the list of peaks
    running.vector <- setdiff(running.vector, n)
    peaks.to.keep <- c(peaks.to.keep, n) # add this peak to the keeping list
    peaks.to.trash <- unique(c(peaks.to.trash, neighbor.peaks.num.discard)) # add neighbors to the list of peaks to discard
  }
  cat('done\n')
  return(peaks.to.keep)
}

load_peaks <- function(sample){
  peaks=fread(paste0(input.path,'/', sample,"/recentered_final.filtered",sample,".tsv"))
  peaks$Sample=sample
  peaks$Cancer=str_split_fixed(sample, '_',2)[,1]
  peaks$new_peak = paste(peaks$seqnames, peaks$recentered_start, peaks$recentered_end, sep = '-') #make sure all peaks are written with -
  total.score.per.mil <- sum(peaks$neg_log10qvalue_summit)/1000000 # this is scaling factor for MACS2 score
  peaks$score.norm <- peaks$neg_log10qvalue_summit / total.score.per.mil # normalize peak score in each sample aka score per million
  return(peaks)
}

getFeatureMatrix <- function (obj, peaks, pro_n) {
  frag <- Fragments(obj@assays$peaksinters)
  cat('Making a large count matrix...\n')
  matrix.counts <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    process_n = pro_n,
    sep = c("-","-"),
    cells = colnames(obj)
  )
  return(matrix.counts)
}

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

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 100 * 1024^3) # for 100 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1zOzl4iLMQ0G0wJo0GDiOzMQF94ApYyG0uVP18a66cc0/edit?usp=sharing", sheet = 1, trim_ws = T)
samples$Keep <- samples$`Include in immune` %>% unlist()
#samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
#samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data folder`, Cell_type)

samples.id <- samples$Sample %>% as.character()
samples.type <- samples$`Data Type` %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))

#########################################################################
#if the object has been merged on 5k random peaks, just open it
if (file.exists(paste0(length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))) {
  cat('opening the object...\n')
  combined <- readRDS(paste0(length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))
  #if not, create merged object on 5k random peaks
} else {
  cat('creating first 5k object \n')
  paths <- NULL
  for (i in 1:length(samples.id)){
    print(samples.id[i])
    p <- list.files(path = input.path, full.names = T, pattern = paste0(str_split_fixed(samples.id[i], '_',2)[2],'.*rds'), all.files = T, recursive = T)
    print(length(p))
    paths <- c(paths, p)
  }
  #stop if not all samples have RDS object
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
    obj$Sample = samples.id[i]
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
      fragments = Fragments(obj[[assay.towork]]),
      features = peaks.use,
      sep = c("-","-"),
      cells = colnames(obj)
    ) 
  }
  stopImplicitCluster()
  
  registerDoParallel(cores=10)
  checking.n.cells <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
    return(ncol(obj)==ncol(co))
    #stopifnot(ncol(obj)==ncol(co))
  }
  names(checking.n.cells) <- samples.id
  print(checking.n.cells[!checking.n.cells])
  stopImplicitCluster()
  
  
  registerDoParallel(cores=10)
  cat ('creating peaksinters and removing useless assays\n')
  atac <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
    obj[['peaksinters']] <- CreateChromatinAssay(counts = co,fragments=Fragments(obj@assays$X500peaksMACS2), min.cells = -1, min.features = -1)
    #obj$dataset=samples.id[i]
    DefaultAssay(obj)<-'peaksinters'
    ###remove another assay
    obj[[assay.towork]]<-NULL
    return(obj)
  }
  stopImplicitCluster()
  
  
  ####Merging on old 5k peaks
  cat ('Merging\n')
  combined <- merge(x = atac[[1]], y = atac[-1], add.cell.ids = samples.id)
  DefaultAssay(combined) <- "peaksinters"
  
  #remove individual objects
  rm(atac)
  gc()
  
  cat('saving the 5k object...\n')
  saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))
  
}
#########################################################################

#########################################################################
# if peaks were already taken care of, load this file
if(file.exists(paste0('peaks/', length(samples.id),'_All_recentered_final.filtered.',add_filename,'.tsv'))) {
  cat('loading ready peaks \n')
  recentered_final.filtered <- fread( paste0('peaks/',length(samples.id),'_All_recentered_final.filtered.',add_filename,'.tsv'), data.table = T)
  #if not create peaks, filter blah blah
} else {
  # read in peaks
  dir.create('peaks/')
  cat('work on peaks...\n')
  registerDoParallel(cores=10)
  all_peaks <- foreach (sample=samples.id) %dopar% {
    sample=str_split_fixed(sample, '_', 2)[2]
    load_peaks (sample)
  }
  stopImplicitCluster()
  
  all_peaks <- rbindlist(all_peaks) # peaks from all samples
  all_peaks <- all_peaks[order(score.norm, decreasing = T), ] # order peaks by normalized scores, this is essential for filtering out overlapping peaks
  fwrite(all_peaks, paste0('peaks/',length(samples.id),"_sample_MACS2_peaks_",add_filename,".tsv"),
         sep='\t',row.names=FALSE)
  
  #all_peaks <- split(all_peaks, by = 'Cancer', keep.by = T)
  cat('Iterative removal \n')
  # filter out overlapping peaks via iterative removal 
  all_peaks <- all_peaks[order(score.norm, decreasing = T), ] # order peaks by normalized scores, this is essential for filtering out overlapping peaks
  recentered_final <- iterative_removal(all_peaks, 'all')
  fwrite(recentered_final,paste0('peaks/',length(samples.id),'_All_recentered_final.nonoverlap.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  
  
  cat('Filtering peaks by chrY and NNN \n')
  #filter out peaks in chrY and by N in peak sequence
  library(BSgenome.Hsapiens.UCSC.hg38)
  recentered_final.filtered <- recentered_final[seqnames!='chrY',]
  recentered_final.filtered <-  filter_N_peaks(recentered_final.filtered)
  fwrite(recentered_final.filtered, paste0('peaks/',length(samples.id),'_All_recentered_final.filtered.',add_filename,'.tsv'),
         sep='\t', row.names=FALSE)
  
  rm(recentered_final, all_peaks)
  gc()
  
  
}
#########################################################################
cat('Number of peaks\n')
print(nrow(recentered_final.filtered))
recentered_p=StringToGRanges(unique(recentered_final.filtered$new_peak), sep = c("-", "-"))
#rm(recentered_final.filtered)
#gc()

###some parallelization-solution from the tutorial:
plan("multicore", workers = 30)
options(future.globals.maxSize = 500 * 1024^3) # for 500 Gb RAM


peak.number <- length(unique(recentered_final.filtered$new_peak))
n.peaks <- round(peak.number/30)
matrix.counts <- getFeatureMatrix(combined, recentered_p, n.peaks)

# extract gene annotations from EnsDb - on congaree fails with database malfunciton WTH??? - didnt fail on lupe
if (file.exists('Annotations.EnsDb.Hsapiens.v86.rds')) {
  annotations <- readRDS('Annotations.EnsDb.Hsapiens.v86.rds')
} else {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,  standard.chromosomes = TRUE)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  saveRDS(annotations, ' Annotations.EnsDb.Hsapiens.v86.rds')
}

cat('Creating chromatin assay...\n')
frag <- Fragments(combined@assays$peaksinters)
combined[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                     annotation = annotations,
                                                     #genome = 'hg38',
                                                     fragments = frag)
# remove ATAC assay
DefaultAssay(combined)<-'X500peaksMACS2'
combined[['peaksinters']] <- NULL

rm(matrix.counts, frag, recentered_p)
gc()


cat('saving the object...\n')
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

combined <- FindClusters( 
  object = combined,
  algorithm = 3,
  resolution = 0.8,
  verbose = FALSE
)

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,".rds"))

### take care of metadata here
#add some more QC stuff
#this is metadata specific, redo if needed
#combined <- readRDS('40_snATAC_Merged_new_peaks_normalized_external_Satpathy_BCC_healthy_immune.rds')
colnames(combined@meta.data)[colnames(combined@meta.data)=='dataset'] <- 'Sample'
obj.meta <- combined@meta.data
rownames(samples) <- samples$Sample

samples <- samples[obj.meta$Sample,]
dim(samples)
to.add <- samples[,c('Disease Type', 'cell_type',	'Case ID',	'Piece_ID', 'Data Type', 'Chemistry')]
to.add <- to.add %>% dplyr::rename(Cancer=`Disease Type`,
                                   Case_ID=`Case ID`,
                                   data.type=`Data Type`)
obj.meta <- cbind(obj.meta,to.add)
combined <- AddMetaData(object = combined, metadata = obj.meta)

# this works always
total_fragments_cell <- combined$passed_filters
peak.counts <- colSums(x = GetAssayData(combined, slot = 'counts'))
frip <- peak.counts *100 / total_fragments_cell
combined <- AddMetaData(object = combined, metadata = frip, col.name = 'pct_read_in_peaks_500MACS2')
combined <- AddMetaData(object = combined, metadata = peak.counts, col.name = 'peak_RF_500MACS2')


cat('saving the object with updated metadata...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_new_peaks_normalized2_",add_filename,".rds"))



library(RColorBrewer)
n <- length(unique(combined$Sample))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cat('plotting...\n')
p1 <- DimPlot(combined, group.by = 'Sample', pt.size = 0.1)
ggsave(paste0(length(samples.id),"_snATAC_Merged_sample_", add_filename, ".pdf"),plot = p1,height=12,width=40, useDingbats = F)

p1 <- DimPlot(combined, group.by = 'Piece_ID', pt.size = 0.1,label = T)
ggsave(paste0(length(samples.id),"_snATAC_Merged_piece_id_", add_filename, ".pdf"),plot = p1,height=12,width=21, useDingbats = F)

p1 <- DimPlot(combined, group.by = 'Piece_ID', split.by = 'Cancer', ncol =3, pt.size = 0.1,label = T)
ggsave(paste0(length(samples.id),"_snATAC_Merged_piece_id_split_", add_filename, ".pdf"),plot = p1,height=15,width=50, useDingbats = F,limitsize = FALSE)

p2 <- DimPlot(combined, pt.size = 0.1,label=T)
ggsave(paste0(length(samples.id),"_snATAC_Merged_cluster0.8_", add_filename, ".pdf"), plot = p2, height=12,width=14, useDingbats = F)


n <- length(unique(combined$cell_type))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p3 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'cell_type', cols = col_vector)
ggsave(paste0( length(samples.id),"_snATAC_Merged_cell_type.harm_", add_filename, ".pdf"), plot = p3,height=12,width=15, useDingbats = F)

p4 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'data.type', cols = 'Paired')
ggsave(paste0(length(samples.id),"_snATAC_Merged_data.type_", add_filename, ".pdf"), plot = p4,height=12,width=13, useDingbats = F)

p6 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'Cancer', cols = 'Paired')
ggsave(paste0(length(samples.id),"_snATAC_Merged_cancer_", add_filename, ".pdf"), plot = p6,height=12,width=13, useDingbats = F)

p1 <- DimPlot(combined, group.by = 'cell_type', split.by = 'Cancer', ncol =3, pt.size = 0.1,label = T, cols = col_vector)
ggsave(paste0(length(samples.id),"_snATAC_Merged_cell_type_split_", add_filename, ".pdf"),plot = p1,height=15,width=50, useDingbats = F,limitsize = FALSE)


write.table(cbind(Embeddings(combined,reduction = 'umap'), combined@meta.data), 
            paste0(length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,"_metaData.txt"),sep="\t",quote=FALSE, row.names = T)
write.table(samples.id,paste0("Samples_snATAC_Merged_",add_filename,".txt"),sep="\t",quote=FALSE)



