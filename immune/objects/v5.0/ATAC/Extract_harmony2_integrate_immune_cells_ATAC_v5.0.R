# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################
library(future)

plan("multicore", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

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
suppressMessages(library(harmony))
################################

#####################################
####### FUNCTIONS ##################
####################################
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

iterative_removal <- function(all_peaks.f) {
  
  #just load existing peaks if any
  recentered_p=StringToGRanges(regions = all_peaks.f$new_peak, sep = c("-", "-"))
  
  cat(paste0('finding overlapping peaks\n'))
  overlapping=as.data.table(x = findOverlaps(query = recentered_p, 
                                             subject = recentered_p)) # find which peaks overlap
  overlapping=overlapping[queryHits!=subjectHits,]
  overlapping.peak.number <- unique(x = overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks.f table
  recentered_non_overlapping=all_peaks.f[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
  # fwrite(recentered_non_overlapping,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_nonOverlapping.',add_filename,'.tsv'),
  #        sep='\t',row.names=FALSE)
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
    # fwrite(best_in_overlapping_cancer,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_Overlapping.',add_filename,'.tsv'),
    #        sep='\t',row.names=FALSE)
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

DoesFilteredExist <- function(add_filename.f) {
  return (file.exists(paste0('recentered_final.filtered',add_filename.f,'.tsv'))) 
}

DoesRawExist <- function(add_filename.f) {
  return(file.exists(paste0('MACS2_peaks.',add_filename.f,'.tsv')))
}

call_peaks <- function(panc.my.f,add_filename.f, macs2_path.f) {
  if (DoesRawExist(add_filename.f)) {
    all_peaks <- fread(paste0('MACS2_peaks.',add_filename.f,'.tsv'))
    return(all_peaks)
  } 
  peaks <- CallPeaks(
    object = panc.my.f,
    group.by = "seurat_clusters",
    macs2.path=macs2_path.f, 
    combine.peaks = F, # important!!!
    format = "BEDPE"
  )
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- peaks[!sapply(peaks, is.null)]
  peaks <- lapply(peaks, keepStandardChromosomes, pruning.mode = "coarse")
  peaks <- lapply(peaks, subsetByOverlaps, ranges = blacklist_hg38_unified, invert = TRUE)
  all_peaks <- lapply(peaks, as.data.table)
  all_peaks <- lapply(all_peaks, function(x) {
    total.score.per.mil <- sum(x$neg_log10qvalue_summit)/1000000 # this is scaling factor for MACS2 score
    x$score.norm <- x$neg_log10qvalue_summit / total.score.per.mil
    return(x)
  })
  all_peaks <- rbindlist(all_peaks)
  return(all_peaks)
}

call_recenter_filter_peaks <- function(panc.my.f,add_filename.f, macs2_path.f ) {
  if (DoesFilteredExist(add_filename.f)) {
    recentered_final <- read.table(paste0('recentered_final.filtered',add_filename.f,'.tsv'),sep='\t',
                                   header = T)
    return (recentered_final)
  }
  #call Macs2 peaks and filter
  all_peaks=call_peaks(panc.my.f,add_filename.f, macs2_path.f)
  fwrite(all_peaks,paste0('MACS2_peaks.',add_filename.f,'.tsv'),sep='\t')
  # recenter peaks
  all_peaks[,peak_center:=start+relative_summit_position]
  all_peaks[,recentered_start:=peak_center-250]
  all_peaks[,recentered_end:=peak_center+250]
  all_peaks[,length:=recentered_end-recentered_start+1]
  all_peaks[,new_peak:=paste0(seqnames,"-", recentered_start, '-',recentered_end)]
  
  ####Now check that new start and end don't go beyond the chromosome boundaries
  chr_size=read.table(chrom.size,sep='\t',header=FALSE)
  colnames(chr_size)=c('seqnames','chr_length')
  all_peaks=merge(all_peaks,chr_size,all.x=TRUE)
  all_peaks=all_peaks[recentered_end<=chr_length & recentered_start>=0,]
  
  ### do iterative removal of overlapping peaks
  all_peaks <- all_peaks[order(neg_log10qvalue_summit, decreasing = T), ]
  recentered_final.f <- iterative_removal(all_peaks)
  recentered_final.f <- filter_N_peaks(recentered_final.f)
  recentered_final.f <- recentered_final.f[seqnames!='chrY',]
  return(recentered_final.f)
}


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
              metavar="character"),
  make_option(c("-m","--macs2_path"),
              type="character",
              default=NULL,
              help = "path to installed MACS2",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character"),
  make_option(c("-s","--chrom_size"),
              type="character",
              default='./hg38.chrom.sizes.txt',
              help = "path to hg38.chrom.sizes.txt",
              metavar="character"),
  make_option(c("-t","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-a","--assay"),
              type="character",
              default='pancan',
              help = "X500peaksMACS2 or peaksinters or pancan",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
cell_column <- opt$cell_type_column
macs2_path <- opt$macs2_path
chrom.size <- opt$chrom_size
meta.path <- opt$metadata.file
assay.towork <- opt$assay
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')

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

if(!file.exists(paste0('PanImmune_merged_object_100K_random_peaks_normalized_', add_filename, '.rds'))) {
  
  
  cancers <- c('BRCA', 'ccRCC', 'GBM', 'CRC', 'HNSCC', 'MM', 'CESC', 'OV', 'UCEC', "PDAC", "SKCM", 'PBMC')
  cat('creating object on old peaks \n')
  paths <- map(cancers, function(c) {
    p <- list.files(path = input.path, 
                    full.names = T, 
                    pattern = paste0(c,'.*rds'), all.files = F, recursive = T)
    print(length(p))
    return(p)
  })
  paths
  
  my.metadata <- fread(meta.path, data.table = F) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F)
  
  cat('opening object...\n')
  atac <- map(paths, function(p) {
    obj=readRDS(p)
    DefaultAssay(obj) <- assay.towork
    if (!file.exists(Fragments(obj)[[1]]@path)) stop("Urgh, this sample object can't locate fragments file")
    obj<- DietSeurat(obj, assay = assay.towork)
    obj <- AddMetaData(obj, my.metadata)
    ct <- obj@meta.data %>% pull(cell_column) %>% unique %>% sort
    print(ct)
    stopifnot(!is.na(ct))
    obj$is_immune <- grepl('Plasma|B-cells|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune|Gran|Mono|MAIT', as.character(unlist(obj[[cell_column]])))
    cat('subsetting\n')
    obj.my <- subset(x = obj, subset = is_immune)
    print(dim(obj.my))
    return(obj.my)
  })
  cat('done\n')
  
  cat ('Reducing peaks\n')
  combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  combined.peaks
  combined.peaks = combined.peaks
  combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
  combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  #peaks.use <- combined.peaks
  peaks.use=sample(combined.peaks, size = 100000, replace = FALSE)
  
  
  registerDoParallel(cores=12)
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
  
  registerDoParallel(cores=12)
  cat ('creating peaksinters and removing useless assays\n')
  atac <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
    obj[['peaksinters']] <- CreateChromatinAssay(counts = co,
                                                 fragments=Fragments(obj[[assay.towork]]), 
                                                 min.cells = -1, min.features = -1)
    #obj$dataset=samples.id[i]
    DefaultAssay(obj)<-'peaksinters'
    ###remove another assay
    obj[[assay.towork]]<-NULL
    return(obj)
  }
  stopImplicitCluster()
  
  
  cat ('Merging on 100k random peaks\n')
  panc.my <- merge(x = atac[[1]], y = atac[-1])
  cat ('done\n')
  DefaultAssay(panc.my) <- 'peaksinters'
  
  #remove individual objects
  rm(atac, matrix.counts)
  gc()
  
  saveRDS(panc.my, paste0('PanImmune_merged_object_100K_random_peaks_', add_filename, '.rds'))
  
  cat ('Normalizing 100k random peak object\n')
  panc.my <- runAllNormalization (panc.my, dims = 30)
  ################
  
  saveRDS(panc.my, paste0('PanImmune_merged_object_100K_random_peaks_normalized_', add_filename, '.rds'))
} else {
  panc.my <- readRDS(paste0('PanImmune_merged_object_100K_random_peaks_normalized_', add_filename, '.rds'))
}

if(!file.exists(paste0('PanImmune_merged_object_new_peaks_', add_filename, '.rds'))) {
  
  print(table(panc.my$seurat_clusters))
  #call peaks
  cat ('Recalling peaks on clusters \n')
  recentered_final <- call_recenter_filter_peaks(panc.my,add_filename, macs2_path )
  fwrite(recentered_final,paste0('recentered_final.filtered',add_filename,'.tsv'),sep='\t')
  
  recentered_p=StringToGRanges(recentered_final$new_peak)
  frag <- Fragments(panc.my@assays[['peaksinters']])
  
  #this will remove fragment objects with just 1 or 0 cells because they fail FeatureMatrix function
  frag.filtered <- frag[do.call( 'c', lapply(frag, function(x) length(x@cells) > 1))]
  cell.count <- table(panc.my$Sample)
  exlude.samples <- names(cell.count[cell.count <= 1])
  
  if (length(exlude.samples)>0) {
    samples <- unique(panc.my$Sample)
    samples.to.keep <- samples[!(samples %in% exlude.samples)]
    #now exclude these samples from the object
    panc.my <- subset(panc.my, subset = Sample %in% samples.to.keep)
  }
  
  peak.number <- length(unique(recentered_final$new_peak))
  n.peaks <- round(peak.number/30)
  
  matrix.counts <- FeatureMatrix(
    fragments = frag.filtered,
    features = recentered_p,
    process_n = n.peaks,
    sep = c("-","-"),
    cells = colnames(panc.my)
  )
  
  panc.my[['ATAC_immune']] <- CreateChromatinAssay(counts = matrix.counts,
                                                   fragments=frag.filtered,
                                                   annotation = annotations)
  DefaultAssay(panc.my)<-'ATAC_immune'
  panc.my <- DietSeurat(panc.my, assays = 'ATAC_immune')
  
  peak.counts <- colSums(GetAssayData(object = panc.my, assay = 'ATAC_immune', slot = "counts"))
  total_fragments_cell <- panc.my$passed_filters
  frip <- peak.counts / total_fragments_cell
  panc.my <- AddMetaData(object = panc.my, metadata = frip, col.name = 'pct_read_in_peaks_ATAC_immune')
  panc.my <- AddMetaData(object = panc.my, metadata = peak.counts, col.name = 'peak_region_fragments_ATAC_immune')
  
  
  #### remove excessive fragment files, artifact after peak calling on cluster level
  cat('removing excessive fragment files\n')
  all.fragment.obj <- Fragments(panc.my)
  all.fragment.obj.cell.count <- map_chr(all.fragment.obj, function(x) length(x@cells))
  all.fragment.obj.upd <- all.fragment.obj[all.fragment.obj.cell.count > 0]
  Fragments(panc.my) <- NULL
  Fragments(panc.my) <- all.fragment.obj.upd
  
  panc.my <- panc.my %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 20) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 )  %>%
    RunUMAP(reduction = "lsi", dims = 2:50)
  
  # save
  saveRDS(panc.my, paste0('PanImmune_merged_object_new_peaks_', add_filename, '.rds'))
} else {
  panc.my <- readRDS(paste0('PanImmune_merged_object_new_peaks_', add_filename, '.rds'))
}
####################################
##### Integration with harmony #######
####################################

panc.my$Data.source <- ifelse(panc.my$Cancer == 'PBMC', '10x', 'DingLab')
panc.my$Batches <- case_when(panc.my$Cancer %in% c('PBMC', 'MM') ~ panc.my$Cancer,
                             TRUE ~ panc.my$Chemistry)

plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024^3)

integrated <-  panc.my %>%
  RunHarmony( group.by.vars = c("Batches"), 
              reduction = 'lsi', 
              assay.use = 'ATAC_immune',  
              project.dim = F, 
              max.iter.harmony = 20) %>%
  FindNeighbors(reduction = "harmony", dims = 2:50) %>%
  FindClusters(verbose = TRUE, resolution = 2) %>%
  RunUMAP(reduction = "harmony", dims = 2:50)

print(integrated@reductions)
saveRDS(integrated, paste0('PanImmune_harmony_integrated_object_new_peaks_', add_filename, '_chemistry_data.source.rds'))

fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_harmony_integrated_object_new_peaks_', add_filename, '_chemistry_data.source_matadata.tsv'),
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


saveRDS(integrated, paste0('PanImmune_harmony_integrated_object_new_peaks_', add_filename, '_chemistry_data.source.rds'))

fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_harmony_integrated_object_new_peaks_', add_filename, '_chemistry_data.source_matadata.tsv'),
       sep='\t', row.names = T)









