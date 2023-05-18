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
#library(SeuratWrappers)

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
    FindTopFeatures(min.cutoff = 500) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 ) %>% 
    FindNeighbors(
      reduction = 'lsi',
      dims = 2:dims ) %>% 
    FindClusters(
      algorithm = 2,
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
              default="boop", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
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
meta.path <- opt$metadata.file
assay.towork <- opt$assay
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')


#out_path <- '/diskmnt/Projects/snATAC_analysis/immune/obj/v1.0'
#add_filename <- '100_sample_obj.v1.0_old_macs2'
dir.create(out_path, showWarnings = F)
setwd(out_path)

if(!file.exists(paste0('PanImmune_merged_object_100K_random_peaks_normalized_', add_filename, '.rds'))) {
  
  if(!file.exists(paste0('PanImmune_merged_object_100K_random_peaks_', add_filename, '.rds'))) {
    
  cancers <- c('BRCA', 'ccRCC', 'GBM', 'CRC', 'HNSCC', 'MM', 'CESC', 'OV', 'UCEC', "PDAC", "SKCM", 'PBMC')
  cat('creating object on old peaks \n')
  paths <- map(cancers, function(c) {
    p <- list.files(path = input.path, 
                    full.names = T, 
                    pattern = paste0(c,'.*rds'), all.files = F, recursive = T)
    print(length(p))
    return(p)
  })
  paths <- unlist(paths)
  print(paths)
  
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
  } 
  else {
    panc.my <- readRDS(paste0('PanImmune_merged_object_100K_random_peaks_', add_filename, '.rds'))
  }
  cat ('Normalizing 100k random peak object\n')
  panc.my <- runAllNormalization (panc.my, dims = 30)
  ################
  
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
  
  saveRDS(panc.my, paste0('PanImmune_merged_object_100K_random_peaks_normalized_', add_filename, '.rds'))
} else {
  panc.my <- readRDS(paste0('PanImmune_merged_object_100K_random_peaks_normalized_', add_filename, '.rds'))
}

####################################
##### Integration with seurat #######
####################################

panc.my$Data.source <- ifelse(panc.my$Cancer == 'PBMC', '10x', 'DingLab')
panc.my$Batches <- case_when(panc.my$Cancer %in% c('PBMC') ~ paste(panc.my$Cancer, panc.my$data.type, sep = '__'),
                             panc.my$Cancer %in% c('MM') ~ panc.my$Cancer,
                             TRUE ~ panc.my$Chemistry)


atac.split <- SplitObject(panc.my, split.by = 'Batches')

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
  anchor.features = rownames(panc.my),
  reduction = "rlsi",
  dims = 2:50
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = panc.my[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)


# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, 
                      reduction = "integrated_lsi", 
                      dims = 2:50)

integrated  <-  integrated %>% 
  FindNeighbors(
    reduction = 'integrated_lsi',
    dims = 2:50
  ) %>% 
  FindClusters( 
    algorithm = 2,
    resolution = 2,
    verbose = FALSE
  )


print(integrated@reductions)
saveRDS(integrated, paste0('PanImmune_seurat_integrated_object_100K_random_peaks_', add_filename, '_chemistry_data.source.rds'))

fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), 
       paste0('PanImmune_seurat_integrated_object_100K_random_', add_filename, '_chemistry_data.source_matadata.tsv'),
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








