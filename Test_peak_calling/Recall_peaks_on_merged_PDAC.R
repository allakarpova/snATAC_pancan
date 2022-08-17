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
library(BSgenome.Hsapiens.UCSC.hg38)

library(googlesheets4)
library(stringr)
suppressMessages(library(doParallel))


iterative_removal <- function(all_peaks.f, cancer.type) {
  print(cancer.type)
  all_peaks.f <- all_peaks.f[order(score.norm, decreasing = T), ]
  recentered_p=StringToGRanges(all_peaks.f$new_peak, sep = c("-", "-"))
  cat(paste0('finding overlapping peaks in ',cancer.type,'\n'))
  overlapping=as.data.table(findOverlaps(recentered_p,recentered_p)) # find which peaks overlap
  print(dim(overlapping))
  overlapping=overlapping[queryHits!=subjectHits,]
  overlapping.peak.number <- unique(overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks.f table
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
      iterative_removal_core (peak.numbers, overlapping.f = overlapping)
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

getFeatureMatrix <- function (obj, peaks, pro_n) {
  frag <- Fragments(obj@assays$X500peaksMACS2)
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


combined <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_by_cancer_upd/PDAC_12_snATAC_Merged_new_peaks_normalized_v4_samples.rds')
table(combined$seurat_clusters)

write.table(combined@meta.data, paste0("139_snATAC_Merged_new_peaks_normalized_v4_samples_metaData2.txt"),sep="\t",quote=FALSE, row.names = T)


peaks <- CallPeaks(
  object = combined,
  group.by = "seurat_clusters",
  macs2.path='/diskmnt/Projects/Users/allakarpova/Tools/anaconda3/envs/signac/bin/macs2', 
  combine.peaks = F # important!!!
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- lapply(peaks, keepStandardChromosomes, pruning.mode = "coarse")
peaks <- lapply(peaks, subsetByOverlaps, ranges = blacklist_hg38_unified, invert = TRUE)
p <- lapply(peaks, as.data.table)
p <- lapply(p, function(x) {
  total.score.per.mil <- sum(x$neg_log10pvalue_summit)/1000000 # this is scaling factor for MACS2 score
  x$score.norm <- x$neg_log10pvalue_summit / total.score.per.mil
  return(x)
})
p <- rbindlist(p)

fwrite(p,paste('/diskmnt/Projects/snATAC_analysis/test_peak_recalling/','MACS2_peaks._by_cluster','PDAC','.tsv',sep=''),sep='\t')

p[,peak_center:= start + relative_summit_position]
p[,recentered_start:=peak_center-250]
p[,recentered_end:=peak_center+250]

####Now check that new start and end don't go beyond the chromosome boundaries
chr_size=read.table(path.to.chrom.size,sep='\t',header=FALSE)
colnames(chr_size)=c('seqnames','chr_length')

p1=merge(p,chr_size,all.x=TRUE)

p1 = p1[recentered_end <= chr_length & recentered_start >= 0,]
p1[,length:=recentered_end - recentered_start + 1]
p1[,new_peak:=paste(seqnames,recentered_start,recentered_end,sep='-')]


recentered_final <- iterative_removal(p1, 'PDAC_test')
recentered_final <- p1[new_peak %in% recentered_final$new_peak,]
recentered_final <- filter_N_peaks(recentered_final)

recentered_final[,.N,by = 'ident']
recentered_final[,.N,by = 'new_peak'][order(N),]

fwrite(recentered_final, paste('/diskmnt/Projects/snATAC_analysis/test_peak_recalling/','MACS2_peaks.recentered_final.filtered_by_cluster','PDAC','.tsv',sep=''))



recentered_p=StringToGRanges(unique(recentered_final$new_peak), sep = c("-", "-"))
plan("multicore", workers = 30)
options(future.globals.maxSize = 500 * 1024^3) # for 500 Gb RAM

peak.number <- length(unique(recentered_final$new_peak))
n.peaks <- round(peak.number/30)

matrix.counts <- getFeatureMatrix(combined, recentered_p, n.peaks)
frag <- Fragments(combined@assays$X500peaksMACS2)

# extract gene annotations from EnsDb - on congaree fails with database malfunciton WTH??? - didnt fail on lupe
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_by_cancer_common_peaks/Annotations.EnsDb.Hsapiens.v86.rds')

cat('Creating chromatin assay...\n')
combined[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                     annotation = annotations,
                                                     genome = 'hg38',
                                                     fragments = frag)
# remove ATAC assay
DefaultAssay(combined)<-'X500peaksMACS2'
saveRDS(combined, paste0('/diskmnt/Projects/snATAC_analysis/test_peak_recalling/',"PDAC_snATAC_Merged_new_peaks_not_normalized.rds"))

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

cat('saving the object with updated metadata...\n')
saveRDS(combined, paste0('/diskmnt/Projects/snATAC_analysis/test_peak_recalling/',"PDAC_snATAC_Merged_new_peaks_normalized.rds"))

metadata <- fread(paste0('/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v4.0/PDAC_12_samples_metadata_data_freeze_v2.0.tsv')) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F) %>%
  dplyr::rename(seurat_clusters_indiv = seurat_clusters)
combined@meta.data <- combined@meta.data[,c("orig.ident",'dataset', 'nCount_X500peaksMACS2','nFeature_X500peaksMACS2', 'seurat_clusters')]
combined <- AddMetaData(object = combined, metadata = metadata)

total_fragments_cell <- combined$passed_filters
peak.counts <- colSums(x = GetAssayData(combined, slot='counts'))
frip <- peak.counts *100 / total_fragments_cell
combined <- AddMetaData(object = combined, metadata = frip, col.name = 'pct_read_in_peaks_500MACS2')
combined <- AddMetaData(object = combined, metadata = peak.counts, col.name = 'peak_RF_500MACS2')

setwd('/diskmnt/Projects/snATAC_analysis/test_peak_recalling/')

n <- length(unique(combined$cell_type.harmonized))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


p3 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'cell_type.harmonized', cols = col_vector)
ggsave(paste0("PDAC_merged_recalled_peaks_cell_type.harmonized.pdf"), plot = p3,height=12,width=15, useDingbats = F)


p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
ggsave(paste0("PDAC_merged_recalled_peaks_dataset.pdf"),plot = p1,height=12,width=17, useDingbats = F)

p2 <- DimPlot(combined, pt.size = 0.1,label=T)
ggsave(paste0("PDAC_merged_recalled_peaks_clusters.pdf"), plot = p2, height=12,width=14, useDingbats = F)

p4 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'data.type', cols = 'Paired')
ggsave(paste0("PDAC_merged_recalled_peaks_data.type.pdf"), plot = p4,height=12,width=13, useDingbats = F)


panc <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_by_cancer_upd/PDAC_12_snATAC_Merged_new_peaks_not_normalized_v4_samples.rds')
peaks.called <- rownames(combined)
peaks.filtered <- rownames(panc)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
peakAnno.called <- annotatePeak(StringToGRanges(peaks.called), tssRegion=c(-3000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,annoDb = 'org.Hs.eg.db')
peakAnno.filtered <- annotatePeak(StringToGRanges(peaks.filtered), tssRegion=c(-3000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')


nnzero(panc@assays$X500peaksMACS2@counts)
nnzero(combined@assays$X500peaksMACS2@counts)


