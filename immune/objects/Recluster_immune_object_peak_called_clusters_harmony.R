# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################
library(future)

plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
library(harmony)
################################

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
    combine.peaks = F # important!!!
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

############################################

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
cell_column <- opt$cell_type_column
macs2_path <- opt$macs2_path
chrom.size <- opt$chrom_size
meta.path <- opt$metadata.file
assay.towork <- opt$assay
  
#   input.path <-'/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v3.0/snATAC/Merged_objects_PanCan_peaks/137Samples_PanCan_merged_obj/137_snATAC_113K_peaks_diffPCs.motifsAdded.chromvar.20210826.rds.gz'
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


#read in input
cat('opening object...\n')
panc <- readRDS(input.path)
cat('done\n')

# add meta data if provided
if (!is.null(meta.path)) {
  my.metadata <- fread(meta.path, data.table = F) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F)
  panc@meta.data <- panc@meta.data[,c('seurat_clusters', 'nCount_pancan_s')]
  panc <- AddMetaData(panc, metadata = my.metadata)
}
panc$tumor_type <- case_when(panc$Cancer=='MM' ~ 'Liquid',
                             TRUE ~ 'Solid')
panc$is_immune <- grepl('Plasma|B-cells|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune|Kuppfer|Lympho|Neutro|Gran', as.character(unlist(panc[[cell_column]])))

my.metadata <- dplyr::filter(my.metadata, grepl('Plasma|B-cells|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune|Kuppfer|Lympho|Neutro|Gran',.data[[cell_column]]))
dim(panc)
#subsetting
cat('subsetting\n')
panc.my <- subset(x = panc, subset = is_immune)
dim(panc.my)

saveRDS(panc.my, paste0('Immune_object_old_peaks_', add_filename, '.rds'))

rm(panc)
gc()


#### run normalization to get initial clusters ###
########
panc.my <- panc.my %>% 
  RunTFIDF() %>% 
  FindTopFeatures(min.cutoff = 20) %>% 
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400)

panc.my <- panc.my %>%
  RunHarmony( c("data.type", 'tumor_type'), reduction = 'lsi', assay.use = assay.towork,  project.dim = F) %>%
  FindNeighbors(reduction = "harmony", dims = 2:30) %>%
  FindClusters(verbose = TRUE, resolution = 2, algorithm = 3) %>%
  RunUMAP(reduction = "harmony", dims = 2:30)

################
saveRDS(panc.my, paste0('Immune_object_old_peaks_', add_filename, '.rds'))

recentered_final <- call_recenter_filter_peaks(panc.my,add_filename, macs2_path )
fwrite(recentered_final,paste0('recentered_final.filtered',add_filename,'.tsv'),sep='\t')

recentered_p=StringToGRanges(recentered_final$new_peak)
frag <- Fragments(panc.my@assays[[assay.towork]])
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
                                                 fragments=frag.filtered)
DefaultAssay(panc.my)<-'ATAC_immune'
panc.my <- DietSeurat(panc.my, assays = 'ATAC_immune')

peak.counts <- colSums(GetAssayData(object = panc.my, assay = 'ATAC_immune', slot = "counts"))
total_fragments_cell <- panc.my$passed_filters
frip <- peak.counts / total_fragments_cell
panc.my <- AddMetaData(object = panc.my, metadata = frip, col.name = 'pct_read_in_peaks_ATAC_immune')
panc.my <- AddMetaData(object = panc.my, metadata = peak.counts, col.name = 'peak_region_fragments_ATAC_immune')


##########################
##### RECLUSTERING #######
##########################

panc.my <- panc.my %>% 
  RunTFIDF() %>% 
  FindTopFeatures( min.cutoff = 20) %>% 
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400)

panc.my <- panc.my %>%
  RunHarmony( c("data.type", 'tumor_type'), reduction = 'lsi', assay.use = 'ATAC_immune',  project.dim = F) %>%
  FindNeighbors(reduction = "harmony", dims = 2:50) %>%
  FindClusters(verbose = TRUE, resolution = 2, algorithm = 3) %>%
  RunUMAP(reduction = "harmony", dims = 2:50)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(panc.my) <- annotations	


#extract gene coordinates and extend them to include the 2 kb upstream region (as promoter accessibility is often correlated with gene expression)
gene.activities <- GeneActivity(panc.my)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
panc.my[['RNA']] <- CreateAssayObject(counts = gene.activities)
panc.my <- NormalizeData(
  object = panc.my,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(panc.my$nCount_RNA)
)

cat('saving the object...\n')
saveRDS(panc.my, paste0("Reclustered_immune_snATAC_Merged_",add_filename,".rds"))


cat('plotting...\n')
p1 <- DimPlot(panc.my, group.by = 'Piece_ID', pt.size = 0.1) + 
  ggplot2::ggtitle("Combined snATAC samples")
p2 <- DimPlot(panc.my, pt.size = 0.1,label=T)
p3 <- DimPlot(panc.my,group.by = cell_column, pt.size = 0.1,label=T)
pdf(paste0("Reclustered_immune_snATAC_Merged_", add_filename, ".pdf"),height=30,width=25, useDingbats = F)
print(p1/p2/p3)
dev.off()

p <- DimPlot(panc.my, group.by = 'Cancer', pt.size = 0.1)
pdf(paste0("Dimplot_cancer_reclustered_immune_snATAC_Merged_", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

p <- DimPlot(panc.my, group.by = 'Cancer', split.by ='Sample_type',  pt.size = 0.1, ncol = 3, cols = 'Paired')
pdf(paste0("Dimplot_Cancer_split_by_type_reclustered_immune_snATAC_Merged_", add_filename, ".pdf"),height=10,width=36, useDingbats = F)
print(p)
dev.off()


p <- DimPlot(panc.my, group.by = 'final.doublet', pt.size = 0.1, cols = c('black', 'yellow'))
pdf(paste0("Dimplot_final.doublet_reclustered_immune_snATAC_Merged_", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()




fwrite(cbind(Embeddings(panc.my, reduction='umap'), panc.my@meta.data), paste0("Reclustered_immune_snATAC_Merged_",add_filename,"metadata.tsv"),
       sep='\t', row.names = T)









