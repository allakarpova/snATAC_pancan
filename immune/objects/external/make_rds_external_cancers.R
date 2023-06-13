suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

#suppressMessages(library(EnsDb.Hsapiens.v75))
suppressMessages(library(EnsDb.Hsapiens.v86))

suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
#suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))


#####################################
####### FUNCTIONS ##################
####################################
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


add_blacklist <- function(obj, tb) {
  frag.gr <-  GRanges(seqnames = tb$seqnames,ranges= IRanges(start = tb$start, end = tb$end))
  blacklist.overlaps <- countOverlaps(frag.gr, blacklist_hg38_unified)
  blacklist.fragments <- data.frame(frag.gr[blacklist.overlaps>0,])
  blacklist.fragments <- paste(blacklist.fragments$seqnames, blacklist.fragments$start, blacklist.fragments$end, sep = '-')
  tb$fragment <- paste(tb$seqnames,
                       tb$start,
                       tb$end, sep = '-')
  bl.counts <- tb %>% 
    dplyr::filter(fragment %in% blacklist.fragments & barcode %in% colnames(obj)) %>% 
    group_by(barcode) %>% 
    tally(name = 'blacklist_region_fragments')
  bl.counts <- bl.counts %>% data.frame(row.names=1)
  obj <- AddMetaData(obj, bl.counts)
  obj$blacklist_region_fragments[is.na(obj$blacklist_region_fragments)] <- 0
  return(obj)
}

option_list = list(
  make_option(c("-c", "--cancer"),
              type="character",
              default=NULL,
              help = "cancer type",
              metavar="character"),
  make_option(c("-d","--data"),
              type="character",
              default=NULL,
              help = "path to Cellranger-arc data folder (e.g. cellranger output's raw matrices folder)",
              metavar="character"),
  make_option(c("-m","--macs2_path"),
              type="character",
              default=NULL,
              help = "path to installed MACS2",
              metavar="character"),
  make_option(c("-o","--output_folder"),
              type="character",
              default=NULL,
              help = "output folder where a sample subfolder will be created",
              metavar="character")

) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

cancer=opt$cancer
input=opt$d
#path.to.chrom.size <- opt$chrom_size
outputpath <- opt$output_folder


dir.create(outputpath, showWarnings = F, recursive = T)
setwd(outputpath)

##### START

###### load google sheet and extract samples from there ########
gs4_deauth()
samples.included <- read_sheet("https://docs.google.com/spreadsheets/d/1zOzl4iLMQ0G0wJo0GDiOzMQF94ApYyG0uVP18a66cc0/edit#gid=0", sheet = 1, trim_ws = T)
samples.included$Keep <- samples.included$`Include in immune` %>% unlist()
samples.included <- samples.included %>% 
  dplyr::filter(Keep == 'TRUE') %>% 
  dplyr::filter(`Disease Type`==cancer)
samples.included <- samples.included %>% pull(`Sample ID`)

frag.files <- list.files(path = input, pattern = '*fragments.tsv.gz$')
samples.found <- str_remove(string = frag.files, pattern = '_fragments.tsv.gz')
samples.needed <- intersect(samples.found, samples.included)
frag.files <- list.files(path = input, pattern = '*fragments.tsv.gz$', full.names = T)
frag.files <- frag.files[samples.found %in% samples.included]

print(frag.files)
print(samples.needed)

#chromsize <- fread('/diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg19.chrom.sizes.txt')
chromsize38 <- fread('/diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt')

cat('Making tiles\n')
gr.chr <-  GRanges(seqnames = chromsize38$V1,ranges= IRanges(start = rep(0, 455), end = chromsize38$V2))
tiles <- GenomicRanges::tile(gr.chr,width = 10^6)
stand.chrom <- c(paste0('chr', 1:22), 'chrX')
tiles <- unlist(tiles)
tiles <- tiles[seqnames(tiles) %in% stand.chrom,]

annotations.hg38 <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
# ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations.hg38),sep=""), pattern="chrMT", replacement="chrM")
# seqlevels(annotations.hg38) <- ucsc.levels


walk2(frag.files, samples.needed, function(cur.fr, sample) {
  tryCatch( {
    setwd(outputpath)
    dir.create(sample)
    setwd(sample)
    print(sample)
    #cur.fr <- frag.files[[1]]
    cat('Read in fragment files \n')
    frag.tb <- fread(cur.fr, col.names = c('seqnames', 'start', 'end', 'barcode', 'N'))
    total.frag <- CountFragments(cur.fr)
    atac.cells <- total.frag[total.frag$frequency_count > 1000, "CB"]
    
    #frag.tb %>% dplyr::filter (barcode =='AGCTATGTCGCCCTTA-1') %>% group_by(barcode) %>% summarize(N= sum(N))
    #frag.tb <- fread(cur.fr)
    
    fragments <- CreateFragmentObject(
      path = cur.fr,
      cells = atac.cells, 
      validate.fragments = TRUE
    )
    
    plan("multicore", workers = 10)
    options(future.globals.maxSize = 100 * 1024^3) # for 100 Gb RAM
    
    matrix.counts <- FeatureMatrix(
      fragments = fragments,
      features = tiles,
      sep = c("-","-"),
      process_n = 2000,
      cells = atac.cells
    )
    
    chromatinassay <- CreateChromatinAssay(counts = matrix.counts, 
                                           #genome = "hg38", 
                                           annotation = annotations.hg38,
                                           fragments = fragments)
    object <- CreateSeuratObject(counts = chromatinassay, assay = 'ATAC')
    
    rm(matrix.counts)
    gc(verbose = F)
    
    # compute TSS enrichment score per cell
    object <- TSSEnrichment(object = object, fast = FALSE)
    # inspecting TSS-enrichment scores
    object$high.tss <- ifelse(object$TSS.enrichment > 2, 'High', 'Low')
    tss_plot=TSSPlot(object, group.by = 'high.tss') + NoLegend()
    
    # add blacklist region fragments
    object <- add_blacklist(object, frag.tb)
    
    rm(frag.tb)
    gc(verbose = F)
    
    DefaultAssay(object) <- 'ATAC'
    
    cat('callling peaks\n')
    peaks <- CallPeaks(
      object = object,
      macs2.path=opt$macs2_path
    )
   
    # 
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    
    all_peaks=as.data.table(peaks)
    fwrite(all_peaks,paste0('MACS2_peaks.',sample,'.tsv'),sep='\t')
    
    # recenter peaks
    all_peaks[,peak_center:=start+relative_summit_position]
    all_peaks[,recentered_start:=peak_center-250]
    all_peaks[,recentered_end:=peak_center+250]
    all_peaks[,length:=recentered_end-recentered_start+1]
    all_peaks[,new_peak:=paste0(seqnames,"-", recentered_start, '-',recentered_end)]
    
    ####Now check that new start and end don't go beyond the chromosome boundaries
    
    colnames(chromsize38)=c('seqnames','chr_length')
    all_peaks=merge(all_peaks,chromsize38,all.x=TRUE)
    all_peaks=all_peaks[recentered_end<=chr_length && recentered_start>=0,]
    
    ### do iterative removal of overlapping peaks
    all_peaks <- all_peaks[order(neg_log10qvalue_summit, decreasing = T), ]
    recentered_final <- iterative_removal(all_peaks)
    fwrite(recentered_final,paste0('recentered_final.filtered',sample,'.tsv'),sep='\t')
    
    recentered_p=StringToGRanges(recentered_final$new_peak)
    cat('making counts\n')
    
    
    
    matrix.counts <- FeatureMatrix(
      fragments = fragments,
      features = recentered_p,
      sep = c("-","-"),
      cells = colnames(object)
    )
    
    object[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                       annotation = annotations.hg38,
                                                       #genome = 'hg38',
                                                       fragments = cur.fr, 
                                                       min.features = -1)
    
    DefaultAssay(object)<-'X500peaksMACS2'
    
    # remove ATAC assay
    object[['ATAC']] <- NULL
    
    rm(matrix.counts)
    gc(verbose = F)
    
    #add some more QC stuff
    peak.counts <- colSums(x = GetAssayData(object = object, assay = 'X500peaksMACS2', slot = "counts"))
    rownames(total.frag) <- total.frag$CB
    total_fragments_cell <- total.frag[atac.cells,'frequency_count']
    object <- AddMetaData(object, total_fragments_cell, col.name = 'passed_filters')
    frip <- peak.counts *100 / total_fragments_cell
    object <- AddMetaData(object = object, metadata = frip, col.name = 'frip_500MACS2')
    object <- AddMetaData(object = object, metadata = peak.counts, col.name = 'peak_RF_500MACS2')
    object$blacklist_ratio <- object$blacklist_region_fragments / object$peak_RF_500MACS2
    
    # compute nucleosome signal score per cell
    DefaultAssay(object = object) <- 'X500peaksMACS2'
    object <- NucleosomeSignal(object = object)
    
    # inspecting fragment length periodicity
    object$nucleosome_group <- ifelse(object$nucleosome_signal > 5, paste('NS >', 5),  paste('NS <', 5))
    fragment_period_plot=FragmentHistogram(object = object, group.by = 'nucleosome_group')	
    
    QCplot <- VlnPlot(object = object, 
                      features = c('frip_500MACS2', 'peak_RF_500MACS2','TSS.enrichment','nucleosome_signal', 'blacklist_ratio'), 
                      ncol =5)
    
    pdf(paste(sample,"_0_ATAC_QC.pdf",sep=""),height=9,width=12)
    print(tss_plot)
    print(fragment_period_plot)
    print(QCplot)
    dev.off()
    
    
    #remove cells that are outliers for these QC metrics
    cat('subsetting cool cells\n')
    object <- subset(
      x = object,
      # ATAC QC
      subset = peak_RF_500MACS2 > 1000 &
        peak_RF_500MACS2 < 20000 &
        frip_500MACS2 > 15 &
        nucleosome_signal < 5 &
        TSS.enrichment > 2 &
        blacklist_ratio < 0.05
    )
    
    #### after QC
    #pdf(paste0("After_QC_in_sample_",sample, "_ATAC_MACS2.pdf"), width=25, height=9)
    VlnPlot(object = object, 
            features = c('frip_500MACS2', 'peak_RF_500MACS2',  'TSS.enrichment', 'nucleosome_signal', "blacklist_ratio"), 
            ncol = 5)
    ggsave(paste0("After_QC_in_sample_",sample, "_ATAC_MACS2.pdf"), width=25, height=9)
    
    #run normalization
    cat('Normalizing...\n')
    object <- object %>% RunTFIDF() %>%
      FindTopFeatures(min.cutoff = 20) %>%
      RunSVD(reduction.key = 'LSI_',
             reduction.name = 'lsi',
             irlba.work = 400) %>% 
      RunUMAP(dims = 2:50, reduction = 'lsi',
              reduction.name = "umap.atac", 
              reduction.key = "atacUMAP_") %>%
      FindNeighbors(
        reduction = 'lsi',
        dims = 2:50) %>% 
      FindClusters(
        algorithm = 3,
        resolution = 0.8,
        verbose = FALSE
      )
    
    cat('saving the object...\n')
    saveRDS(object, paste0(sample, "_processed.rds"))
    
    p2 <- DimPlot(object, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
    pdf(paste0(sample,"_2_Dimplots.pdf"),height=6,width=7)
    print(p2)
    dev.off()
  },
  error=function(cond) {
    message(paste("sample failed", sample))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  })
})




