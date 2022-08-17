



all_peaks.split <-  split(all_peaks, by = 'Sample', keep.by = T)

all_peaks.split <- lapply(all_peaks.split, function(x) {
  total.score.per.mil <- sum(x$neg_log10pvalue_summit)/1000000 # this is scaling factor for MACS2 score
  x$score.norm_pval <- x$neg_log10pvalue_summit / total.score.per.mil
  total.score.per.mil <- sum(x$score)/1000000 # this is scaling factor for MACS2 score
  x$score.norm_score <- x$score / total.score.per.mil
  total.score.per.mil <- sum(x$neg_log10qvalue_summit)/1000000 # this is scaling factor for MACS2 score
  x$score.norm_qval <- x$neg_log10qvalue_summit / total.score.per.mil
  x$quantile_qval <- trunc(rank(x$neg_log10qvalue_summit))/length(x$neg_log10qvalue_summit)
  return(x)
})

all_peaks <- rbindlist(all_peaks.split)
all_peaks.split <-  split(all_peaks, by = 'Cancer', keep.by = T)
all_peaks.split <- lapply(all_peaks.split, function(x) {
  x <- x[order(score.norm_pval, decreasing = T),]
  return(x)
})


add_filename <- 'test_qval'
recentered_final.cancer <- foreach(ap=all_peaks.split, cancer.type=names(all_peaks.split))  %do% {
  ap <- ap[order(score.norm_qval, decreasing = T), ] # order peaks by normalized scores, this is essential for filtering out overlapping peaks
  x <- iterative_removal(ap, cancer.type)
  fwrite(x,paste0(cancer.type, '_recentered_final.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  return(x)
}

library(BSgenome.Hsapiens.UCSC.hg38)
recentered_final.cancer.filtered <- lapply(recentered_final.cancer, function(x) {
  x <- x[seqnames!='chrY',]
  filtered.x <- filter_N_peaks(x)
  cancer.type <- x$Cancer[1]
  fwrite(filtered.x, paste0(cancer.type, '_recentered_final.filtered.',add_filename,'.tsv'),
         sep='\t', row.names=FALSE)
  return(filtered.x)
})

iterative_removal <- function(all_peaks.f, cancer.type) {
  print(cancer.type)
  
  #all_peaks.f <- all_peaks.f[order(score.norm, decreasing = T), ]
  #just load existing peaks if any
  # if (file.exists(paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'))) {
  #   recentered_final.f <- fread(paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'))
  # } else {
    
    recentered_p=StringToGRanges(regions = all_peaks.f$new_peak, sep = c("-", "-"))
    
    cat(paste0('finding overlapping peaks in ',cancer.type,'\n'))
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
  #}
  return(recentered_final.f)
}


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

all_peaks.brca.g <- all_peaks.brca.g[order(score.norm_qval, decreasing = T), ]
recentered_p=StringToGRanges(regions = all_peaks.brca.g$new_peak, sep = c("-", "-"))
overlapping=as.data.table(x = findOverlaps(query = recentered_p, 
                                           subject = recentered_p))

overlapping=overlapping[queryHits!=subjectHits,]
overlapping.peak.number <- unique(x = overlapping$queryHits)
all_peaks.brca.g[overlapping.peak.number,]
recentered_non_overlapping=all_peaks.brca.g[-overlapping.peak.number,]
fwrite(recentered_non_overlapping,paste0(cancer.type, '_recentered_nonOverlapping.',add_filename,'.tsv'),
              sep='\t',row.names=FALSE)
tmp <- data.table(chr = all_peaks.brca.g$seqnames[overlapping.peak.number], 
                  num = overlapping.peak.number)
overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) 
overlapping.peak.number.split.chr17 <- overlapping.peak.number.split$chr17
chr = overlapping.peak.number.split.chr17$chr[1]
running.vector <- overlapping.peak.number.split.chr17$num
peaks.to.trash <- NULL
peaks.to.keep <- NULL
while (length(running.vector) != 0) {
  n <- running.vector[1] # this is the first and the best peak since peaks are sorted by scores
  cat(paste('Keep this one:', n, '\n'))
  neighbor.peaks.num.discard <- overlapping[queryHits==n, subjectHits] #find positions of other peaks overlapping with the first one 
  cat(paste('Discard these:', paste(neighbor.peaks.num.discard, collapse = ','), '\n'))
  running.vector <- setdiff(running.vector, neighbor.peaks.num.discard) # remove them from the list of peaks
  running.vector <- setdiff(running.vector, n)
  peaks.to.keep <- c(peaks.to.keep, n) # add this peak to the keeping list
  peaks.to.trash <- unique(c(peaks.to.trash, neighbor.peaks.num.discard))
}
all_peaks.brca.g.keep<- all_peaks.brca.g[peaks.to.keep,]

all_peaks.brca.g[neighbor.peaks.num.discard,]
all_peaks.brca.g[c(n,neighbor.peaks.num.discard),]
length(unique(all_peaks.brca.g$new_peak[peaks.to.keep]))

best_in_overlapping_cancer <- all_peaks.brca.g[peaks.to.keep,]

all_peaks.brca.g[,row_num:=rownames(all_peaks.brca.g)]

p1 <- all_peaks.brca.g
olap1 <- overlapping
pairs=cbind(p1[olap1$queryHits,c('seqnames', 'start', 'end', 'score.norm_qval')],olap1$queryHits,p1[olap1$subjectHits,c('seqnames', 'start', 'end', 'score.norm_qval')],olap1$subjectHits)

colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')
#pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs=pairs[order(-pairs$score_1),]
pairs_all=pairs
p1$row=as.numeric(rownames(p1))
p1=as.data.table(p1)
setkey(p1,row)

pairs=pairs_all[pairs_all$chr_1=='chr17',]
pairs=pairs[,c(4,5,9,10)]
key_list=unique(pairs$row_1)
#pairs=as.data.table(pairs)
setkey(pairs,row_1)
all_st_chr=NULL
for (i in 1:length(key_list)){
  if (length(key_list)>0){
    p_add=p1[key_list[1]]
    all_st_chr=rbind(all_st_chr,p_add)
    
    p_del=pairs[.(key_list[1])]
    del_row_1=c(p_del$row_1,p_del$row_2)
    key_list=key_list[!(key_list %in% del_row_1)]
  }
  #    print(paste(chr,length(key_list),sep=' '))
}


pairs.row_1.less <- pairs[pairs$score_1<pairs$score_2,]
peaks.more.me <- setdiff(peaks.to.keep, key_list.no.row_1.less)
fwrite (all_peaks.brca.g[peaks.more.me,], 'BRCA_peaks.found.by.me_but_not_Nadya.tsv', sep = '\t')

wrong.p <- StringToGRanges(all_st_chr_wrong$new_peak)
overlapping[queryHits %in% peaks.more.me | subjectHits %in% peaks.more.me,]
over.test <- overlapping[queryHits == 255444 | subjectHits == 255444,]

all_peaks.test <- all_peaks.brca.g[unique(c(over.test[['queryHits']], over.test[['subjectHits']])), ]


overlapping[queryHits == 242172 | subjectHits == 242172,]













