

#data table version of code
all_peaks <- fread('139_sample_MACS2_peaks_v4_samples.tsv')
recentered_p=StringToGRanges(all_peaks$new_peak, sep = c("-", "-"))

overlapping=as.data.table(findOverlaps(recentered_p,recentered_p))
overlapping=overlapping[queryHits!=subjectHits,]
fwrite(overlapping, 'tmp_files/139_overlapping_MACS2_peak_number.tsv', sep = '\t')
overlapping.peak.number <- unique(overlapping$queryHits)


overlap.counts <- countOverlaps(recentered_p, recentered_p)
length(overlap.counts[overlap.counts>1])
recentered_non_overlapping=all_peaks[overlap.counts==1,]
fwrite(recentered_non_overlapping,paste0(length(samples.id),'_recentered_nonOverlapping.filtered.',add_filename,'.tsv'),
       sep='\t',row.names=FALSE)
rm(overlap.counts, recentered_p)
gc()

#split overlapping peaks by chromosome and process 500k peaks at a time
tmp <- data.table(chr = all_peaks$seqnames[overlapping.peak.number], num = overlapping.peak.number)
overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) # overlapping peak number split by chromosome

scores <- all_peaks$score # save scores for each peak
best_in_overlapping <- lapply (X = overlapping.peak.number.split, 
                               FUN = function(overlapping.peak.number) { # apply function on each chromosome peaks, overlapping.peak.number contains positions of overlapping peaks
  best_in_overlapping_chr=NULL
  
  chr = overlapping.peak.number$chr[1] # just extarct which chromosome we are working with
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



score1 = all_peaks$score[overlapping$queryHits]
score2 = all_peaks$score[overlapping$subjectHits]
which.to.keep <- score1>=score2
peaks.num.to.keeps <- unique(overlapping$queryHits[which.to.keep])
best_in_overlapping <- all_peaks[peaks.num.to.keeps,]
fwrite(best_in_overlapping,paste0(length(samples.id),'_recentered_Overlapping.mine.',add_filename,'.tsv'),sep='\t',
       row.names=FALSE)

pair_first <- all_peaks[overlapping$queryHits[which.to.keep],c('seqnames','start','end', 'score'), with = F]
pair_first[,row_1:= overlapping$queryHits[which.to.keep]]
dim(pair_first)
setnames(pair_first, old =c('seqnames','start','end', 'score'), 
         new = c('chr_1','st_1','en_1','score_1'))

pair_second <- all_peaks[overlapping$subjectHits[which.to.keep],c('seqnames','start','end', 'score'), with = F]
pair_second[,row_2:= overlapping$subjectHits[which.to.keep]]
dim(pair_second)
setnames(pair_second, old =c('seqnames','start','end', 'score'), 
         new = c('chr_2','st_2','en_2','score_2'))

pairs_all <- cbind(pair_first, pair_second)
rm(pair_first, pair_second, overlapping)
gc()

dim(pairs_all)
pairs_all <- pairs_all[score_1>=score_2,]
setorder(pairs_all, -score_1)
pairs_all <- pairs_all[,c('chr_1', 'score_1', 'row_1', 'score_2','row_2')]

#get non overlapping peaks
best_in_overlapping <- all_peaks[unique(pairs_all$row_1),]
best_in_overlapping <- best_in_overlapping[!duplicated(best_in_overlapping),]
fwrite(best_in_overlapping,paste0(length(samples.id),'_recentered_Overlapping.filtered.',add_filename,'.tsv'),sep='\t',
       row.names=FALSE)


#################
# first do iterative removal of overlapping peaks, keep the best score, discard ALL overlapping by cancer type
#
all_peaks.cancer <- split(all_peaks, by = 'Cancer', keep.by = T)
# filter out overlapping peaks via iterative removal 
recentered_final.cancer <- foreach(all_peaks=all_peaks.cancer)  %do% {
  iterative_removal(all_peaks)
}
# filter out peaks found in only 1 sample
recentered_final.cancer.reproducible <- foreach(recentered_final=recentered_final.cancer, all_peaks=all_peaks.cancer) %do% {
  filter_reproducible (recentered_final, all_peaks)
  }
#filter out peaks in chrY and by N in peak sequence
library(BSgenome.Hsapiens.UCSC.hg38)
recentered_final.cancer.reproducible.filtered <- lapply(recentered_final.cancer.reproducible, function(x) {
  x <- x[seqnames!='chrY',]
  filtered.x <- filter_N_peaks(x)
  return(filtered.x)
})
  



filter_N_peaks <- function(peak.dt) {
  gr <- StringToGRanges(peak.dt$new_peak, sep = c("-", "-")) #get GRanges object from peaks
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr) #extract fasta sequence
  names(seq) <- peak.dt$new_peak
  peaks.match.pattern <- vmatchPattern("N", seq) #match peak sequence with N in them
  peaks.withN <- names(peaks.match.pattern)[elementNROWS(peaks.match.pattern)>0] # these are peaks that contain N in their sequence
  toreturn <- peak.dt[! new_peak %in% peaks.withN,]
  fwrite(toreturn,paste0(length(samples.id),'_',cancer.type, '_recentered_final.reproducible.filtered.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  return()
}

filter_reproducible <- function(recentered_final, all_peaks) {
  cancer.type <- all_peaks$Cancer[1]
  print(cancer.type)
  current_peaks <- unique(recentered_final$new_peak)
  peaks.tokeep <- NULL
  all_peaks.current <- all_peaks[new_peak %in% current_peaks,]
  updated_peaks <- all_peaks.current[score.norm>=5, .N, by='new_peak'][N>=2,][['new_peak']] # this counts peak occurance across samples and selects peaks found in at least 2 samples
  toreturn <- all_peaks[new_peak %in% updated_peaks,] 
  fwrite(toreturn,paste0(length(samples.id),'_',cancer.type, '_recentered_final.reproducible.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  return(toreturn)
}

iterative_removal <- function(all_peaks) {
  cancer.type <- all_peaks$Cancer[1]
  print(cancer.type)
  recentered_p=StringToGRanges(all_peaks$new_peak, sep = c("-", "-"))
  cat(paste0('finding overlapping peaks in ',cancer.type,'\n'))
  overlapping=as.data.table(findOverlaps(recentered_p,recentered_p)) # find which peaks overlap
  overlapping=overlapping[queryHits!=subjectHits,]
  overlapping.peak.number <- unique(overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks table
  recentered_non_overlapping=all_peaks[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
  fwrite(recentered_non_overlapping,paste0(length(samples.id),'_',cancer.type, '_recentered_nonOverlapping.',add_filename,'.tsv'),
         sep='\t',row.names=FALSE)
  if (length(overlapping.peak.number)>0) {
    tmp <- data.table(chr = all_peaks$seqnames[overlapping.peak.number], 
                      num = overlapping.peak.number)
    overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) 
    registerDoParallel(cores=25)
    #this is where iterative removal of peaks is done
    best_in_overlapping_num <- foreach(overlapping.peak.number=overlapping.peak.number.split) %dopar% {
      cat('removing overlapping peaks\n')
      iterative_removal_core (overlapping.peak.number)
    }
    stopImplicitCluster()
    best_in_overlapping_num <- do.call('c', best_in_overlapping_num)
    best_in_overlapping_cancer <- all_peaks[best_in_overlapping_num,]
    fwrite(best_in_overlapping_cancer,paste0(length(samples.id),'_',cancer.type, '_recentered_Overlapping.',add_filename,'.tsv'),
           sep='\t',row.names=FALSE)
    recentered_final=rbindlist(list(recentered_non_overlapping,best_in_overlapping_cancer))
  } else {
    recentered_final=recentered_non_overlapping
  }
  fwrite(recentered_final,paste0(length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  return(recentered_final)
}

# this works like a charm
iterative_removal_core <- function(overlapping.peak.number) {
  chr = overlapping.peak.number$chr[1]
  
  running.vector <- overlapping.peak.number$num
  peaks.to.trash <- NULL
  peaks.to.keep <- NULL
  #i=1
  while (length(running.vector) != 0) {
    n <- running.vector[1] # this is the first and the best peak since peaks are sorted by scores
    # if (any(peaks.to.trash == n)) {
    #   i=i+1
    #   running.vector <- running.vector[running.vector != n]
    #   next()
    # } else {
    neighbor.peaks.num.discard <- overlapping[queryHits==n, subjectHits] #find positions of other peaks overlapping with the first one 
    running.vector <- setdiff(running.vector, neighbor.peaks.num.discard) # remove them from the list of peaks
    running.vector <- setdiff(running.vector, n)
    peaks.to.keep <- c(peaks.to.keep, n) # add this peak to the keeping list
    peaks.to.trash <- unique(c(peaks.to.trash, neighbor.peaks.num.discard)) # add neighbors to the list of peaks to discard
    #print(length(running.vector))
    #i=i+1
    #}
  }
  return(peaks.to.keep)
}




iterative_removal_test <- function(all_peaks) {
  cancer.type <- all_peaks$Cancer[1]
  print(cancer.type)
  all_peaks_unique <- distinct(all_peaks[,c('seqnames', 'new_peak', 'Cancer', 'score.norm')])
  recentered_p=StringToGRanges(all_peaks_unique, sep = c("-", "-"))
  cat(paste0('finding overlapping peaks in ',cancer.type,'\n'))
  overlapping=as.data.table(findOverlaps(recentered_p,recentered_p)) # find which peaks overlap
  overlapping=overlapping[queryHits!=subjectHits,]
  
  overlapping.peak.number <- unique(overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks table
  overlapping.peak <- all_peaks_unique[overlapping.peak.number]
  recentered_non_overlapping=all_peaks[! new_peak %in% overlapping.peak,] # select peaks that are not overlapping as non-overlapping peaks
  fwrite(recentered_non_overlapping,paste0(length(samples.id),'_',cancer.type, '_recentered_nonOverlapping.filtered.',add_filename,'.tsv'),
         sep='\t',row.names=FALSE)
  
  if (length(overlapping.peak.number)>0) {
    tmp <- data.table(chr = distinct(all_peaks[,c('seqnames', 'new_peak')])      all_peaks$seqnames[overlapping.peak.number], 
                      num = overlapping.peak.number)
    overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) 
    registerDoParallel(cores=25)
    best_in_overlapping_num <- foreach(overlapping.peak.number=overlapping.peak.number.split) %dopar% {
      cat('removing overlapping peaks\n')
      iterative_removal_core (overlapping.peak.number)
    }
    stopImplicitCluster()
    best_in_overlapping_num <- do.call('c', best_in_overlapping_num)
    best_in_overlapping_cancer <- all_peaks[best_in_overlapping_num,]
    fwrite(best_in_overlapping_cancer,paste0(length(samples.id),'_',cancer.type, '_recentered_Overlapping.filtered.',add_filename,'.tsv'),
           sep='\t',row.names=FALSE)
    recentered_final=rbindlist(list(recentered_non_overlapping,best_in_overlapping_cancer))
  } else {
    recentered_final=recentered_non_overlapping
  }
  fwrite(recentered_final,paste0(length(samples.id),'_',cancer.type, '_recentered_final.filtered.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  return(recentered_final)
}

