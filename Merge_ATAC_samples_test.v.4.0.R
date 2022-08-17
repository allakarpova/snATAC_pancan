## merge ATAC samples which can be either regular ATAC seq or combo ATAC seq
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)
library(optparse)
library(data.table)

library(googlesheets4)
library(stringr)

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
              default="peaks", 
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
plan("multiprocess", workers = 16)
options(future.globals.maxSize = 1000000 * 1024^2) # for 50 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE') %>% dplyr::filter(`Disease Type` == 'BRCA')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data folder`)
#samples <- fread(sample.path, data.table = F, header = F)
samples.id <- samples$Sample %>% as.character()
samples.type <- samples$`Data Type` %>% as.character()
samples$`Disease Type`
#########
cat (paste("Samples found:" ,length(samples.id), '\n'))


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
cat ('Reading in objects\n')
atac=vector(mode = "list", length = length(samples.id))
for (i in 1:length(samples.id)){
  print(samples.id[i])
  atac[[i]]=readRDS(list.files(path = input.path, full.names = T, pattern = paste0(str_split_fixed(samples.id[i], '_',2)[2],'.*rds'), all.files = T, recursive = T))
  DefaultAssay(atac[[i]]) <- assay.towork
  atac[[i]] <- DietSeurat(atac[[i]], assay = assay.towork)
}

#####To obtain the best results - use ALL peaks!
cat ('Combining peaks\n')

combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
peaks.use=combined.peaks

#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.
cat ('creating matrix counts\n')
matrix.counts=vector(mode = "list", length = length(samples.id))

for (i in 1:length(samples.id)){
  matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
  ) 
}

cat ('creating peaksinters and removing useless assays\n')
for (i in 1:length(samples.id)){
  atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],fragments=Fragments(atac[[i]]@assays$X500peaksMACS2))
  atac[[i]]$dataset=samples.id[i]
  DefaultAssay(atac[[i]])<-'peaksinters'
  ###remove original assay
  atac[[i]][['X500peaksMACS2']]<-NULL
}


####Merging and reduction on old peaks
cat ('Merging\n')
combined <- merge(x = atac[[1]], y = atac[2:length(samples.id)], add.cell.ids = samples.id)
DefaultAssay(combined) <- "peaksinters"
#remove individual objects
rm(atac)
gc()

cat ('Normalizing and reduction\n')
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
  resolution = 1,
  verbose = FALSE
)

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_",add_filename,".rds"))


####2021-03-20: Change Peaks to MACS2
###########################################################
############MACS2 peak calling#############################
###########################################################
DefaultAssay(combined) <- 'peaksinters'
peaks <- CallPeaks(
  object = combined,
  macs2.path=opt$m
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

p=as.data.frame(peaks)
write.table(p,paste0('MACS2_peaks.',add_filename,'.tsv'),sep='\t',quote=FALSE,
            row.names=FALSE)

# some peaks happen to be in a wierd chromosome pieces
#p <- p %>% filter(seqnames %in% standardChromosomes(grange.counts))
p$peak_center=p$start+p$relative_summit_position
p$recentered_start=p$peak_center-250
p$recentered_end=p$peak_center+250

####Now check that new start and end don't go beyond the chromosome boundaries
chr_size=read.table(path.to.chrom.size,sep='\t',header=FALSE)
colnames(chr_size)=c('seqnames','chr_length')
p1=merge(p,chr_size,all.x=TRUE)

p1=p1[p1$recentered_end<=p1$chr_length && p1$recentered_start>=0,]
p1$length=p1$recentered_end-p1$recentered_start+1
p1$new_peak=paste0(p1$seqnames,":", p1$recentered_start, '-',p1$recentered_end)

recentered_p=StringToGRanges(p1$new_peak, sep = c(":", "-"))

olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
olap1=olap[olap$queryHits!=olap$subjectHits,]

recentered_non_olap=p1[-olap1$queryHits,]
recentered_olap=p1[olap1$queryHits,]

pairs=cbind(p1[olap1$queryHits,c(1:3,7)],olap1$queryHits,p1[olap1$subjectHits,c(1:3,7)],olap1$subjectHits)
colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

pairs=pairs[pairs$score_1>=pairs$score_2,]
all_st=NULL
for (i in 1:nrow(pairs)){
  if (nrow(pairs)>0){
    pairs=pairs[order(-pairs$score_1),]
    if (pairs[1,4]>=pairs[1,9]){
      all_st=rbind(all_st,p1[rownames(p1)==pairs[1,5],])
      pairs=pairs[-1,]
      pairs=pairs[pairs$row_1!=pairs[1,10],]
    }else{
      all_st=rbind(all_st,p1[rownames(p1)==pairs[1,10],])
      pairs=pairs[-1,]
      pairs=pairs[pairs$row_1!=pairs[1,5],]
    }
  }
}
all_st=as.data.frame(all_st)
all_st=all_st[!duplicated(all_st),]

recentered_final=rbind(recentered_non_olap,all_st)
write.table(recentered_final,paste0('recentered_final.filtered',add_filename,'.tsv'),sep='\t',
            quote=FALSE,row.names=FALSE)

recentered_p=StringToGRanges(recentered_final$new_peak, sep = c(":", "-"))
matrix.counts <- FeatureMatrix(
  fragments = Fragments(combined@assays$peaksinters),
  features = recentered_p,
  sep = c("-","-"),
  cells = colnames(combined)
)

combined[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                     annotation = annotations,
                                                     genome = 'hg38',
                                                     fragments = Fragments(combined@assays$peaksinters))

DefaultAssay(combined)<-'X500peaksMACS2'
# remove previous assay
combined[['peaksinters']] <- NULL

#add some more QC stuff
peak.data <- GetAssayData(object = combined, assay = 'X500peaksMACS2', slot = "counts")
total_fragments_cell <- combined$atac_fragments
peak.counts <- colSums(x = peak.data)
frip <- peak.counts *100 / total_fragments_cell
combined <- AddMetaData(object = combined, metadata = frip, col.name = 'pct_read_in_peaks_500MACS2')
combined <- AddMetaData(object = combined, metadata = peak.counts, col.name = 'peak_region_fragments_500MACS2')

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
  resolution = 1,
  verbose = FALSE
)

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples.id),"_snATAC_Merged_recalled_peaks_",add_filename,".rds"))

cat('plotting...\n')
p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1) + 
  ggplot2::ggtitle("Combined snATAC samples.id")+
  scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),"grey"))

p2 <- DimPlot(combined, pt.size = 0.1,label=T)
pdf(paste0(length(samples.id),"_snATAC_Merged_", add_filename, ".pdf"),height=6,width=16, useDingbats = F)
p1+p2
dev.off()

write.table(combined@meta.data, paste0(length(samples.id),"_snATAC_Merged_recalled_peaks_",add_filename,"_metaData.txt"),sep="\t",quote=FALSE, row.names = T)
write.table(samples.id,paste0("Samples_snATAC_Merged_",add_filename,".txt"),sep="\t",quote=FALSE)
