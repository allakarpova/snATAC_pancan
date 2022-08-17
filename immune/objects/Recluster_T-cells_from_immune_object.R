# Subset T-cells and NK cells and recluster from immune object
library(future)
###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

library(Signac)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(tibble)
library(reshape)
library(data.table)

library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicRanges)
library(optparse)
library(patchwork)

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
              default='X500peaksMACS2',
              help = "X500peaksMACS2 or peaksinters",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
cell_column <- opt$cell_type_column
macs2_path <- opt$macs2_path
chrom.size <- opt$chrom_size
meta.path <- opt$metadata.file
assay.towork <- opt$assay

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
  panc <- AddMetaData(panc, metadata = my.metadata)
}

panc$is_immune <- grepl('T-cell|CD8|CD4|NK|Treg|Lympho', as.character(unlist(panc[[cell_column]])))

dim(panc)
#subsetting
cat('subsetting\n')
panc.my <- subset(x = panc, subset = is_immune)
dim(panc.my)

rm(panc)
gc()

if (file.exists(paste0('recentered_final.filtered_T-cells_',add_filename,'.tsv'))) {
  recentered_final <- read.table(paste0('recentered_final.filtered_T-cells_',add_filename,'.tsv'),sep='\t',
                                 header = T)
} else {
  
  ###########################################################
  ############MACS2 peak calling#############################
  ###########################################################
  peaks <- CallPeaks(
    object = panc.my,
    macs2.path=macs2_path
  )
  
  p=as.data.frame(peaks)
  write.table(p,paste0('MACS2_peaks.T-cells_',add_filename,'.tsv'),sep='\t',quote=FALSE,
              row.names=FALSE)
  
  p$peak_center=p$start+p$relative_summit_position
  p$recentered_start=p$peak_center-250
  p$recentered_end=p$peak_center+250
  
  ####Now check that new start and end don't go beyond the chromosome boundaries
  chr_size=read.table(chrom.size,sep='\t',header=FALSE)
  colnames(chr_size)=c('seqnames','chr_length')
  p1=merge(p,chr_size,all.x=TRUE)
  
  p1=p1[p1$recentered_end <= p1$chr_length && p1$recentered_start>=0,]
  p1$length=p1$recentered_end - p1$recentered_start+1
  p1$new_peak=paste(p1$seqnames,p1$recentered_start,p1$recentered_end,sep='-')
  
  recentered_p=StringToGRanges(p1$new_peak, sep = c("-", "-"))
  
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
  write.table(recentered_final,paste0('recentered_final.filtered_T-cells_',add_filename,'.tsv'),sep='\t',
              quote=FALSE,row.names=FALSE)
}


recentered_p=StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
frag <- Fragments(panc.my@assays[[assay.towork]])
#this will remove fragment objects with just 1 or 0 cells because they fail FeatureMatrix function
frag.filtered <- frag[do.call( 'c', lapply(frag, function(x) length(x@cells) > 1))]

cell.count <- panc.my$dataset %>% count()
exlude.samples <- cell.count$x[cell.count$freq <= 1]
samples <- unique(panc.my$dataset)
samples.to.keep <- samples[!(samples %in% exlude.samples)]
#now exclude these samples from the object
panc.my <- subset(panc.my, subset = dataset %in% samples.to.keep)


matrix.counts <- FeatureMatrix(
  fragments = frag.filtered,
  features = recentered_p,
  sep = c("-","-"),
  cells = colnames(panc.my),
)

panc.my[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                    fragments=frag.filtered)
DefaultAssay(panc.my)<-'X500peaksMACS2'

peak.data <- GetAssayData(object = panc.my, assay = 'X500peaksMACS2', slot = "counts")
total_fragments_cell <- panc.my$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
panc.my <- AddMetaData(object = panc.my, metadata = frip, col.name = 'frip_500MACS2')
panc.my <- AddMetaData(object = panc.my, metadata = peak.counts, col.name = 'peak_RF_500MACS2')

##########################
##### RECLUSTERING #######
##########################

panc.my <- RunTFIDF(panc.my)
panc.my <- FindTopFeatures(panc.my, min.cutoff = 20)
panc.my <- RunSVD(
  panc.my,
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work = 400
)

panc.my <- RunUMAP(panc.my, dims = 2:30, reduction = 'lsi')
panc.my <- FindNeighbors(
  object = panc.my,
  reduction = 'lsi',
  dims = 2:30
)
panc.my <- FindClusters(
  object = panc.my,
  algorithm = 3,
  resolution = 1,
  verbose = FALSE
)

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
saveRDS(panc.my, paste0("Reclustered_T-cells_snATAC_Merged_",add_filename,".rds"))


cat('plotting...\n')
p1 <- DimPlot(panc.my, group.by = 'dataset', pt.size = 0.1) + 
  ggplot2::ggtitle("Combined snATAC samples")
p2 <- DimPlot(panc.my, pt.size = 0.1,label=T)
p3 <- DimPlot(panc.my,group.by = cell_column, pt.size = 0.1,label=T)
pdf(paste0("Reclustered_T-cells_snATAC_Merged_", add_filename, ".pdf"),height=30,width=17, useDingbats = F)
print(p1/p2/p3)
dev.off()

p <- DimPlot(panc.my, group.by = 'Disease', pt.size = 0.1)
pdf(paste0("Dimplot_disease_reclustered_T-cells_snATAC_Merged_", add_filename, ".pdf"),height=10,width=12, useDingbats = F)
print(p)
dev.off()

