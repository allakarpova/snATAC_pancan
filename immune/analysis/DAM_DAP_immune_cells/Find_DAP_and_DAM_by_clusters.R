# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################

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
suppressMessages(library(org.Hs.eg.db))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

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
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
panc <- readRDS(input.path)

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

if(!file.exists(paste0('DAM_findAllMarkers_by_clusters_',add_filename,'.txt'))) {
  Idents(panc) <- 'seurat_clusters'
  DefaultAssay(panc) <- 'chromvar'
  unique(Idents(panc))
  
  da_motif <- FindAllMarkers(
    object = panc,
    only.pos = F,
    test.use = 'LR',
    latent.vars = 'nCount_peaks'
  )
  
  fwrite(da_motif, paste0('DAM_findAllMarkers_by_clusters_',add_filename,'.txt'), sep = '\t', row.names = T)
}


# annotate peaks
all_m <- data.frame (panc@assays$ATAC_immune@ranges)
all_m$m_coords=paste0(all_m$seqnames,"-",all_m$start,"-",all_m$end)
m_all_m <- StringToGRanges(all_m$m_coords)
###Annotate the Motifs with the closest gene, to get the info about promoter regions:
peakAnno <- annotatePeak(m_all_m, tssRegion=c(-1500, 500),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
peak.annotation <- data.frame (peakAnno)
peak.annotation$coord <- paste0(peak.annotation$seqnames,"-",peak.annotation$start,"-",peak.annotation$end)

Idents(panc) <- 'seurat_clusters'
DefaultAssay(panc) <- 'ATAC_immune'
unique(Idents(panc))

da_peaks <- FindAllMarkers(
  object = panc,
  only.pos = F,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
da_peaks <- merge(da_peaks, peak.annotation[c('SYMBOL', 'annotation' ,'GENENAME', 'coord')], by.x = 0, by.y = 'coord', all.x = T)
fwrite(da_peaks, paste0('DAP_findAllMarkers_by_clusters_',add_filename,'.txt'), sep = '\t', row.names = F)










