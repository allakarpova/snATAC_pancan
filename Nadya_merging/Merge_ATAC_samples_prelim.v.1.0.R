#####author: Nadezhda V. Terekhanova
#####2020-09-26: v.3.0 - this version of the script was adaapted for the latest Signac v.1.0.0
#####2020-07-01##################
#####v.2.0 - "auto", now just need to provide the sample-ids. Script will create the list of Seurat objects, and will store the provided sample ids in the $dataset (meta-data).
#####Output - Seurat object of the merged datasets and list of samples' ids used.

###############
#IMPORTANT: do prior to running the script `export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp` - otherwise it crushes (probably because not enough space to write everything)
###############

library(Signac)
library(Seurat)
library(GenomeInfoDb)
###library(EnsDb.Hsapiens.v75) ###Now using newer version of the annotation:
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

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM
date='20210706'

n_samples=19
samples_atac=c('1408-06', 'TWCE-HT029B1-XBa1', 'TWCE-HT035B-XBa1', 'TWCE-HT088B1-S1H1A2K2Y2N1-ATAC',
 'TWCE-HT088B1-S1H2A2Y1N1-ATAC', 'TWCE-HT128B1-XBa2_1-HT128B1-XBa2', 'TWCE-HT137B1-XBa1', 'TWCE-HT141B1-XBa1',
 'TWCE-HT163B1-S1H6A3-XBa1', 'TWCE-HT206B1-XBa1_1-HT206B1-XBa1', 'TWCE-HT214B1-XBa1', 'TWCE-HT217B1-XBa1')
samples_multiome=c('BRCA_HT235B1-S1H1Fc2A2N1Bma1_1', 'BRCA_HT243B1-S1H4Fc2A2N1Bma1_1',
'BRCA_HT263B1-S1H1A3N1Ba1_1', 'BRCA_HT271B1-S1H3Fc2A5N1Bma1_1', 'BRCA_HT305B1-S1H1Fc2A2_1N1Bma1_1',
 'TWCE-HT268B1-Th1H3XBmn1_1', 'TWCE-HT297B1-S1H1Fc2A2XBmn1_1')


tab=read.table('Case_id_table.20210707.txt',sep='\t',header=T)
rownames(tab)=tab$Sample
tab1=tab[samples_atac,]
piece_ids_atac=tab1$Case

tab1=tab[samples_multiome,]
piece_ids_multiome=tab1$Case

samples=c(samples_atac,samples_multiome)
piece_ids=c(piece_ids_atac,piece_ids_multiome)

atac=vector(mode = "list", length = length(samples))

for (i in 1:length(samples_atac)){
    atac[[i]]=readRDS(paste("../1.Create_rds/out/",samples[i],"/",samples[i],"_processed_atac.rds",
sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['ATACGeneActivity']]<-NULL
    atac[[i]][['peaks']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}
for (i in (length(samples_atac)+1):(length(samples_atac)+length(samples_multiome))){
    atac[[i]]=readRDS(paste("../../../../Multiome/BR_HTAN/1.Create_rds/out/",samples[i],"/",
samples[i],"_processed_multiome.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['ATACGeneActivity']]<-NULL
    atac[[i]][['ATAC']]<-NULL
    atac[[i]][['RNA']]<-NULL
    atac[[i]][['SCT']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}

#####To obtain the best results - use ALL peaks!

combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#peaks.use=combined.peaks
####Now using MACS2-peak calling:
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#For testing purposes only:
#peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.

matrix.counts=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
    ) 
}


for (i in 1:length(samples)){
atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],
fragments=Fragments(atac[[i]]@assays$X500peaksMACS2))
atac[[i]]$dataset=samples[i]
DefaultAssay(atac[[i]])<-'peaksinters'
###remove other assay
#atac[[i]][['X500peaksMACS2']]<-NULL
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
saveRDS(combined, paste('19_BR_snATAC.',date,'.rds',sep=''))





