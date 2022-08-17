#####author: Nadezhda V. Terekhanova
#####2020-09-26: v.3.0 - this version of the script was adapted for the latest Signac v.1.0.0
###do prior to running the script `export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp` - otherwise it crushes

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

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2) # for 50 Gb RAM

samples=c('TWCE-HT-027-B-Slice1fresh_ATAC-lib1_V2', 'TWCE-HT029B1-XBa1', 'TWCE-HT035B-XBa1', 'TWCE-HT088B1-S1H1A2K2Y2N1-ATAC', 'TWCE-HT088B1-S1H2A2Y1N1-ATAC', 'TWCE-HT137B1-XBa1', 'TWCE-HT141B1-XBa1', 'TWCE-HT163B1-S1H6A3-XBa1', 'TWCE-HTAN_1408-06-ATAC-lib1_V2','TWCE-HT206B1-XBa1_1-HT206B1-XBa1','TWCE-HT128B1-XBa2_1-HT128B1-XBa2','TWCE-HT214B1-XBa1', 'TWCE-HT217B1-XBa1')

atac=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    atac[[i]]=readRDS(paste("inputs/",samples[i],"_cellTyped.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'peaks'

}

#####To obtain the best results - use ALL peaks!

combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
peaks.use=combined.peaks


#For testing purposes only:
#peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.

matrix.counts=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$peaks),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
    ) 
}


for (i in 1:length(samples)){
atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],fragments=Fragments(atac[[i]]@assays$peaks))
atac[[i]]$dataset=samples[i]
DefaultAssay(atac[[i]])<-'peaksinters'
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)

DefaultAssay(combined) <- "peaksinters"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(
  combined,
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work = 400
)
	
#combined <- RunUMAP(combined, dims = 1:30, reduction = 'lsi')
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1) + 
ggplot2::ggtitle("Combined snATAC samples")+
scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),"grey"))

#combined$predicted.gr=combined$predicted.id
#combined$predicted.gr=ifelse(combined$predicted.id %in% c("B","DC","NK","Macrophage","Mast","CD4_T","CD8_T","Plasma","Treg"),"Immune",combined$predicted.id)
#combined$predicted.gr=ifelse(combined$predicted.id %in% c("Endothelial","Fibroblast"),"Stroma",combined$predicted.id)
#combined$predicted.gr=ifelse(combined$predicted.id %in% c("Endothelial","Fibroblast"),"Immune",combined$predicted.id)


p3 <- DimPlot(combined, group.by = 'predicted.id', pt.size = 0.1,label=FALSE) + 
ggplot2::ggtitle("Cell Types snATAC")+
scale_color_manual(values=c(brewer.pal(n = 12, name = "Set3"),"grey"))

pdf(paste(length(samples),"_snATAC_Merged_Breast_HTAN.pdf",sep=""),height=6,width=16)
p1+p3
dev.off()

write.table(samples,"Samples_snATAC_Merged_Breast_HTANrds.txt",sep="\t",quote=FALSE)

###Run ChromVar using RunChromVar.R

saveRDS(combined, paste(length(samples),"_snATAC_Merged_Breast_HTAN.rds",sep=""))
