####Important!!! --to limit number of cores:
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(30)) #number of cores to use - otherwise it crushes

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

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
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

#atac=readRDS('26_ccRCC_snATAC.selectedPeaks.chromvar.CICERo.v6.20210512.rds')
DefaultAssay(atac)='peaksMACS2'
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(atac,"26_ccRCC_snATAC.selectedPeaks.chromvar.CICERo.MotifsAdded.v7.20210520.rds")

atac=readRDS('26_ccRCC_snATAC.selectedPeaks.chromvar.v3.20210602.rds')

x=Motifs(atac[['peaksMACS2']])
mot=as.data.frame(x@positions)
peaks=as.data.frame(rownames(atac@assays$peaksMACS2))
peaks$Peaks_souce="Merged_snATAC"


colnames(peaks)[1]='Peak'
mot$motif_cooord=paste(mot$seqnames,mot$start,mot$end,sep="-")
in_peaks=StringToGRanges(peaks$Peak, sep = c("-", "-"))
out_peaks=StringToGRanges(mot$motif_cooord, sep = c("-", "-"))
olap=as.data.frame(findOverlaps(in_peaks,out_peaks))
#pairs=cbind(da_p[olap$queryHits,],olap$queryHits,peaks[olap$subjectHits,],olap$subjectHits)
pairs=cbind(peaks[olap$queryHits,],mot[olap$subjectHits,])
pairs_1=pairs %>% dplyr::select ('Peak','group_name','strand','score','motif_cooord')
colnames(pairs_1)[5]= 'motif_coord'
jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',header=TRUE)
colnames(jaspar)[1]='group_name'
pairs_2=merge(pairs_1,jaspar,all.x=TRUE)
write.table(pairs_2,paste('Motifs_matched.26_snATAC_merged.object.20210531.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(mot,paste('All_motifs_before_matching.26_snATAC_merged.object.20210531.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
peaks_sel=read_delim('../../6.DA_motifs/Match_motifs/20210520/DEG_associated_Peaks.20210517.v1.tsv',delim='\t')
peaks_sel=as.data.frame(peaks_sel)
peaks_sel_1=merge(peaks_sel,pairs_2,all.x=TRUE)
write.table(peaks_sel_1, 'Motifs_matched.DEG_associated_Peaks.20210517.v1.tsv',sep="\t",row.names=FALSE,quote=FALSE)
