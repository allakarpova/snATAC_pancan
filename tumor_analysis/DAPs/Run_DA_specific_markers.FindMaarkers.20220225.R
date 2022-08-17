#Run for BRCA and BRCA_Basal separate tests:
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(presto)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(future)
plan("multiprocess", workers =4)
options(future.globals.maxSize = 50 * 1024^3)

atac=readRDS(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/',
'Merge.vers.20220207/PanCan_object/Tumor_Normal.v.20220220/159_snATAC_129K_peaks_TumorNormal.motifsAdded.chromvar.20220221.rds.gz',
sep=''))

atac$Disease=ifelse(atac$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3"), "BRCA_Basal",
atac$Disease)


ATAC=subset(atac, cell_type.harmonized.cancer=='Tumor')


###Try with FindMarkers:
peak.data <- GetAssayData(object = ATAC, assay = 'pancan_s', slot = "counts")
total_fragments_cell <- ATAC$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
ATAC <- AddMetaData(object = ATAC, metadata = frip, col.name = 'frip_pancan_s')
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_pancan_s')
DefaultAssay(ATAC)='pancan_s'


###calculate in parallel:
cell_types=unique(ATAC$Disease)
Idents(ATAC)=ATAC$Disease

cell_types=c("CRC","BRCA","MM","GBM")
cell_types=c("HNSCC", "UCEC","PDAC","CESC")
cell_types=c("ccRCC","OV","BRCA_Basal")


all_da_peaks=NULL
for (cell_type in cell_types){
da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_type,
#  ident.2='', 
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'peak_RF_pancan_s'
)

da_peaks$Disease=cell_type
da_peaks$peak=rownames(da_peaks)
all_da_peaks=rbind(all_da_peaks,da_peaks)
print(cell_type)
}

write.table(all_da_peaks, paste("out/da_peaks_oneCances_vs_Others.minPct0.1.part3.20220226.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
write.table(all_da_peaks, paste("out/da_peaks_oneCances_vs_Others.minPct0.1.part2.20220226.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
write.table(all_da_peaks, paste("out/da_peaks_oneCances_vs_Others.minPct0.1.part1.20220226.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)

###ENDS HERE FOR NOW####

#################################################################
###Also try 1 sample vs pooled all other cancer types' samples###
#################################################################

###Now try running presto! using normalized counts:
###Maybe complicated...
###try filtering instead.







annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atac) <- annotations

Idents(atac)=factor(Idents(atac),levels=c('CESC','HNSCC','BRCA','PDAC','CRC','ccRCC','GBM','MM','UCEC','OV'))
new_peak="chr8-127650000-127750000"
p=CoveragePlot(
  object = atac,
  region = new_peak,
  annotation = TRUE,
   peaks = FALSE,
  links=FALSE
)

pdf(paste("plots/CoverageByDisease_",new_peak,".v4.pdf",sep=""),width=8,height=5.5,useDingbats=FALSE)
print(p)
dev.off()