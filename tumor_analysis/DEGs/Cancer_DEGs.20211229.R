#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
RhpcBLASctl::blas_set_num_threads(50)
library(future)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50 * 1024^3)

date='20211228'

rna=readRDS(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/10.snRNA/1.Merge_snRNA/',
'Merged_132_snRNA.SCTrandformed.3KVarFeatures.200PCs.v2.20211228.rds',sep=''))
rna_all=rna
rna=subset(rna_all, cell_type.harmonized.cancer=='Tumor')

rna$Piece_ID=paste(rna$Disease,rna$Piece_ID, sep='_')
rna$Disease=ifelse(rna$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3"), "BRCA_Basal",
rna$Disease)
dis_cols_ed=c("BRCA"= "#E9967A","BRCA_Basal"="#C70039", "CESC"="#FFFF00", "CRC"="#FF8C00", "GBM"="#6A3D9A",
"HNSCC"="#FF69B4", "MM"="#A65628", "OV"="#57C785", "PDAC"="#80B1D3","UCEC"="#1F78B4", "ccRCC"="#C0C0C0")

cancers=unique(rna$Disease)
Idents(rna)=rna$Disease
DefaultAssay(rna)='SCT'
all_degs=NULL

for (disease in cancers){
degs <- FindMarkers(
  object = rna,
  ident.1 = disease,
#  ident.2='',
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  assay='SCT',
  logfc.threshold=0
)

degs$Disease=disease
degs$Gene=rownames(degs)
all_degs=rbind(all_degs,degs)
print(disease)
}
write.table(all_degs, paste("out/degs_oneCances_vs_Others.minPct0.1.all.20211229.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)

write.table(all_degs, paste("out/degs_oneCances_vs_Others.minPct0.1.part1.20211229.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
