###We want to re-size peaaks relative to summits +/-250bp (the same as in the ChromVar and CGA ATAC paper)

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 100 * 1024 ^ 2)

#######################################################
#############Now try the merged object, all samples:###
#######################################################

date='20210706'
n_samples=19
samples_atac=c('1408-06', 'TWCE-HT029B1-XBa1', 'TWCE-HT035B-XBa1', 'TWCE-HT088B1-S1H1A2K2Y2N1-ATAC',
 'TWCE-HT088B1-S1H2A2Y1N1-ATAC', 'TWCE-HT128B1-XBa2_1-HT128B1-XBa2', 'TWCE-HT137B1-XBa1', 'TWCE-HT141B1-XBa1',
 'TWCE-HT163B1-S1H6A3-XBa1', 'TWCE-HT206B1-XBa1_1-HT206B1-XBa1', 'TWCE-HT214B1-XBa1', 'TWCE-HT217B1-XBa1')
samples_multiome=c('BRCA_HT235B1-S1H1Fc2A2N1Bma1_1', 'BRCA_HT243B1-S1H4Fc2A2N1Bma1_1', 
'BRCA_HT263B1-S1H1A3N1Ba1_1', 'BRCA_HT271B1-S1H3Fc2A5N1Bma1_1', 'BRCA_HT305B1-S1H1Fc2A2_1N1Bma1_1',
 'TWCE-HT268B1-Th1H3XBmn1_1', 'TWCE-HT297B1-S1H1Fc2A2XBmn1_1')

tab=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/',
'snRNA_processed_Yige/meta_data.20200505.v1.tsv',sep=''),sep='\t',header=TRUE)
tab=tab[tab$Aliquot.snRNA %in% samples,]
samples=tab$Aliquot.snRNA
piece_ids=tab$Aliquot.snRNA.WU

all_peaks=NULL
for (sample in samples_atac){
    peaks=read.table(paste("../1.Create_rds/out/",sample,"/recentered_final.filtered",sample,".tsv",
sep=""),sep='\t',header=TRUE)
    peaks$Sample=sample
    all_peaks=rbind(all_peaks,peaks)
}
for (sample in samples_multiome){
    peaks=read.table(paste("../../../../Multiome/BR_HTAN/1.Create_rds/out/",sample,
"/recentered_final.filtered",sample,".tsv",
sep=""),sep='\t',header=TRUE)
    peaks$Sample=sample
    all_peaks=rbind(all_peaks,peaks)
}
            
p=all_peaks
write.table(p,paste('peaks/MACS2_peaks.19BR_samples.BySample.',date,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)
p=read.table(paste('peaks/MACS2_peaks.19BR_samples.BySample.',date,'.tsv',sep=''),sep='\t',header=TRUE)


p1=p
#p1$length=p1$recentered_end-p1$recentered_start+1
#p1$new_peak=paste(p1$seqnames,p1$recentered_start,p1$recentered_end,sep='-')

recentered_p=StringToGRanges(p1$new_peak, sep = c("-", "-"))

olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
olap1=olap[olap$queryHits!=olap$subjectHits,]
#write.table(recentered_non_olap,'peaks/recentered_non_overlapping.tsv',sep='\t',quote=FALSE,row.names=FALSE)
#write.table(recentered_olap,'peaks/recentered_overlapping.tsv',sep='\t',quote=FALSE,row.names=FALSE)

recentered_non_olap=p1[-olap1$queryHits,]
#recentered_olap=p1[olap1$queryHits,]

pairs=cbind(p1[olap1$queryHits,c(1:3,7)],olap1$queryHits,p1[olap1$subjectHits,c(1:3,7)],olap1$subjectHits)
colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs=pairs[order(-pairs$score_1),]
pairs_all=pairs

library(doParallel)
registerDoParallel(cores=30)

all_st=NULL
all_st<-foreach(chr_n=c(1:22,"X","Y")) %dopar% {
chr=paste("chr",chr_n,sep='')
pairs=pairs_all[pairs_all$chr_1==chr,]
pairs=pairs[,c(4,5,9,10)]
all_st_chr=NULL
for (i in 1:nrow(pairs)){
    if (nrow(pairs)>0){
    p_del=pairs[pairs$row_1==pairs[1,2],]
    all_st_chr=rbind(all_st_chr,p1[rownames(p1)==pairs[1,2],])
    pairs=pairs[!(pairs$row_1 %in% c(p_del$row_1[1],p_del$row_2)),]
    }
#    print(paste(chr,nrow(pairs),sep=' '))
}
return(all_st_chr)
}

all_st_f=NULL
for (i in 1:24){
    all_st_1=as.data.frame(all_st[[i]])
    all_st_1=all_st_1[!duplicated(all_st_1),]
    all_st_f=rbind(all_st_f,all_st_1)
}


recentered_final=rbind(recentered_non_olap,all_st_f)
write.table(recentered_final,paste('peaks/recentered_final.filtered.',date,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(recentered_non_olap,paste('peaks/recentered_nonOverlapping.filtered.',date,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)
write.table(all_st_f,paste('peaks/recentered_Overlapping.filtered.',date,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)


##########################################################
####Now create FeatureMatrix with the new set of peaks:###
##########################################################
recentered_final=read.table(paste('peaks/recentered_final.filtered.',date,'.tsv',sep=''),
sep='\t',header=T)
atac=readRDS(paste('19_BR_snATAC.',date,'.rds',sep=''))

recentered_p=StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
matrix.counts <- FeatureMatrix(
    fragments = Fragments(atac@assays$peaksinters),
    features = recentered_p,
    sep = c("-","-"),
    cells = colnames(atac)
)

###next is here:
atac[['peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
fragments=Fragments(atac@assays$peaksinters))
DefaultAssay(atac)<-'peaksMACS2'

#atac[['X500peaksMACS2']]<-NULL
atac[['peaksinters']]<-NULL

###Overlapping ranges supplied -- check this later:

###Remove some assays
#atac[['RNA']]<-NULL
#atac[['chromvar']]<-NULL
saveRDS(atac, '26_ccRCC_snATAC.selectedPeaks.20210531.rds')

###################
#ends here for now#
###################


atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:50)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:50)


###this helps:
options(future.globals.maxSize= 891289600)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)

#warning can be ignored, see here: https://github.com/timoast/signac/issues/266

###Run ChromVar and then save:
saveRDS(atac,'19_BR_snATAC.selectedPeaks.chromvar.20210707.rds.gz',compress=TRUE)


p1 <- DimPlot(atac, group.by = 'Piece_ID', pt.size = 0.1) +
ggplot2::ggtitle("Combined snATAC samples")+
scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 7, name = "Dark2")))


pdf(paste("19_snATAC_Merged_BR_CPTAC.MACS2peaks_Sample.20210503.pdf",sep=""),height=6,width=8)
p1
dev.off()


#####ENDS HERE###


annot=read.table('Update_annotation/26_snATAC_ccRCC_samples.20210505.tsv',sep='\t',header=TRUE)
orig_1=as.data.frame(atac$Piece_ID)
annot=annot[rownames(orig_1),]
atac$cell_type_manual_5=annot$cell_type_manual_5

s_info=read.table('Sample_categories.20210503.txt',sep='\t',header=TRUE)
orig_1=as.data.frame(atac$Piece_ID)
colnames(orig_1)='Aliquot.WU'
orig_1$Barcode=rownames(orig_1)

s_inf=merge(orig_1,s_info,all.x=TRUE)
rownames(s_inf)=s_inf$Barcode
s_inf=s_inf[rownames(orig_1),]

atac$Category=s_inf$Category
meta=atac@meta.data
write.table(meta,'26_ccRCC_snATAC.selectedPeaks.20210531.tsv',sep='\t',quote=FALSE)
saveRDS(atac,'26_ccRCC_snATAC.selectedPeaks.chromvar.v2.20210531.rds'

########################
###Saving the object:###
########################

combined=atac

p1 <- DimPlot(combined, group.by = 'Piece_ID', pt.size = 0.1) +
ggplot2::ggtitle("Combined snATAC samples")+
scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 8, name = "Dark2"),
brewer.pal(n = 6, name = "Set1")))

#n_colors=length(unique(combined$predicted.id))
p3 <- DimPlot(combined, group.by = 'cell_type_manual_5', pt.size = 0.1,label=FALSE) +
ggplot2::ggtitle("Cell Types snATAC")+
scale_color_manual(values=c(brewer.pal(n = 12, name = "Set3"),brewer.pal(n = 8, name = "Dark2"),"grey"))

pdf(paste("26_snATAC_Merged_ccRCC_CPTAC.MACS2peaks_Sample.20210503.pdf",sep=""),height=6,width=16)
p1+p3
dev.off()


p3 <- DimPlot(combined, group.by = 'Category', pt.size = 0.1,label=FALSE) +
ggplot2::ggtitle("Cell Types snATAC")+scale_color_manual(values=c(brewer.pal(n = 5, name = "Dark2")))

pdf(paste("26_snATAC_Merged_ccRCC_CPTAC.Category.MACS2peaks_Sample.20210503.pdf",sep=""),height=6,width=16)
p1+p3
dev.off()

p3 <- DimPlot(combined, group.by = 'seurat_clusters', pt.size = 0.1,label=TRUE) +
ggplot2::ggtitle("Cell Types snATAC")
#+scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 4, name = "Dark2")))

pdf(paste("26_snATAC_Merged_ccRCC_CPTAC.Clusters.MACS2peaks_Sample.20210503.pdf",sep=""),height=6,width=16)
p1+p3
dev.off()




###############
#ENDS here#####
###############



pdf('Pct_in_peaks_11_GBM_snATAC_samples.pdf',width=8, height=6,useDingbats=FALSE)
FeaturePlot(atac,features='pct_reads_in_peaks')+scale_color_viridis()
dev.off()

pdf('nCount_peaks_11_GBM_snATAC_samples.pdf',width=8, height=6,useDingbats=FALSE)
FeaturePlot(atac,features='nCount_peaks')+scale_color_viridis()
dev.off()

pdf('Passed_filters_11_GBM_snATAC_samples.pdf',width=8, height=6,useDingbats=FALSE)
FeaturePlot(atac,features='passed_filters')+scale_color_viridis()
dev.off()




