system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

###For some samples (with many cells >6K) python-package used in chromVar doesn't work properly; Need to use this 2 commands:
###export OMP_NUM_THREADS=1
###export USE_SIMPLE_THREADED_LEVEL3=1
###export OPENBLAS_NUM_THREADS=1

library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)
library(RColorBrewer)




data_dir='/diskmnt/Projects/cptac_scratch_4/CPTAC3_GBM/snATAC_seq/signac_nadya/DATA/'

subt=read.table(paste(data_dir,'gbm_all_subtype_collections.v5.1.tsv',sep=''),sep='\t',header=TRUE)

ATAC=readRDS(paste(data_dir,'13_snATAC_GBM.chromvar.v20210108.rds',sep=''))

meta=read.table(paste(data_dir,'13_samples_Annotated.v.20210201.ManuallyReiewed.tsv',sep=''),
sep='\t',header=TRUE)

orig_1=as.data.frame(ATAC$dataset)
orig_1$i_barc=rownames(orig_1)

meta=meta[orig_1$i_barc,]

ATAC$cell_type_manual=meta$cell_type_manual
ATAC$Piece_ID=meta$Piece_ID


subt1=subt[subt$case %in% as.character(ATAC$Piece_ID),]
subt1=subt1 %>% select ('case','rna_wang_cancer_cell_2017','multiomic')
subt1=subt1[!is.na(subt1$rna_wang_cancer_cell_2017),]

###############################################
#1.Wang subtypes###############################
###############################################


#table(ATAC@meta.data %>% dplyr::select ('cell_type_manual','dataset'))

atac=ATAC

ATAC=subset(atac, cell_type_manual=="Tumor")

ATAC$test=as.character(ATAC$Piece_ID)


mesench=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="Mesenchymal"])
classic=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="Classical"])
pron=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="Proneural"])
idh=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="IDH mutant"])

ATAC$test=ifelse(ATAC$Piece_ID %in% mesench, "Mesenchymal", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% classic, "Classical", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% pron, "Proneural", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% idh, "IDH_mutant", ATAC$test)

atac=ATAC
ATAC=subset(atac, test %in% c('Mesenchymal','Classical','Proneural','IDH_mutant'))

DefaultAssay(ATAC) <- 'chromvar'
chromv= GetAssayData(object = ATAC)

cell_types=unique(ATAC$test)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]

ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

cell_types=unique(ATAC$test)

for (cell_t1 in cell_types){

    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types!=cell_t1]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=(all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2)
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]
        all_wilcoxon_stat$p_val_adj=p.adjust(as.numeric(as.character(unlist(all_wilcoxon_stat$V5))),
method="fdr")
colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
final_wilcoxon_stat=rbind(final_wilcoxon_stat, all_wilcoxon_stat)
print (cell_t1)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat=final_wilcoxon_stat[order(final_wilcoxon_stat$pvalue),]
tab=final_wilcoxon_stat
colnames(tab)[1]='cell_group_1'
colnames(tab)[3:4]=c('mean_score.group1','mean_score.group2')
tab$cell_t2="Other"

write.table(tab,paste("out/Score_difference_Wang_Subtypes_vs_Others.20210202.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)


###############################################
#1.Multiomic###################################
###############################################


#table(ATAC@meta.data %>% dplyr::select ('cell_type_manual','dataset'))

atac=ATAC

ATAC=subset(atac, cell_type_manual=="Tumor")

ATAC$test=as.character(ATAC$Piece_ID)


nmf1=as.character(subt1$case[subt1$multiomic=="nmf1"])
nmf2=as.character(subt1$case[subt1$multiomic=="nmf2"])
nmf3=as.character(subt1$case[subt1$multiomic=="nmf3"])
idh=as.character(subt1$case[subt1$multiomic=="IDH mutant"])

ATAC$test=ifelse(ATAC$Piece_ID %in% nmf1, "nmf1", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% nmf2, "nmf2", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% nmf3, "nmf3", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% idh, "IDH_mutant", ATAC$test)

atac=ATAC
ATAC=subset(atac, test %in% c('nmf1','nmf2','nmf3','IDH_mutant'))

DefaultAssay(ATAC) <- 'chromvar'
chromv= GetAssayData(object = ATAC)

cell_types=unique(ATAC$test)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]

ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

cell_types=unique(ATAC$test)

for (cell_t1 in cell_types){

    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types!=cell_t1]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=(all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2)
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]
        all_wilcoxon_stat$p_val_adj=p.adjust(as.numeric(as.character(unlist(all_wilcoxon_stat$V5))),
method="fdr")
colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
final_wilcoxon_stat=rbind(final_wilcoxon_stat, all_wilcoxon_stat)
print (cell_t1)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat=final_wilcoxon_stat[order(final_wilcoxon_stat$pvalue),]
tab=final_wilcoxon_stat
colnames(tab)[1]='cell_group_1'
colnames(tab)[3:4]=c('mean_score.group1','mean_score.group2')
tab$cell_t2="Other"

write.table(tab,paste("out/Score_difference_Multiomic_Subtypes_vs_Others.20210202.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)









####################################
###For_the_sample_level_comprison###
####################################

atac=ATAC

ATAC=subset(atac, cell_type_manual=="Tumor")

ATAC$test=as.character(ATAC$Piece_ID)

nmf1=as.character(subt1$case[subt1$multiomic=="nmf1"])
nmf2=as.character(subt1$case[subt1$multiomic=="nmf2"])
nmf3=as.character(subt1$case[subt1$multiomic=="nmf3"])
idh=as.character(subt1$case[subt1$multiomic=="IDH mutant"])

ATAC$test=ifelse(ATAC$Piece_ID %in% nmf1, "nmf1", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% nmf2, "nmf2", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% nmf3, "nmf3", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% idh, "IDH_mutant", ATAC$test)

atac=ATAC
ATAC=subset(atac, test %in% c('nmf1','nmf2','nmf3','IDH_mutant'))


ATAC$test=ATAC$Piece_ID
DefaultAssay(ATAC) <- 'chromvar'
chromv= GetAssayData(object = ATAC)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]

ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

cell_types=unique(as.character(ATAC$Piece_ID))
cell_types_1=nmf1
cell_types_2=cell_types[!(cell_types %in% cell_types_1)]

for (cell_t1 in cell_types_1){
    for (cell_t2 in cell_types_2){

    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t2]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value,cell_t2)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=(all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2)
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]
        all_wilcoxon_stat$p_val_adj=p.adjust(as.numeric(as.character(unlist(all_wilcoxon_stat$V5))),
method="fdr")
colnames(all_wilcoxon_stat)[1:2]=c('cell_t2','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
final_wilcoxon_stat=rbind(final_wilcoxon_stat, all_wilcoxon_stat)
print(cell_t2)
}
print (cell_t1)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat=final_wilcoxon_stat[order(final_wilcoxon_stat$pvalue),]
tab=final_wilcoxon_stat
colnames(tab)[1]='cell_group_1'
colnames(tab)[3:4]=c('mean_score.group1','mean_score.group2')

write.table(tab,paste("out/Score_difference_nmf3_vs_Others.20210202.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)









#########################################
###For_the_sample_level_comprison_Wang###
#########################################

atac=ATAC

ATAC=subset(atac, cell_type_manual=="Tumor")

ATAC$test=as.character(ATAC$Piece_ID)

mesench=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="Mesenchymal"])
classic=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="Classical"])
pron=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="Proneural"])
idh=as.character(subt1$case[subt1$rna_wang_cancer_cell_2017=="IDH mutant"])

ATAC$test=ifelse(ATAC$Piece_ID %in% mesench, "Mesenchymal", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% classic, "Classical", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% pron, "Proneural", ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% idh, "IDH_mutant", ATAC$test)

atac=ATAC
ATAC=subset(atac, test %in% c('Mesenchymal','Classical','Proneural','IDH_mutant'))


ATAC$test=ATAC$Piece_ID
DefaultAssay(ATAC) <- 'chromvar'
chromv= GetAssayData(object = ATAC)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]

ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

cell_types=unique(as.character(ATAC$Piece_ID))
cell_types_1=mesench
cell_types_2=cell_types[!(cell_types %in% cell_types_1)]

for (cell_t1 in cell_types_1){
    for (cell_t2 in cell_types_2){

    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t2]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value,cell_t2)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=(all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2)
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]
        all_wilcoxon_stat$p_val_adj=p.adjust(as.numeric(as.character(unlist(all_wilcoxon_stat$V5))),
method="fdr")
colnames(all_wilcoxon_stat)[1:2]=c('cell_t2','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
final_wilcoxon_stat=rbind(final_wilcoxon_stat, all_wilcoxon_stat)
print(cell_t2)
}
print (cell_t1)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat=final_wilcoxon_stat[order(final_wilcoxon_stat$pvalue),]
tab=final_wilcoxon_stat
colnames(tab)[1]='cell_group_1'
colnames(tab)[3:4]=c('mean_score.group1','mean_score.group2')

write.table(tab,paste("out/Score_difference_mesench_vs_Others.20210202.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)


