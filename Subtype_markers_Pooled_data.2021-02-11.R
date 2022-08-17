###Author: Nadezhda V. Terekhanova
###Date:2021/02/11
##################################

library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/GBM/Analysis/7.Subtype_markers/20210128_Subtypes/20210202/')

score=read.table('out/Score_difference_Wang_Subtypes_vs_Others.20210202.tsv',header=TRUE,sep="\t")
score=score[,1:3]
colnames(score)=c('cell_type','TF_Name','mean_score')



#Plotting all scores in a heatmap
score_diff=read.table('out/Score_difference_Wang_Subtypes_vs_Others.20210202.tsv',sep='\t',header=TRUE)
colnames(score_diff)[1:4]=c('cell_t2','TF_Name','mean_score1','mean_score2')
score_diff$diff=score_diff$mean_score1-score_diff$mean_score2
score_diff=score_diff[order(-score_diff$mean_score1),]
score_diff=score_diff[score_diff$p_val_adj<0.05,]

cell_types=unique(as.character(score_diff$cell_t2))

score_diff_1=score_diff[score_diff$diff>0 & score_diff$mean_score1>0,]
score_diff_1=dcast(score_diff_1, TF_Name~cell_t2,value.var='diff')
score_diff_1$na=apply(score_diff_1,1,function(x) sum(is.na(x)))
score_diff_1=score_diff_1[score_diff_1$na==3,]
score_diff_1$max=rowMeans(score_diff_1[,2:5],na.rm=TRUE)
score_diff_1=score_diff_1[order(-score_diff_1$max),]

score_diff_2=melt(score_diff_1[,1:5], by=list(score_diff_1$TF_Name))
score_diff_2=score_diff_2[!is.na(score_diff_2$value),]
score_diff_2=score_diff_2[order(-score_diff_2$value),]

mesench=as.character(score_diff$TF_Name[score_diff$cell_t2=="Mesenchymal"])
pron=as.character(score_diff$TF_Name[score_diff$cell_t2=="Proneural"])
class=as.character(score_diff$TF_Name[score_diff$cell_t2=="Classical"])
idh=as.character(score_diff$TF_Name[score_diff$cell_t2=="IDH_mutant"])


all_top=rbind(score_diff_2[score_diff_2$variable=="Mesenchymal",][1:20,], score_diff_2[score_diff_2$variable=="Proneural",][1:20,], score_diff_2[score_diff_2$variable=="Classical",][1:20,], score_diff_2[score_diff_2$variable=="IDH_mutant",][1:20,])
all_top=all_top[!is.na(all_top$TF_Name),]


#####Now grabbing info for the selected TFs for all subtypes
markers=score[score$TF_Name %in% all_top$TF_Name,]
enr=all_top[,1:2]
colnames(enr)[2]="Cell_type_enriched"
markers=merge(markers,enr,all.x=TRUE)



cols <- brewer.pal(9, "YlOrRd")
getPalette= colorRampPalette(cols)

all_merged=markers
all_merged$Cell_type_enriched=as.character(unlist(all_merged$Cell_type_enriched))
all_merged$cell_type=as.character(unlist(all_merged$cell_type))


all_merged$cell_type=factor(all_merged$cell_type,levels=c('Mesenchymal', 'Proneural', 'Classical', 'IDH_mutant'))
all_merged=all_merged[!is.na(all_merged$cell_type),]
all_merged=all_merged[!is.na(all_merged$Cell_type_enriched),]
all_merged$Cell_type_enriched=factor(all_merged$Cell_type_enriched,levels=c('Mesenchymal', 'Proneural', 'Classical', 'IDH_mutant'))


####Heatmap-version:
all_merged$z_mean_score=scale(all_merged$mean_score)
cols <- c("#377EB8","#377EB8","white","#E41A1C","#E41A1C")
getPalette= colorRampPalette(cols)

print(min(all_merged$mean_score))
print(max(all_merged$mean_score))

#all_merged$mean_score=ifelse(all_merged$mean_score>4,4,all_merged$mean_score)
#all_merged$mean_score=ifelse(-all_merged$mean_score>4,-4,all_merged$mean_score)

p = ggplot()
p = p + facet_grid(.~Cell_type_enriched,drop=T,scales = "free", space = "free")
p = p + geom_tile(data=all_merged,aes(x=TF_Name, y=cell_type, fill= mean_score), linetype="blank",width=0.9, height=0.9) + scale_fill_gradientn(name= "Motif Score", colours=getPalette(100), na.value="grey", limit=c(-2.2,2.2))
#p = p + geom_tile(data=all_merged,aes(x=TF_Name, y=cell_type, fill= mean_score), linetype="blank",width=0.9, height=0.9) + scale_fill_gradientn(name= "Motif Score", colours=getPalette(100), na.value="grey")
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5,hjust=0.95), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())+
  theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5,hjust=0.95), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
#p=p+theme(strip.background = element_blank(), strip.text.x=element_text(size=12)
p=p+theme(strip.text.x=element_text(size=12)
)


pdf(paste("Wang_subtypes.bulk.2021-02-02.pdf",sep=""),width=14, height=2.9,useDingbats=FALSE)
p
dev.off()


################################################
#####Now checking multiomic-subtype-markers:####
################################################


score=read.table('out/Score_difference_Multiomic_Subtypes_vs_Others.20210202.tsv',header=TRUE,sep="\t")
score=score[,1:3]
colnames(score)=c('cell_type','TF_Name','mean_score')

#Plotting all scores in a heatmap
score_diff=read.table('out/Score_difference_Multiomic_Subtypes_vs_Others.20210202.tsv',sep='\t',header=TRUE)
colnames(score_diff)[1:4]=c('cell_t2','TF_Name','mean_score1','mean_score2')
score_diff$diff=score_diff$mean_score1-score_diff$mean_score2
score_diff=score_diff[order(-score_diff$mean_score1),]
score_diff=score_diff[score_diff$p_val_adj<0.05,]

cell_types=unique(as.character(score_diff$cell_t2))

score_diff_1=score_diff[score_diff$diff>0 & score_diff$mean_score1>0,]
score_diff_1=dcast(score_diff_1, TF_Name~cell_t2,value.var='diff')
score_diff_1$na=apply(score_diff_1,1,function(x) sum(is.na(x)))
score_diff_1=score_diff_1[score_diff_1$na==3,]
score_diff_1$max=rowMeans(score_diff_1[,2:5],na.rm=TRUE)
score_diff_1=score_diff_1[order(-score_diff_1$max),]

score_diff_2=melt(score_diff_1[,1:5], by=list(score_diff_1$TF_Name))
score_diff_2=score_diff_2[!is.na(score_diff_2$value),]
score_diff_2=score_diff_2[order(-score_diff_2$value),]


all_top=rbind(score_diff_2[score_diff_2$variable=="nmf1",][1:20,], score_diff_2[score_diff_2$variable=="nmf2",][1:20,], score_diff_2[score_diff_2$variable=="nmf3",][1:20,], score_diff_2[score_diff_2$variable=="IDH_mutant",][1:20,])
all_top=all_top[!is.na(all_top$TF_Name),]


#####Now grabbing info for the selected TFs for all subtypes
markers=score[score$TF_Name %in% all_top$TF_Name,]
enr=all_top[,1:2]
colnames(enr)[2]="Cell_type_enriched"
markers=merge(markers,enr,all.x=TRUE)



cols <- brewer.pal(9, "YlOrRd")
getPalette= colorRampPalette(cols)

all_merged=markers
all_merged$Cell_type_enriched=as.character(unlist(all_merged$Cell_type_enriched))
all_merged$cell_type=as.character(unlist(all_merged$cell_type))


all_merged$cell_type=factor(all_merged$cell_type,levels=c('nmf1', 'nmf2', 'nmf3', 'IDH_mutant'))
all_merged=all_merged[!is.na(all_merged$cell_type),]
all_merged=all_merged[!is.na(all_merged$Cell_type_enriched),]
all_merged$Cell_type_enriched=factor(all_merged$Cell_type_enriched,levels=c('nmf1', 'nmf2', 'nmf3', 'IDH_mutant'))


####Heatmap-version:
all_merged$z_mean_score=scale(all_merged$mean_score)
cols <- c("#377EB8","#377EB8","white","#E41A1C","#E41A1C")
getPalette= colorRampPalette(cols)

print(min(all_merged$mean_score))
print(max(all_merged$mean_score))

#all_merged$mean_score=ifelse(all_merged$mean_score>4,4,all_merged$mean_score)
#all_merged$mean_score=ifelse(-all_merged$mean_score>4,-4,all_merged$mean_score)

p = ggplot()
p = p + facet_grid(.~Cell_type_enriched,drop=T,scales = "free", space = "free")
p = p + geom_tile(data=all_merged,aes(x=TF_Name, y=cell_type, fill= mean_score), linetype="blank",width=0.9, height=0.9) + scale_fill_gradientn(name= "Motif Score", colours=getPalette(100), na.value="grey", limit=c(-1.7,1.7))
#p = p + geom_tile(data=all_merged,aes(x=TF_Name, y=cell_type, fill= mean_score), linetype="blank",width=0.9, height=0.9) + scale_fill_gradientn(name= "Motif Score", colours=getPalette(100), na.value="grey")
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5,hjust=0.95), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())+
  theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5,hjust=0.95), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
#p=p+theme(strip.background = element_blank(), strip.text.x=element_text(size=12)
p=p+theme(strip.text.x=element_text(size=12)
)


pdf(paste("Mutiomic_subtypes.bulk.2021-02-02.pdf",sep=""),width=14, height=2.9,useDingbats=FALSE)
p
dev.off()
