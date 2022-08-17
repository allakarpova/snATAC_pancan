library(spatstat)
library(Signac)
library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(paletteer)
library(readxl)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


panc <- readRDS ('~/lab_Ding/work/single_cell/senescence/snATAC_RNA/objects/SM001_merged_both_combo_hepatocytes_v1.rds')
panc <- NormalizeData(panc)
DefaultAssay(panc) <- 'peaksinters'
cell.type <- 'Hepatocytes'

# annotate peaks
all_m <- data.frame (panc@assays$peaksinters@ranges)
all_m$m_coords=paste(all_m$seqnames,":",all_m$start,"-",all_m$end,sep='')
m_all_m <- StringToGRanges(all_m$m_coords, sep = c(":", "-"))
###Annotate the Motifs with the closest gene, to get the info about promoter regions:
peakAnno <- annotatePeak(m_all_m, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# annotate object
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(panc) <- annotations	


DimPlot(panc)
#saveRDS(panc, '~/lab_Ding/work/single_cell/senescence/snATAC_combo/objects/mouse_liver/2_snATAC_Merged_mouse_liver_young_old_chromvar_macs2_chromvar_annot.rds')

wd <- paste0('~/lab_Ding/work/single_cell/senescence/snATAC_RNA/DAMs/Hepatocytes')
dir.create(wd, recursive = T)
setwd(wd)
add_filename <- 'mouse_liver'


Idents(panc) <- 'Cell_type_macs2'

for (ct in unique(panc$Cell_type_macs2)) {
  da_peaks <- FindMarkers(
    object = panc,
    subset.ident = ct,
    group.by = 'Age',
    ident.1 = 'Old',
    #ident.2 = 'Sst',
    only.pos = F,
    test.use = 'LR',
    latent.vars = 'nCount_peaks'
  )
  fwrite(da_peaks, paste0('Marker_macs2_peaks_old_vs_young_', ct, '.txt'), sep = '\t', row.names = T)
  
  VlnPlot(
    object = panc,ncol = 5,
    features = rownames(da_peaks)[1:10],
    pt.size = 0.1,group.by = 'Age', 
    idents = ct
  )
  ggsave(paste0('Vlnplot_top10_differential_peaks_', ct,'.pdf'), useDingbats = F, width = 20, height = 8)
  
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0, ])
  enriched.motifs <- FindMotifs(
    object = panc,
    features = top.da.peak
  )
  fwrite(enriched.motifs, paste0('Enriched_motifs_in_old_vs_young_', ct, '.txt'), sep = '\t', row.names = T)
  
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$avg_log2FC < 0, ])
  enriched.motifs <- FindMotifs(
    object = panc,
    features = top.da.peak
  )
  fwrite(enriched.motifs, paste0('Enriched_motifs_in_young_vs_old_', ct, '.txt'), sep = '\t', row.names = T)
}

# use peak annotation to see which genes have differentially accessible peaks
peak.annotation <- data.frame (peakAnno)
peak.annotation$coord <- paste0(peak.annotation$seqnames,"-",peak.annotation$start,"-",peak.annotation$end)

for (ct in unique(panc$Cell_type_macs2)) {
  da_peaks <- fread(input = paste0('Marker_macs2_peaks_old_vs_young_', ct,'.txt'), data.table = F)
  da_peaks <- merge(subset(da_peaks, p_val_adj < 0.005) , peak.annotation[c('SYMBOL', 'annotation' ,'GENENAME', 'coord')], by.x = 'V1', by.y = 'coord', all.x = T)
}


markers.from.paper <- list (Sen.core = sen.core.mouse, Sen.effector = sen.effector.mouse, sasp = sasp.mouse)

# annotate DA peaks
ct <- 'Hepathocytes'
da_peaks <- fread(paste0('~/lab_Ding/work/single_cell/senescence/snATAC_combo/DAMs/Marker_macs2_peaks_old_vs_young_', ct, '.txt'), data.table = F)
peak.an.df <- data.frame(peakAnno@anno)
peak.an.df$m_coords=paste(peak.an.df$seqnames,"-",peak.an.df$start,"-",peak.an.df$end,sep='')
da_peaks <- merge(subset(da_peaks, p_val_adj < 0.005) , peak.an.df[c('SYMBOL', 'annotation' ,'GENENAME', 'm_coords')], by.x = 'V1', by.y = 'm_coords', all.x = T)

da_peaks %>% 
  filter (avg_log2FC > 0 & p_val_adj < 0.05) %>% 
  filter( grepl('Intron', annotation)) %>% dplyr::select (SYMBOL) %>% unlist %>% as.character


## find DEGs between hepatocyte clusters
deg <- FindAllMarkers(panc, assay = 'RNA', test.use = 'LR', min.pct = 0.05, logfc.threshold = 0.3)
fwrite(deg , paste0('DEGs_clusters', cell.type, '.txt'), sep = '\t')

deg <- fread(paste0('DEGs_clusters', cell.type, '.txt'), data.table = F)
top10.deg <- deg %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, -p_val_adj)
top10.deg %>% dplyr::group_by(cluster) %>% tally()
fwrite(top10.deg, 'Top10_deg_clusters_hepato.txt', sep = '\t')

# first compute the GC content for each peak
panc <- RegionStats(panc, genome = BSgenome.Mmusculus.UCSC.mm10)
panc <- LinkPeaks(panc, peak.assay = 'peaksinters', expression.assay = 'SCT', genes.use = top10.deg$gene)
data.frame(panc@assays$peaksinters@links)$gene %>% unique %>% length 

linked.peaks <- GetLinkedPeaks(panc, top10.deg$gene, min.abs.score = 0.2)
fwrite (data.frame(linked.peaks), 'Linked_peaks_for_cluster_deg.txt',sep = '\t')
linked.genes <- GetLinkedGenes(panc, linked.peaks, min.abs.score = 0.2)

Idents(panc)
for (gene in (data.frame(panc@assays$peaksinters@links)$gene %>% unique)) {
  print(gene)
  CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'wknn_res.0.5', region =gene, features = gene, expression.assay = "SCT", extend.upstream = 2000, extend.downstream = 1000)
  ggsave(paste0('CoveragePlot_sen.core_', gene, '_', cell.type, '_linked.pdf'), useDingbats = F, width = 10, height = 7)
}

panc <- LinkPeaks(panc, peak.assay = 'peaksinters', expression.assay = 'RNA', genes.use = sasp.mouse)

panc$sasp.status <- ifelse (panc$sasp > 0.1, 'SASP.high', 'SASP.low')
for (gene in (data.frame(panc@assays$peaksinters@links)$gene %>% unique)) {
  print(gene)
  CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'sasp.status', region =gene, features = gene, expression.assay = "RNA", extend.upstream = 2000, extend.downstream = 1000)
  ggsave(paste0('CoveragePlot_sasp_', gene, '_RNA_', cell.type, '_sasp.status_linked.pdf'), useDingbats = F, width = 10, height = 5)
}

for (gene in (data.frame(panc@assays$peaksinters@links)$gene %>% unique)) {
  print(gene)
  CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'age', region =gene, features = gene, expression.assay = "RNA", extend.upstream = 2000, extend.downstream = 1000)
  ggsave(paste0('CoveragePlot_sasp_', gene, '_RNA_', cell.type, '_age_linked.pdf'), useDingbats = F, width = 10, height = 5)
}

panc@meta.data[c('sasp.status', 'age')]

panc$sasp.status <- ifelse (panc$sasp > 0.05, 'SASP.high', 'SASP.low')
xtabs(~sasp.status + age, data=panc@meta.data[c('sasp.status', 'age')])
fisher.test(xtabs(~sasp.status + age, data=panc@meta.data[c('sasp.status', 'age')]))

VlnPlot(panc, features = 'Il18', assay = 'RNA', group.by = 'age', cols = c('#a6cee3', '#1f78b4'), log = F) +
  ggpubr::stat_compare_means()


FeatureScatter(panc,  feature1 = 'Cdkn1a', feature2 = 'Cdkn2a', group.by = 'age')



