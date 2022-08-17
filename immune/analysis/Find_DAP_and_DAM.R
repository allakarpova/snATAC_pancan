suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
#suppressMessages(library(org.Hs.eg.db))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

#find differentially accessible peaks
panc <- readRDS('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/objects/v1.2/Reclustered_immune_snATAC_Merged_94_sample_obj.v1.2_new_macs2.chromvar.rds')
meta <- fread('~/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/cell_typing/v1.2/Reclustered_immune_snATAC_Merged_94_sample_obj.v1.2_new_macs2_upd.metaData', data.table = F) %>%
  data.frame(row.names = 1, check.rows = F, check.names = F)
panc <- AddMetaData(panc, meta)

wd <- '/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/immune_analysis/DAM/v1.2'
dir.create(wd, recursive = T)
setwd(wd)
add_filename <- 'snATAC_Merged_94_sample_obj.v1.2_new_macs2'


# annotate peaks
all_m <- data.frame (panc@assays$X500peaksMACS2@ranges)
all_m$m_coords=paste0(all_m$seqnames,":",all_m$start,"-",all_m$end)
m_all_m <- StringToGRanges(all_m$m_coords, sep = c(":", "-"))
###Annotate the Motifs with the closest gene, to get the info about promoter regions:
peakAnno <- annotatePeak(m_all_m, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
peak.annotation <- data.frame (peakAnno)
peak.annotation$coord <- paste0(peak.annotation$seqnames,"-",peak.annotation$start,"-",peak.annotation$end)

Annotation(panc)

Idents(panc) <- 'cell_type_upd'
unique(Idents(panc))

da_peaks <- FindMarkers(
  object = panc,
  ident.1 = 'Exhausted CD8 T-cells',
  ident.2 = 'CD8 CTL',
  only.pos = F,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
da_peaks <- merge(subset(da_peaks, p_val_adj < 0.005) , peak.annotation[c('SYMBOL', 'annotation' ,'GENENAME', 'coord')], by.x = 0, by.y = 'coord', all.x = T)

fwrite(da_peaks, paste0('Marker_macs2_peaks_exhausted_vs_CTL.txt'), sep = '\t', row.names = T)

DoHeatmap(panc, slot = 'data', 
          features = subset(da_peaks, p_val_adj < 0.005)$Row.names, cells = WhichCells(panc, idents = c('Exhausted CD8 T-cells', 'CD8 CTL')) , )

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0, ])
enriched.motifs <- FindMotifs(
  object = panc,
  features = top.da.peak
)
fwrite(enriched.motifs, paste0('Enriched_motifs_in_old_vs_young_', ct, '.txt'), sep = '\t', row.names = T)

DefaultAssay(panc) <- 'chromvar'
dif.active.motifs <-  FindMarkers(
  object = panc,
  ident.1 = 'Exhausted CD8 T-cells',
  ident.2 = 'CD8 CTL',
  only.pos = T,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

names <- ConvertMotifID( object = panc,assay = 'chromvar', id = rownames(dif.active.motifs))
fwrite(dif.active.motifs, paste0('Dif.active.motifs_exhausted_vs_CTL.txt'), sep = '\t', row.names = T)


rownames(names)
