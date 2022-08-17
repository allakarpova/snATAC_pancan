library(Signac)
library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(paletteer)
library(readxl)
library(ggplot2) 
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(cowplot)


panc <- readRDS ('~/lab_Ding/work/single_cell/senescence/snATAC_combo/objects/mouse_liver/2_snATAC_Merged_mouse_liver_young_old_chromvar_macs2_chromvar_annot.rds')
DefaultAssay(panc) <- 'peaksMACS2'


all_m <- data.frame (panc@assays$peaksMACS2@ranges)
all_m$m_coords=paste(all_m$seqnames,":",all_m$start,"-",all_m$end,sep='')
m_all_m <- StringToGRanges(all_m$m_coords, sep = c(":", "-"))
###Annotate the Motifs with the closest gene, to get the info about promoter regions:
peakAnno <- annotatePeak(m_all_m, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")



wd <- '~/lab_Ding/work/single_cell/senescence/snRNA_combo/module_signatures/all.cells'
dir.create(wd, recursive = T)
setwd(wd)
cell.type = 'mouse_liver'

source('~/R_working_dir/scripts/senescence/markers.R')
markers.from.paper <- list (Sen.core = sen.core.mouse, Sen.effector = sen.effector.mouse, sasp = sasp.mouse)

total.sen <- markers.from.paper

cat('calculate modules\n')
panc <- AddModuleScore(
  object = panc,
  assay = 'RNA',
  features = total.sen,
  ctrl = 10,
  name = 'senescence'
)
cat('done\n')

n.list <- length(total.sen)
colnames(panc@meta.data)[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)] <- names(total.sen)


for(s in names(total.sen)) {
  #s <- 'Sen.core'
  VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'Age',
          assay = 'RNA', cols =c( '#fc8d62', '#66c2a5'))
  ggsave(paste0( 'Vlnplot_', cell.type, '_scores_', s, '.pdf'), useDingbats = F, width = 18, height = 8)
  
  VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'seurat_clusters', ncol = 1, split.by = 'Age',
          assay = 'RNA',cols =c( '#fc8d62', '#66c2a5'))
  ggsave(paste0( 'Vlnplot_', cell.type, '_clusters_scores_', s, '.pdf'), useDingbats = F, width = 18, height = 8)
  
  FeaturePlot(panc, features = s, order = T, min.cutoff = 0)
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '.pdf'), useDingbats = F, width = 6, height = 5)
  
  FeaturePlot(panc, features = s, order = T, min.cutoff = 0,split.by = 'Age')
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '_splitted.pdf'), useDingbats = F, width = 18, height = 5.8)
  
  FeaturePlot(panc, features = total.sen[[s]], order = T, min.cutoff = 0,split.by = 'Age')
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '_gene_level_splitted.pdf'), useDingbats = F, width = 10, height = 3*length(total.sen[[s]]), limitsize = F)
  
  ggplot (panc@meta.data, aes_string (x = s, fill = 'Age')) +
    geom_density(alpha = 0.5) +
    theme_cowplot() +
    facet_wrap(~Cell_type, scales = 'free') +
    ggtitle(s) +
    scale_fill_manual(values =c( '#fc8d62', '#66c2a5'))
  ggsave(paste0('Densityplot_',cell.type, '_scores_', s, '_filled.pdf'), useDingbats = F, width = 8, height = 5.8)
}


saveRDS(panc, '~/lab_Ding/work/single_cell/senescence/snRNA_combo/objects/mouse_liver/2_snRNA_Merged_mouse_liver_young_old_annot.rds')
