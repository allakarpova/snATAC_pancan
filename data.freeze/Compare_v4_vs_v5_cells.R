library(tidyverse)
library(data.table)
library(Signac)
library(Seurat)
library(doParallel)

m <- '/diskmnt/Projects/snATAC_analysis/immune/misc/Old_new_objects_location.txt'
locations <- fread(m)
locations.few <- locations[1:3,]

registerDoParallel(cores=10)
stats <- foreach::foreach(new=locations$new, old=locations$old) %dopar% {
  cat('open obj\n')
  new.obj <- readRDS(list.files(path = new, full.names = T, pattern = '*rds'))
  print(old)
  old.obj <- readRDS(list.files(path = old, full.names = T, pattern = '*rds'))
  cells.new <- colnames(new.obj)
  cells.old <- colnames(old.obj)
  cat('working with cells\n')
  common.cells <- intersect(cells.new, cells.old)
  in.new.only <- setdiff(cells.new, cells.old)
  in.old.only <- setdiff(cells.old,cells.new)
  if (sum(grepl('RNA', colnames(new.obj@meta.data))>0 )) {
    not.cells <- rownames(subset(new.obj@meta.data, is_cell!=1))
  } else {
    not.cells <- rownames(subset(new.obj@meta.data, is__cell_barcode!=1))
  }
  in.new.only.not.cells <- intersect(in.new.only, not.cells)
  to.return <- c('new' = new,
                 'common.cells' = length(common.cells),
                 'in.new.only' = length(in.new.only),
                 'in.old.only' = length(in.old.only), 
                 'not.cells' = length(not.cells),
                 'in.new.only.not.cells'  = length(in.new.only.not.cells))
  print(to.return)
  return(to.return)
}

stopImplicitCluster()

stats <- do.call('rbind', stats)

fwrite(cbind(locations, stats), '/diskmnt/Projects/snATAC_analysis/immune/misc/New_objects_location_with_stats.txt', sep = '\t')

raw.obj <- readRDS('~/lab_Ding/work/single_cell/snATAC_pancan/test_raw_filtered/UCEC_TWHG-CPT4096DU-XBa1_1_processed_atac.rds')
filt.obj <- readRDS('~/lab_Ding/work/single_cell/snATAC_pancan/test_raw_filtered/UCEC_TWHG-CPT4096DU-XBa1_1_processed_atac_filtered.rds')

head(raw.obj@meta.data)
ggplot() +
  geom_density(data=raw.obj@meta.data, aes(x = log10(passed_filters)), color = 'red') +
  geom_density(data=filt.obj@meta.data, aes(x = log10(passed_filters)), color = 'black')

  

ggplot() +
  
  geom_jitter(data=raw.obj@meta.data, aes(x = log10(passed_filters), y = frip_500MACS2), color = 'red', size = 0.5) +
  geom_jitter(data=filt.obj@meta.data, aes(x = log10(passed_filters), y = frip_500MACS2), color = 'black', size = 0.1)

common.cells <- intersect(rownames(raw.obj@meta.data), rownames(filt.obj@meta.data))
toplot <- cbind(raw.obj@meta.data[common.cells, 'peak_region_fragments'], filt.obj@meta.data[common.cells, 'peak_RF_500MACS2'])
colnames(toplot) <- c('Raw', 'Filtered')
toplot <- data.frame(toplot)
ggplot(data=toplot, aes(x = Raw, y = Filtered)) +
  geom_point(color = 'black') 
  

toplot <- cbind(raw.obj@meta.data['peak_region_fragments'], raw.obj@meta.data['peak_RF_500MACS2'])
colnames(toplot) <- c('Cellranger_peak_region_fragments', 'MACS2_peak_RF_500MACS2')
toplot <- data.frame(toplot)
ggplot(data=toplot, aes(x = Cellranger_peak_region_fragments, y = MACS2_peak_RF_500MACS2)) +
  geom_point(color = 'black')

ggplot() +
  geom_jitter(data=raw.obj@meta.data, aes(x = peak_region_fragments, y = peak_RF_500MACS2), color = 'red', size = 0.5) +
  geom_jitter(data=filt.obj@meta.data, aes(x = peak_region_fragments, y = peak_RF_500MACS2), color = 'black', size = 0.1) +
  xlim(1,20000) + ylim(1,20000)

  
ggplot(data=raw.obj@meta.data, aes(x = TSS.enrichment, fill = high.tss)) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~high.tss)


ggplot() +
  geom_jitter(data=raw.obj@meta.data, aes(x = peak_region_fragments, y = TSS.enrichment), color = 'red', size = 0.5) +
  geom_jitter(data=filt.obj@meta.data, aes(x = peak_region_fragments, y = TSS.enrichment), color = 'black', size = 0.1)

ggplot() +
  geom_jitter(data=raw.obj@meta.data, aes(x = peak_RF_500MACS2, y = TSS.enrichment), color = 'red', size = 0.5) +
  geom_jitter(data=filt.obj@meta.data, aes(x = peak_RF_500MACS2, y = TSS.enrichment), color = 'black', size = 0.1)


VlnPlot(raw.obj, features = 'peak_RF_500MACS2') +
  VlnPlot(filt.obj, features = 'peak_RF_500MACS2')

VlnPlot(raw.obj, features = 'passed_filters') +
  VlnPlot(filt.obj, features = 'passed_filters')


DimPlot(raw.obj, label = T) /
  DimPlot(filt.obj, label = T)

toplot <- cbind(Embeddings(filt.obj, reduction = 'umap'), 
                seurat_clusters = raw.obj@meta.data[rownames(filt.obj@meta.data), 'seurat_clusters']) %>%
  data.frame()
toplot$seurat_clusters <- factor(toplot$seurat_clusters, levels = as.character(0:19))

ggplot(data = data.frame(toplot), aes (x=UMAP_1, y=UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.5) +
  theme_cowplot()

toplot <- cbind(Barcodes = rownames(raw.obj@meta.data), 
                peak_region_fragments = as.numeric(raw.obj$peak_region_fragments))  %>%
  data.frame() %>%
  mutate(peak_region_fragments=as.numeric(peak_region_fragments)) %>%
  arrange(-peak_region_fragments)
toplot$Barcodes <- factor(toplot$Barcodes, levels = toplot$Barcodes)

ggplot() +
  geom_jitter(data=toplot, aes(x = Barcodes, y = peak_region_fragments), color = 'red', size = 0.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

toplot <- cbind(Barcodes = rownames(filt.obj@meta.data), 
                  peak_region_fragments = as.numeric(filt.obj$peak_RF_500MACS2))  %>%
    data.frame() %>%
    mutate(peak_region_fragments=as.numeric(peak_region_fragments)) %>%
    arrange(-peak_region_fragments)
toplot$Barcodes <- factor(toplot$Barcodes, levels = toplot$Barcodes)
  
ggplot(data=toplot, aes(x = Barcodes, y = peak_region_fragments)) +
    geom_jitter(color = 'black', size = 0.5) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 











