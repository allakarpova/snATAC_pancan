# Integrate with public PBMC data
###libraries
##################
library(future)

plan("multicore", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
################################

###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

###################################

runAllNormalization <- function(obj, dims) {
  #### run normalization to get initial clusters ###
  ########
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 20) %>% 
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400 ) %>% 
    FindNeighbors(
      reduction = 'lsi',
      dims = 2:dims ) %>% 
    FindClusters(
      algorithm = 3,
      resolution = 1,
      verbose = FALSE
    ) %>% 
    RunUMAP(dims = 2:dims,
            reduction = 'lsi')
  return(obj)
}

############################################

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
annotations <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')


dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
int <- readRDS(input.path)
print(head(rownames(int)))

cat('opening pbmc object \n')
pbmc.obj <- readRDS('/diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v4.0/PBMC_granulocyte_sorted_10k/PBMC_granulocyte_sorted_10k_processed_multiomic.rds')

pbmc.meta <- fread('/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v5.0_for_immune/PBMC_1_samples_metadata_data_freeze_v3.0.tsv', header = T) %>% 
  data.frame(row.names = 1)
pbmc.meta <- pbmc.meta %>% 
  dplyr::select(!starts_with("prediction"))
pbmc.meta$Chemistry <- 'Multiome'
pbmc.meta$cell_type.harmonized.cancer <- pbmc.meta$cell_type.harmonized

pbmc.obj <- RenameCells(object = pbmc.obj, new.names = rownames(pbmc.meta))
pbmc.obj <- AddMetaData(pbmc.obj, pbmc.meta)

print(head(pbmc.obj@meta.data))

cat('pbmc matrix on my peaks \n')
peak.number <- nrow(int)
n.peaks <- round(peak.number/30)

frag <- Fragments(pbmc.obj@assays[['X500peaksMACS2']])
print(frag)

matrix.counts <- FeatureMatrix(
  fragments = frag,
  features = StringToGRanges(rownames(int)),
  process_n = n.peaks,
  sep = c("-","-"),
  cells = colnames(pbmc.obj)
)

pbmc.obj[['ATAC_immune']] <- CreateChromatinAssay(counts = matrix.counts,
                                                  fragments=frag,
                                                  annotation = annotations)
DefaultAssay(pbmc.obj)<-'ATAC_immune'
pbmc.obj <- DietSeurat(pbmc.obj, assays = 'ATAC_immune')

#################
cat('opening merged Satpathy object \n')
satpathy.obj <- readRDS('/diskmnt/Projects/snATAC_primary/05_merged_rds/external_Satpathy/v2/40_snATAC_Merged_new_peaks_normalized_external_Satpathy_BCC_healthy_immune_celltyped.rds')
satpathy.obj$Chemistry <- 'snRNA'
satpathy.obj$cell_type.harmonized.cancer <- satpathy.obj$Cell_type
satpathy.obj <- subset(satpathy.obj, Cancer == 'BCC', invert = TRUE)
print(head(satpathy.obj@meta.data))

cat('Satpathy matrix on my peaks \n')
peak.number <- nrow(int)
n.peaks <- round(peak.number/30)

frag <- Fragments(satpathy.obj@assays[['X500peaksMACS2']])
print(frag)

matrix.counts <- FeatureMatrix(
  fragments = frag,
  features = StringToGRanges(rownames(int)),
  process_n = n.peaks,
  sep = c("-","-"),
  cells = colnames(satpathy.obj)
)

satpathy.obj[['ATAC_immune']] <- CreateChromatinAssay(counts = matrix.counts,
                                                  fragments=frag,
                                                  annotation = annotations)
DefaultAssay(satpathy.obj)<-'ATAC_immune'
satpathy.obj <- DietSeurat(satpathy.obj, assays = 'ATAC_immune')


cat('Merging \n')
merged <- merge(x = int, y = list(pbmc.obj, satpathy.obj))


#normalize original object
merged <- merged %>% FindTopFeatures(min.cutoff = 10) %>%
  RunTFIDF() %>%
  RunSVD() %>%
  RunUMAP(reduction = "lsi", dims = 2:50)

DimPlot(merged, group.by = "Cancer", cols = 'Spectral') + 
  DimPlot(merged, group.by = "cell_type.harmonized.cancer") +
  DimPlot(merged, group.by = "Chemistry")
ggsave(paste0(add_filename, '_merged_dimplots.pdf'), width = 17, height = 4.5)

saveRDS(merged, paste0(add_filename, '_merged.rds'))


##########################
##### Integration #######
##########################

atac.split <- SplitObject(merged, split.by = 'Chemistry')

atac.split <- map(atac.split, function(obj) {
  obj <- FindTopFeatures(obj, min.cutoff = 10) %>%
    RunTFIDF() %>%
    RunSVD()
  return(obj)
})

#######integration############
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100 * 1024^3)
integration.anchors <- FindIntegrationAnchors(
  object.list = atac.split,
  anchor.features = rownames(merged),
  reduction = "rlsi",
  dims = 2:50
)


# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)


# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, 
                      reduction = "integrated_lsi", 
                      dims = 2:50)

p2 <- DimPlot(integrated, group.by = "Chemistry")
p1 <- DimPlot(merged, group.by = "Chemistry")
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_Chemistry.pdf'), width = 13, height = 4.5)
#saveRDS(integrated, paste0(add_filename, '_Chemistry.rds'))

p2 <- DimPlot(integrated, group.by = "cell_type.harmonized.cancer")
p1 <- DimPlot(merged, group.by = "cell_type.harmonized.cancer")
(p1 + ggtitle("Merged") + NoLegend()) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_cell_type.harm.cancer.pdf'), width = 15, height = 4.5)

p2 <- DimPlot(integrated, group.by = "Cancer", cols = 'Spectral')
p1 <- DimPlot(merged, group.by = "Cancer", cols = 'Spectral')
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(paste0(add_filename, '_Cancer.pdf'), width = 12, height = 4.5)


integrated  <-  integrated %>% 
  FindNeighbors(
    reduction = 'integrated_lsi',
    dims = 2:50
  ) %>% 
  FindClusters( 
    algorithm = 3,
    resolution = 2,
    verbose = FALSE
  )


saveRDS(integrated, paste0(add_filename, '_integrated_chemistry.rds'))

fwrite(cbind(Embeddings(integrated, reduction='umap'), integrated@meta.data), paste0(add_filename,"_chemistry_metadata.tsv"),
       sep='\t', row.names = T)

