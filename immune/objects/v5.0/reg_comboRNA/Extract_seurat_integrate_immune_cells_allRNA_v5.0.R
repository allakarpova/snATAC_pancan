## Alla Karpova
### merge regular RNA data for immune cells and integrate with seurat

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



runAllNormalization <- function(obj, dims=30) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims, assay = 'SCT')
  
  obj <- NormalizeData(obj, assay = 'RNA')
  obj <- FindNeighbors(object = obj,  dims = 1:dims)
  obj <- FindClusters(object = obj,resolution = 2, verbose = FALSE)
  return(obj)
  
}

filter <- dplyr::filter
select <- dplyr::select


option_list = list(
  # make_option(c("-i", "--input.folder"),
  #             type="character",
  #             default=NULL,
  #             help="path to folder with cancer level merged rds objects",
  #             metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-t","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type.harmonized',
              help = "cell_type_column",
              metavar="character"),
  make_option(c("--comboRNA"),
              type="character",
              default='/diskmnt/Projects/snATAC_analysis/immune/obj/v5.0_comboRNA/no_doublets/PanImmune_merged_object_snRNA_combo_v5.0_df_no_doublets.rds',
              help = "path to combo RNA immune object",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
#input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column
comboRNA.path <- opt$comboRNA

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


plan("multicore", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in immune` %>% unlist()

samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
#samples <- samples %>% dplyr::filter( `Data Type` != '10x_SC_Multi_ATAC_SEQ')
#samples <- samples %>% dplyr::select(Sample,Piece_ID, `Data Type`, `Data folder`)

cancers <- samples %>% dplyr::filter( `Data Type` != '10x_SC_Multi_ATAC_SEQ') %>% 
  pull(`Disease Type`) %>% unique()

print(cancers)
merged.obj <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", trim_ws = T, sheet = "matched snRNA merged objects")
merged.obj.path <- merged.obj$`Merged object RDS`
names(merged.obj.path) <- merged.obj$`Disease Type`

merged.obj.path <- merged.obj.path[names(merged.obj.path) %in% cancers]

print(merged.obj.path)

if(!file.exists(paste0('PanImmune_merged_regularRNA_normalized_', add_filename, '.rds'))) {
  
  cat('creating regular RNA object \n')

  
  my.metadata <- fread(meta.path, data.table = F) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F)
  my.metadata$is_immune <- grepl('Plasma|B-cells|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune|Gran|Mono|MAIT', as.character(unlist(my.metadata[[cell_column]])))
  cancers.with.immune <- my.metadata %>% filter(is_immune & data.type == 'snATAC') %>% pull(Cancer)
  cancers.with.immune <- cancers.with.immune[-which(cancers.with.immune=='OV')] #remove ovarian cancer because the only snRNA sample with immune cells is weird VF031V1-Tm1Y1
  merged.obj.path <- merged.obj.path[names(merged.obj.path) %in% cancers.with.immune]
  
  
  cat('opening cancer object...\n')
  rna <- map(merged.obj.path, function(p) {
    print(names(merged.obj.path[merged.obj.path==p]))
    obj=readRDS(p)
    DefaultAssay(obj) <- 'RNA'
    obj<- DietSeurat(obj, assay = 'RNA', counts = TRUE, data = TRUE)
    
    if(names(merged.obj.path[merged.obj.path==p])=='CRC') {
      if (!grepl('CRC_', colnames(obj)[1])) {
        obj <- RenameCells(object = obj, new.names = paste('CRC', colnames(obj), sep = '_')) # I hate that it is different
      }
      #print(head(obj@meta.data))
    }
    
    obj <- AddMetaData(obj, my.metadata)
    print(head(obj@meta.data))
    ct <- obj@meta.data %>% pull(cell_column) %>% unique %>% sort
    print(ct)
    stopifnot(!is.na(ct))
    
    cat('subsetting\n')
    
    obj.my <- subset(x = obj, subset = (is_immune & data.type == 'snATAC'))
    print(unique(obj.my$Case_ID))
    if(c('HT029B1') %in% obj.my$Case_ID) {
      cat('removing HT029B1 \n')
      obj.my <- subset(x = obj.my, subset = Case_ID %in% c('HT029B1'), invert = TRUE)
    }
    if(c('MMY34600') %in% obj.my$Case_ID) {
      cat('removing MMY34600 \n')
      obj.my <- subset(x = obj.my, subset = Case_ID %in% c( 'MMY34600'), invert = TRUE)
    }
    if(c( 'C3N-01816') %in% obj.my$Case_ID) {
      cat('removing C3N-01816\n')
      obj.my <- subset(x = obj.my, subset = Case_ID %in% c('C3N-01816'), invert = TRUE)
    }
    if (c('VF031V1-Tm1Y1') %in% obj.my$Piece_ID) {
      cat('removing VF031V1-Tm1Y1\n')
      obj.my <- subset(x = obj.my, subset = Piece_ID %in% c('VF031V1-Tm1Y1'), invert = TRUE)
    }
    
    print(dim(obj.my))
    return(obj.my)
  })
  cat('done\n')
  
  cat ('Merging regular RNA\n')
  reg.rna <- merge(x = rna[[1]], y = rna[-1])
  cat ('done\n')
  
  #remove individual objects
  rm(rna)
  gc()
  
  saveRDS(reg.rna, paste0('PanImmune_merged_regularRNA_', add_filename, '.rds'))
  
  cat ('Normalizing merged regular RNA object\n')
  reg.rna <- runAllNormalization (reg.rna, dims = 30)
  ################
  
  saveRDS(reg.rna, paste0('PanImmune_merged_regularRNA_normalized_', add_filename, '.rds'))
} else {
  reg.rna <- readRDS(paste0('PanImmune_merged_regularRNA_normalized_', add_filename, '.rds'))
}


DefaultAssay(reg.rna) <- 'RNA'
combo.rna <- readRDS(comboRNA.path) 
DefaultAssay(combo.rna) <- 'RNA'
combo.rna <- combo.rna %>%  DietSeurat(assay = 'RNA', counts = TRUE, data = TRUE)
reg.rna <- reg.rna %>% DietSeurat(assay ='RNA', counts = TRUE, data = TRUE)


all.rna <- merge(x=reg.rna, y=combo.rna)

remove(combo.rna, reg.rna)
gc()

all.rna$Data.source <- ifelse(all.rna$Cancer == 'PBMC', '10x', 'DingLab')
all.rna$Batches <- case_when(all.rna$Cancer %in% c('PBMC') ~ paste(all.rna$Cancer, all.rna$Chemistry, sep = '__'),
                             all.rna$Cancer %in% c('MM') ~ all.rna$Cancer,
                             TRUE ~ all.rna$Chemistry)

cat ('Integrate regular RNA and combo RNA by batches \n')
all.rna.list <- SplitObject(all.rna, split.by = 'Batches')

all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- CellCycleScoring(x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
  x <- x %>% SCTransform(
    assay = 'RNA',
    vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
    conserve.memory = T,
    return.only.var.genes = T
  )
  return(x)
})

features <- SelectIntegrationFeatures(object.list = all.rna.list, nfeatures = 3000)
all.rna.list <- PrepSCTIntegration(object.list = all.rna.list, anchor.features = features)
all.rna.list <- lapply(X = all.rna.list, FUN = RunPCA, features = features)

rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
int <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:30)
int <- RunPCA(int, verbose = FALSE)
int <- RunUMAP(int, reduction = "pca", dims = 1:30)
int <- FindNeighbors(int, reduction = "pca", dims = 1:30)
int <- FindClusters(int, resolution = 2)

tb <- table(int$seurat_clusters, int$cell_type.harmonized.cancer)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
int$cell_type.immune <- cluster.match.celltype[as.character(int$seurat_clusters)]

cat('saving the object...\n')
saveRDS(int,   paste0("PanImmune_integrated_object_regularRNA_multiome_",add_filename,".rds"))

dimplot=DimPlot(object = int, label = TRUE) + NoLegend()
pdf(paste0("Dimplot_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "Piece_ID")
pdf(paste0("Dimplot_Piece_ID_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "data.type")
pdf(paste0("Dimplot_data.type_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=14, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized.cancer")
pdf(paste0("Dimplot_cell_type.harmonized.cancer_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.immune")
pdf(paste0("Dimplot_cell_type.immune_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "Cancer")
pdf(paste0("Dimplot_Cancer_integrated_object_regularRNA_multiome_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

