# Alla Karpova create pancan combo object with RNA and ATAC

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))



######################
### FUNCTIONS #####

normalize_multiome <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, 
                     vars.to.regress = c("nCount_RNA","percent.mito"),
                     return.only.var.genes = T, verbose = FALSE, conserve.memory = TRUE) %>% 
    RunPCA(npcs = 50, verbose = FALSE) %>% 
    RunUMAP(dims = 1:50, reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  
  DefaultAssay(obj) <- "pancan"
  
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures( min.cutoff = 20) %>%
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400) %>%
    RunUMAP(dims = 2:50,reduction = 'lsi', reduction.name = "atac.umap", reduction.key = "atacUMAP_")
  
  obj <- FindMultiModalNeighbors(obj, 
                                 reduction.list = list("pca", "lsi"), 
                                 dims.list = list(1:50, 2:50))
  obj <- RunUMAP(obj, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = T)
  
  return(obj)
}


###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default="/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/Multiome/Merged_objects_PanCan_peaks/Cancer_level/out/", 
              help="input folder path",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--cell.group"),
              type="character",
              default='all', 
              help="all or tumor")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output
cell.group = opt$cell.group

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')
dir.create('plots')

annotation <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/Colors_panatac_v1.0.rds')

cat('opening objects \n')
all.objects <- list.files(input_path, pattern = '*.rds')
print(all.objects)
need.objects <- all.objects[grepl(pattern = cell.group, all.objects)]
print(need.objects)
need.objects <- paste0(input_path, '/', need.objects)

obj.list <- need.objects %>% lapply(function(x) {
  obj <- readRDS(x)
  DefaultAssay(obj) <- 'RNA'
  obj <- DietSeurat(obj, assays = c('RNA', 'pancan'))
  return(obj)
  })
cat('done \n')
  
cat('merging \n')
obj.merged <- merge(obj.list[[1]], obj.list[-1])
DefaultAssay(obj.merged) <- 'pancan'
cat('done \n')
  
cat('normalizing all \n')
obj.merged <- normalize_multiome(obj.merged)
cat('done \n')

Annotation(obj.merged) <-annotation 

saveRDS(obj.merged, paste0("out/Pancan_",cell.group ,"_cells_multiome_obj.", format(Sys.Date(), format="%Y%m%d"),".rds"))
 
p1 <- DimPlot(obj.merged, reduction = "rna.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(obj.merged, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(obj.merged, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(paste0("plots/Pancan_",cell.group ,"_cell_types_seurat_clusters.pdf"),height=8,width=24)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

DimPlot(obj.merged, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.rna", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
ggsave(paste0("plots/Pancan_",cell.group ,"_cell_types_cell_type.harmonized.cancer.rna.pdf"), height=8,width=11)

DimPlot(obj.merged, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.atac", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
ggsave(paste0("plots/Pancan_",cell.group ,"_cell_types_cell_type.harmonized.cancer.atac.pdf"), height=8,width=11)

DimPlot(obj.merged, reduction = "wnn.umap", group.by = "Piece_ID", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("plots/Pancan_",cell.group ,"_cell_types_Piece_ID.pdf"), height=8,width=11)

DimPlot(obj.merged, reduction = "wnn.umap", group.by = "Cancer", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$Cancer)
ggsave(paste0("plots/Pancan_",cell.group ,"_cell_types_Cancer.pdf"), height=8,width=11)
