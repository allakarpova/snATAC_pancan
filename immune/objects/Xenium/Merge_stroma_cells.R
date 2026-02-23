## Alla Karpova
### merge Xenium objects from PanImmune project but only stromal cells
# v1.0 02.23.2025  

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(harmony))
library(doParallel)
options(future.globals.maxSize= 1610612736) #1.5Gb

############## FUNCTIONS #####################

option_list = list(
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
  make_option(c("-v","--panel.version"),
              type="character",
              default="5K",
              help = "5K",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments

out_path <- opt$output
add_filename <- opt$extra

ver <- opt$panel.version

select <- dplyr::select

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1-YafCxG0eDfLF0vNISZPsx1lqsM74BI9kWED5MFrOLg/edit?gid=0#gid=0", 
                      sheet = "Xenium T-cells", trim_ws = T)

samples <- samples %>% dplyr::filter(`Can_use` == 'Yes')


if(ver=='5K') {
  samples <- samples %>% dplyr::filter(`5k`)
} else {
  samples <- samples %>% dplyr::filter(!`5k`)
}
print(samples)
samples.id <- samples$Sample %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))

options(Seurat.object.assay.version = "v3")
if (!file.exists(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))) {
  cat('creating object \n')
  paths <- samples$`Object path`
  meta.paths <- samples$`Cell type path`
  
  print(paths)
  print(meta.paths)
  missing_paths <- paths[!file.exists(paths)]
  missing_paths2 <- meta.paths[!file.exists(meta.paths)]
  
  if (length(c(missing_paths, missing_paths2)) > 0) {
    message("These paths do not exist:")
    print(c(missing_paths,missing_paths2))
  } else {
    message("All paths exist ✅")
  }
  
  # make the list of objects
  registerDoParallel(cores=10)
  cat ('Reading in objects\n')
  
  obj <- foreach (s=samples.id, 
                  p = paths, 
                  pid = samples$Case, 
                  c = samples$Cancer, 
                  subtype = samples$Subtype, 
                  meta = samples$`Cell type path`,
                  v=samples$Panel, .combine=c) %dopar% {
                    print(s)
                    obj=readRDS(p) 
                    obj@images[[1]] <- NULL
                    obj@images[[1]] <- NULL
                    print(paste('opened', s))
                    meta.df <- fread(meta, header = T) %>% column_to_rownames('cell_id')
                    
                    DefaultAssay(obj) <- 'Xenium'
                    obj <- DietSeurat(obj, assays = 'Xenium')
                    obj <- AddMetaData(obj, meta.df)
                    obj <- subset(obj, group %in% c('Endothelial cells', 'Fibroblasts', 'vSMCs', 'Pericytes', 'Portal fibroblasts'))
                    obj$Sample = s
                    obj$Case = pid
                    obj$Cancer = c
                    obj$Cancer_subtype = subtype
                    obj$Panel_version = v
                    return(obj)
                  }
  stopImplicitCluster()
  
  combined <- merge(x = obj[[1]], y = obj[-1], add.cell.ids = samples.id)  
  saveRDS(combined,  paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))
  
} else {
  combined <- readRDS(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))
}

cat('normalizing Xenium\n')
DefaultAssay(combined) <- 'Xenium'
options(Seurat.object.assay.version = "v3")

combined <- combined %>%
  SCTransform(
    ncells = 20000,
    assay = 'Xenium',
    return.only.var.genes = TRUE, 
    verbose = T) %>%
  RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
  RunUMAP(dims = 1:30, verbose = T) %>%
  RunHarmony('Sample', reduction = 'pca', assay.use = 'SCT') %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5, verbose = FALSE) %>%
  RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:30)


cat('saving the object...\n')
saveRDS(combined,  paste0(length(samples.id),"_Merged_normalized_",add_filename,".rds"))


DimPlot(combined,  group.by = "seurat_clusters", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_clusters_", add_filename, ".pdf"),height=10,width=11)
DimPlot(combined,  group.by = "Sample", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_Sample_", add_filename, ".pdf"),height=10,width=12)
DimPlot(combined,  group.by = "Cancer", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_Cancer_", add_filename, ".pdf"),height=10,width=11)
DimPlot(combined,  group.by = "group", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_cell.type_", add_filename, ".pdf"),height=10,width=11)

fwrite(combined@meta.data, 
       paste0(length(samples.id),"_Merged_normalized_",add_filename,".metadata.tsv"), sep = '\t', row.names = T)



