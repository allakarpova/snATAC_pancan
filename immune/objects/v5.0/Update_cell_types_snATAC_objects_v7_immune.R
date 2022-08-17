# after labeltransfer give each cluster cell type that is major
library(Seurat)
library(Signac)
library(data.table)
library(dplyr)
library(ggplot2)
library(optparse) 
library(googlesheets4)
library(stringr)
library(RColorBrewer)

option_list = list(
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-d", "--do.all.over"),
              type="logical",
              default=T, 
              help="Do you want to redo samples that already have their metadata updated?",
              metavar="character"),
  make_option(c("-i", "--meta.input"),
              type="character",
              default="/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v4.0", 
              help="location of metadata that needs to be updated",
              metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
out_path <- opt$output
do.all.over <- opt$do.all.over
input.meta <- opt$meta.input

dir.create(out_path, showWarnings = F)
setwd(out_path)

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in immune` %>% unlist()
#samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data folder`)

for (i in 1:nrow(samples)) {
  input.obj.path <- list.files (path = samples$`Data folder`[i], '*rds', full.names = T)
  s <- samples$Sample[i]
  print(s)
  if (!do.all.over) {
    if (file.exists(paste0(s,'_upd_v7_cellTyped.meta.data'))) { # dont do updating if its already updated
      next()
    }
  }
  #use updated metadata first
  possible.files <- list.files(path = input.meta, pattern = s, full.names = T, recursive = T)
  possible.files <- possible.files[!grepl('PanCancer',possible.files) & !grepl('archived', possible.files)]
  input.neta.path <- possible.files
  
  
  print(input.neta.path)
  if (length(input.neta.path)==0) {
    next()
  }
  #is.updated.version <- grepl('v3', input.neta.path)
  input.metadata <- fread(input.neta.path, data.table = F) %>% 
    data.frame(row.names = 1, check.names = F)
  input.metadata.celltype <- input.metadata %>% dplyr::select(cell_type) # keep cell type only
  input.obj <- readRDS(input.obj.path)
  input.obj <- AddMetaData(input.obj, input.metadata.celltype) # add existing cell type annotation
  my.reduction <- ifelse(samples$`Data Type`[i]=='snATAC' | samples$`Data Type`[i]=='scATAC', 'atac.umap', 'wnn.umap')
  DimPlot(input.obj, label = T,reduction = my.reduction) + DimPlot(input.obj,reduction = my.reduction, group.by = 'cell_type', label = T, cols = 'Paired') # plot
  ggsave(paste0(s,'_before_update_to_v7.pdf'), width = 12, height = 6, useDingbats = F)
  
  obj.meta <- input.obj@meta.data
  tb <- table(obj.meta$seurat_clusters, obj.meta$cell_type)
  cluster.match.celltype <- apply(tb, 1, function(x) {
    colnames(tb)[which.max (x)]
  })
  obj.meta$cell_type <- cluster.match.celltype[as.character(obj.meta$seurat_clusters)]
  input.obj <- AddMetaData(input.obj, obj.meta) # add updated cell type metadata
  DimPlot(input.obj, label = T,reduction = my.reduction) + DimPlot(input.obj,reduction = my.reduction, group.by = 'cell_type', label = T, cols = 'Paired')
  ggsave(paste0(s,'_upd_celltypes.pdf'), width = 12, height = 6, useDingbats = F)
  
  fwrite(input.obj@meta.data, paste0(s,'_upd_v7_cellTyped.meta.data'), sep = '\t', row.names = T)
  
}
