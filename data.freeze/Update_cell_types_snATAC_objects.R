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
              metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
out_path <- opt$output
do.all.over <- opt$do.all.over

dir.create(out_path, showWarnings = F)
setwd(out_path)

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data folder`)

for (i in 1:nrow(samples)) {
  input.obj.path <- list.files (path = samples$`Data folder`[i], '*rds', full.names = T)
  s <- samples$Sample[i]
  print(s)
  if (!do.all.over) {
    if (file.exists(paste0(s,'_upd_cellTyped.meta.data'))) { # dont do updating if its already updated
      next()
    }
  }
  possible.files <- list.files(path = '/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v.3.0', pattern = s, full.names = T, recursive = T)
  possible.files <- possible.files[!grepl('PanCancer',possible.files) & !grepl('archived', possible.files)]
  input.neta.path <- possible.files[grepl('meta.data',possible.files)]
  if (s == "PDAC_TWCE-HT-018-P-Slice2_ATAC-lib1") { #handle special case
    input.neta.path <- input.neta.path[1]
  }
  print(input.neta.path)
  if (length(input.neta.path)==0) {
    next()
  }
  input.metadata <- fread(input.neta.path, data.table = F) %>% data.frame(row.names = 1, check.names = F)
  if (samples$`Data Type`[i]=='snATAC') {
    
    input.obj <- readRDS(input.obj.path)
    tb <- table(input.metadata$seurat_clusters, input.metadata$predicted.id)
    cluster.match.celltype <- apply(tb, 1, function(x) {
      colnames(tb)[which.max (x)]
    })
    input.metadata$cell_type <- cluster.match.celltype[as.character(input.metadata$seurat_clusters)]
    if (s =='BRCA_1408-06') { # special case for this sample, i know its tumor
      input.metadata$cell_type[input.metadata$seurat_clusters==2] <- 'Tumor'
      input.metadata$cell_type[input.metadata$seurat_clusters==7] <- 'Tumor'
    }
    if (s =='CESC_TWDD-CE334E1-N1Ba1_1') { # special case for this sample, i know its tumor
      input.metadata$cell_type[input.metadata$seurat_clusters==1] <- 'Fibroblast'
      input.metadata$cell_type[input.metadata$seurat_clusters==9] <- 'Fibroblast'
    }
    if (s =='UCEC_TWHG-CPT4427DU-XBa2_1') { # special case for this sample, i know its tumor
      input.metadata$cell_type[input.metadata$seurat_clusters %in% c(2,3)] <- 'Tumor'
    }
    if (s =='HNSCC_TWAH-WU-HTA44-OP16_CD45_sort') { # special case for this sample, i know its tumor
      input.metadata$cell_type[input.metadata$seurat_clusters %in% c(5,11)] <- 'Tumor'
    }
    if (s =='GBM_C3N-02784_CPT0206000015_2020-06-10') { # special case for this sample, i know its tumor
      input.metadata$cell_type[input.metadata$seurat_clusters %in% c(7)] <- 'Tumor'
    }
    input.obj <- AddMetaData(input.obj, input.metadata)
    DimPlot(input.obj, label = T) + DimPlot(input.obj, group.by = 'predicted.id', label = T, cols = col_vector) + DimPlot(input.obj, group.by = 'cell_type', label = T, cols = 'Paired')
    ggsave(paste0(s,'_upd_celltypes.pdf'), width = 18, height = 6, useDingbats = F)
    fwrite(input.metadata, paste0(s,'_upd_cellTyped.meta.data'), sep = '\t', row.names = T)
  } else { # basically dont update combo data metadata but sill save it
    fwrite(input.metadata, paste0(s,'_upd_cellTyped.meta.data'), sep = '\t', row.names = T)
  }
}
