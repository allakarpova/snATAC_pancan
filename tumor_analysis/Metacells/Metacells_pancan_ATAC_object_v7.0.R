suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
library(hdWGCNA)
library(WGCNA)


SelectFractionGenes <- function(
    seurat_obj,
    fraction=0.05,
    group.by=NULL # should be a column in the Seurat object, eg clusters
    
){
  # get assay
  assay <- DefaultAssay(seurat_obj)
  tryCatch(
    tmp <- dim(seurat_obj[[assay]]@counts),
    error = function(e){
      print("Assay must contain counts slot.")
    })
  
  # handle different selection strategies
  
  # binarize counts matrix in chunks to save memory
  expr_mat <- GetAssayData(seurat_obj, slot='counts')
  expr_mat <- expr_mat > 0
  # n_chunks <- ceiling(ncol(expr_mat) / 5000)
  # 
  # if(n_chunks == 1){
  #   chunks <- factor(rep(1), levels=1)
  # } else{
  #   chunks <- cut(1:nrow(expr_mat), n_chunks)
  # }
  # expr_mat <- do.call(rbind, lapply(levels(chunks), function(x){
  #   cur <- expr_mat[chunks == x,]
  #   cur[cur > 0] <- 1
  #   cur
  # }))
  
  group_gene_list <- list()
  if(!is.null(group.by)){
    # identify genes that are expressed
    groups <- unique(seurat_obj@meta.data[[group.by]])
    for(cur_group in groups){
      
      # subset expression matrix by this group
      cur_expr <- expr_mat[,seurat_obj@meta.data[[group.by]] == cur_group]
      print(dim(cur_expr))
      
      gene_filter <- Matrix::rowSums(cur_expr) >= round(fraction*ncol(cur_expr));
      group_gene_list[[cur_group]] <- rownames(seurat_obj)[gene_filter]
    }
    gene_list <- unique(unlist(group_gene_list))
    
  } else{
    # identify genes that are expressed in at least some fraction of cells
    gene_filter <- Matrix::rowSums(expr_mat) >= round(fraction*ncol(seurat_obj));
    gene_list <- rownames(seurat_obj)[gene_filter]
  }
  
  return(gene_list)
  
}

###options###
######################
option_list = list(

  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-f", "--fraction"),
              type="double",
              default=0.01, 
              help="fraction cutoff",
              metavar="double")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

frac <- opt$fraction
out_path <- opt$output
dir.create(out_path, showWarnings = F)
setwd(out_path)


obj <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/PanCan_merged_obj/Pancan_225Samples_1000cellsPerCluster_TumorNormal.200PCs.chromVar.20230118.rds.gz')

atac.meta <- fread('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/All_225_samples_metadata_data_freeze_v7.0_new_columns_added.tsv') %>%
  data.frame(row.names = 1)
obj <- AddMetaData(obj, atac.meta)



DefaultAssay(obj) <- 'pancan'
genes.oi <- SelectFractionGenes(obj,
                                fraction=frac, # fraction of cells that a gene needs to be expressed in order to be included
                                group.by = c("Piece_cell_type.normal"))
atac.obj <- SetupForWGCNA(
  obj,
  gene_list=genes.oi,
  gene_select = "custom", # the gene selection approach
  wgcna_name = 'metacells' # the name of the scWGCNA experiment
)

atac.obj <- MetacellsByGroups(
  reduction = 'umap',
  seurat_obj = obj,
  assay = 'pancan',
  group.by = c("Piece_cell_type.normal"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  ident.group = 'Piece_cell_type.normal' # set the Idents of the metacell seurat object
)

meta.atac.obj <- obj@misc$metacells$wgcna_metacell_obj

meta.atac.obj <- meta.atac.obj %>% 
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 500) %>%
  RunSVD(
    reduction.key = 'LSI_',
    reduction.name = 'lsi',
    irlba.work = 400) %>%
  RunUMAP(dims = 2:30,reduction = 'lsi')


saveRDS('Metacells_Pancan_225Samples_1000cellsPerCluster_TumorNormal.200PCs.chromVar.20230118.rds')
DimPlot(meta.atac.obj, label = TRUE)
ggsave('Dimplot_metacell.obj.pdf', width=10, height = 10)




