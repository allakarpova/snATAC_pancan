## Alla Karpova
### merge combo objects for ATAC pancancer

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


option_list = list(
  make_option(c("-i", "--input.obj"),
              type="character",
              default=NULL,
              help="path to input object",
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
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.obj
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in immune` %>% unlist()
samples <- samples %>% dplyr::filter(Keep == 'TRUE') %>% 
  dplyr::filter(`Cellranger version` == 'v2.0') %>% 
  dplyr::filter( `Data Type` == '10x_SC_Multi_ATAC_SEQ')
samples.id <- samples$Sample %>% as.character()

plan("multicore", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM


int <- readRDS(input.path)
# add meta data if provided
if (!is.null(meta.path)) {
  my.metadata <- fread(meta.path, data.table = F) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F)
  int <- AddMetaData(int, metadata = my.metadata)
}



library(harmony)
#NOW extract lineages

int$is_Tcell_lin <- grepl('T-cell|NK|Treg|MAIT', as.character(unlist(int[[cell_column]])))
int$is_Bcell_lin <- grepl('B-cell|Plasma', as.character(unlist(int[[cell_column]])))
int$is_Myeloid_lin <- grepl('Macro|DC|Mast|Microgl|Mono', as.character(unlist(int[[cell_column]])))
int[["percent.ribo"]] <- PercentageFeatureSet(int, pattern = "^RP[S,L][1-9]+$")


conditions <- c('is_Tcell_lin', 'is_Bcell_lin', 'is_Myeloid_lin') 

conditions %>% walk (function(column) {
  int.sub <- subset(x = int, cells = rownames(dplyr::filter(int@meta.data, .data[[column]])))
  int.sub <- int.sub %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score", 'percent.ribo'),
      conserve.memory = T,
      return.only.var.genes = T
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:30, assay = 'SCT') %>%
    NormalizeData(assay = 'RNA') %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(verbose = T, resolution = 2)
  cat('saving the object...\n')
  saveRDS(int.sub, paste0("snRNA_combo_Merged_",column,"_",add_filename,".rds"))
  
  
  dimplot=DimPlot(object = int.sub, label = TRUE) + NoLegend()
  pdf(paste0('Dimplot_', column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = F, group.by = "Sample")
  pdf(paste0('Dimplot_sample_',column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = F, group.by = "data.type")
  pdf(paste0('Dimplot_data.type_',column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=14, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = TRUE, group.by = cell_column)
  pdf(paste0('Dimplot_',cell_column, '_', column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = TRUE, group.by = "cell_type.harmonized.cancer")
  pdf(paste0('Dimplot_cell_type.harmonized.cancer_', column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = TRUE, group.by = "Cancer")
  pdf(paste0('Dimplot_Cancer_',column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
  print(dimplot)
  dev.off()
})
