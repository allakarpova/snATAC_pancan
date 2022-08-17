## Alla Karpova
### recluster T-cells in separate objects

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))


suppressMessages(library(future))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))


option_list = list(
  make_option(c("-i", "--input.path"),
              type="character",
              default=NULL, 
              help="path to folder with rds objects",
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
              default='cell_type_cluster',
              help = "cell_type_column",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.path
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


plan("multicore", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM

int <- readRDS(input.path)

my.metadata <- fread(meta.path, data.table = F, header = T) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)
colnames(my.metadata)[which(colnames(my.metadata)==cell_column)] <- paste(cell_column, 'atac',sep = '.')
head(my.metadata)
my.metadata <- my.metadata %>% select(ends_with('atac'))
head(my.metadata)
int <- AddMetaData(int, metadata = my.metadata)

int$touse <- int@meta.data[, paste(cell_column, 'atac',sep = '.')]
#NOW extract lineages

int$lymphocyte.type <- case_when(grepl('CD8-', int$touse) ~ 'CD4/CD8- T-cells',
                                 grepl('CD8', int$touse) ~ 'CD8 T-cells',
                                 grepl('CD4|Treg', int$touse) ~ 'CD4 T-cells',
                                 grepl('NK', int$touse) ~ 'NK cells',
                                 grepl('NKT', int$touse) ~ 'NKT cells',
                                 TRUE ~ as.character(int$touse))
table(int$lymphocyte.type)
conditions <- c( 'CD4 T-cells', 'NK cells', 'CD8 T-cells') 

conditions %>% walk (function(column) {
  int.sub <- subset(x = int, lymphocyte.type == column)
  int.sub <- int.sub %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:30, assay = 'SCT') %>%
    NormalizeData(assay = 'RNA') %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(verbose = T, resolution = 1)
  cat('saving the object...\n')
  saveRDS(int.sub, paste0("snRNA_combo_Merged_",make.names(column),"_",add_filename,".rds"))
  
  
  dimplot=DimPlot(object = int.sub, label = TRUE) + NoLegend()
  pdf(paste0('Dimplot_', make.names(column),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=7,width=8, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = F, group.by = "Piece_ID")
  pdf(paste0('Dimplot_Piece_ID_',make.names(column),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=7,width=12, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = F, group.by = "data.type")
  pdf(paste0('Dimplot_data.type_',make.names(column),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=7,width=8, useDingbats = F)
  print(dimplot)
  dev.off()
  
  dimplot=DimPlot(object = int.sub, label = TRUE, group.by =  paste(cell_column, 'atac',sep = '.'))
  pdf(paste0('Dimplot_cell_type.atac_', make.names(column),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=7,width=8, useDingbats = F)
  print(dimplot)
  dev.off()
  
  
  dimplot=DimPlot(object = int.sub, label = TRUE, group.by = "Cancer")
  pdf(paste0('Dimplot_Cancer_',make.names(column),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=7,width=8, useDingbats = F)
  print(dimplot)
  dev.off()
})


