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
              default='cell_type.harmonized',
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

int <- readRDS(input.path)

my.metadata <- fread(meta.path, data.table = F, header = T) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)
int <- AddMetaData(int, metadata = my.metadata)
int$is_Singlet <- !grepl('Doublet', as.character(unlist(int[[cell_column]])))

head(int@meta.data, n=1)

column <- 'no_doublet'
int.sub <- subset(x = int, subset = is_Singlet)
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


dimplot=DimPlot(object = int.sub, label = TRUE, group.by = cell_column)
pdf(paste0('Dimplot_',cell_column,'_', column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int.sub, label = TRUE, group.by = "Cancer", cols = 'Spectral')
pdf(paste0('Dimplot_Cancer_',column,"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()