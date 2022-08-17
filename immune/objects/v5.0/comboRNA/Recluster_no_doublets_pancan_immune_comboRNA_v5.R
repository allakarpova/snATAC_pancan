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

int$to_remove <- grepl('oublet', as.character(unlist(int[[cell_column]])))

int <- subset(x = int, subset = to_remove, invert = TRUE)

int[["percent.mt"]] <- PercentageFeatureSet(int, pattern = "^MT-")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

int <- CellCycleScoring(int, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

int <- int %>%
  SCTransform(
    assay = 'RNA',
    vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
    conserve.memory = T,
    return.only.var.genes = T
  ) %>%
  RunPCA(assay = 'SCT', do.print = FALSE) %>%
  RunUMAP(dims = 1:30, assay = 'SCT')

int <- NormalizeData(int, assay = 'RNA')
int <- FindNeighbors(object = int,  dims = 1:30)
int <- FindClusters(object = int, resolution = 1, verbose = FALSE)


cat('saving the object...\n')
saveRDS(int,   paste0("PanImmune_merged_object_snRNA_combo_",add_filename,".rds"))

dimplot=DimPlot(object = int, label = TRUE) + NoLegend()
pdf(paste0('Dimplot_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "Sample")
pdf(paste0('Dimplot_sample_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "data.type")
pdf(paste0('Dimplot_data.type_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=14, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized.cancer")
pdf(paste0('Dimplot_cell_type.harmonized.cancer_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = cell_column)
pdf(paste0('Dimplot_',cell_column,'_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "Cancer")
pdf(paste0('Dimplot_Cancer_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

