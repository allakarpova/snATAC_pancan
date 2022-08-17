### merge combo objects for ATAC pancancer

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
  make_option(c("-i", "--input.folder"),
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
              metavar="character"),
  make_option(c("--meta.rownames"),
              type="character",
              default='Barcodes_cancer',
              help = "which column of metadata contains correct rownames",
              metavar="character"),
  make_option(c("--filter_cells"),
              type="logical",
              default=TRUE,
              help = "should I filter out non-cells by cellranger?",
              metavar="logical")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


plan("multiprocess", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in immune` %>% unlist()
#samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples$Cancer_piece = paste(samples$`Disease Type`, samples$Piece_ID, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::filter( `Data Type` == '10x_SC_Multi_ATAC_SEQ')
samples <- samples %>% dplyr::select(Sample,Cancer_piece, `Data Type`, `Data folder`)

samples.id <- samples$Sample %>% as.character()
samples.type <- samples$`Data Type` %>% as.character()
samples.piece <- samples$Cancer_piece %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))
if (!file.exists(paste0(length(samples.id),"_snRNA_combo_Merged_all_",add_filename,".rds"))) {
  paths <- NULL
  for (i in 1:length(samples.id)){
    print(samples.id[i])
    p <- list.files(path = input.path, full.names = T, pattern = paste0(str_split_fixed(samples.id[i], '_',2)[2],'.*rds'), all.files = T, recursive = T)
    print(length(p))
    paths <- c(paths, p)
  }
  
  # make the list of atac objects
  registerDoParallel(cores=10)
  cat ('Reading in objects\n')
  
  atac <- foreach (i=1:length(samples.id), p = paths, .combine=c) %dopar% {
    print(samples.id[i])
    obj=readRDS(p)
    DefaultAssay(obj) <- 'RNA'
    #if (!file.exists(Fragments(obj)[[1]]@path)) stop("Urgh, this sample object can't locate fragments file")
    obj<- DietSeurat(obj, assay = 'RNA')
    obj$dataset = samples.id[i]
    obj$Data.type = samples.type[i]
    return(obj)
  }
  stopImplicitCluster()
  
  cat ('Merging\n')
  int <- merge(x = atac[[1]], y = atac[-1], add.cell.ids = samples.piece)  
  saveRDS(int,   paste0(length(samples.id),"_snRNA_combo_Merged_all_",add_filename,".rds"))
  
} else {
  int <- readRDS(paste0(length(samples.id),"_snRNA_combo_Merged_all_",add_filename,".rds"))
}

# add meta data if provided
if (!is.null(meta.path)) {
  my.metadata <- fread(meta.path, data.table = F) %>% data.frame(row.names = 'Barcodes_cancer', check.rows = F, check.names = F)
  #rownames(my.metadata) <- my.metadata[opt$meta.rownames]
  int <- AddMetaData(int, metadata = my.metadata)
}

int$is_immune <- grepl('B-cell|Plasma|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune|Kuppfer|Lympho|Neutro|Gran', as.character(unlist(int[[cell_column]]))) |
  int$Cancer == 'PBMC'

cat ('Subsetting\n')
int <- subset(x = int, subset = is_immune)

if (opt$filter_cells) {
  cat ('Subsetting cells only\n')
  int <- subset(x = int, subset = is__cell_barcode==1)
}

samples.id <- samples.id[samples.id %in% int$Sample]
int[["percent.mt"]] <- PercentageFeatureSet(int, pattern = "^MT-")
int <- CellCycleScoring(int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)

int <- int %>%
  SCTransform(
    assay = 'RNA',
    vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
    conserve.memory = T,
    return.only.var.genes = T
  ) %>%
  RunPCA(assay = 'SCT', do.print = FALSE) %>%
  RunUMAP(dims = 1:50, assay = 'SCT')

int <- NormalizeData(int, assay = 'RNA')
int <- FindNeighbors(object = int,  dims = 1:50)
int <- FindClusters(object = int, resolution = 2, verbose = FALSE)


cat('saving the object...\n')
saveRDS(int,   paste0(length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".rds"))

dimplot=DimPlot(object = int, label = TRUE) + NoLegend()
pdf(paste0('Dimplot_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "dataset")
pdf(paste0('Dimplot_dataset_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized")
pdf(paste0('Dimplot_cell_type.harmonized_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized.cancer")
pdf(paste0('Dimplot_cell_type.harmonized.cancer_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "Cancer", cols = 'Spectral')
pdf(paste0('Dimplot_Cancer_', length(samples.id),"_snRNA_combo_Merged_immune_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()
