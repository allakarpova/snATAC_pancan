# Alla Karpova - partially based on Nadya's code. But everything is rewritten.
## merge ATAC samples which can be either regular ATAC seq or combo ATAC seq
#v4.1 Alla implemented checking that all fragments files exist before starting merging
#v5.5 removed filtering by reproducibility, only iterative removal and removal by chrY and NNNN regions remained
# added removal of doublets
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


runAllNormalization <- function(obj, dims=30) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims, assay = 'SCT')
  
  obj <- NormalizeData(obj, assay = 'RNA')
  obj <- FindNeighbors(object = obj,  dims = 1:dims)
  obj <- FindClusters(object = obj,resolution = 1, verbose = FALSE)
  return(obj)
  
}

filter <- dplyr::filter
select <- dplyr::select



###################################
option_list = list(
  # make_option(c("-i", "--input.folder"),
  #             type="character",
  #             default=NULL, 
  #             help="path to folder with rds objects",
  #             metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="foo", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-c", "--cancer.type"),
              type="character",
              default='SKCM',
              help="cancer type to merge",
              metavar="character"),
  make_option(c("--scrublet.table"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/snATAC/scrublet/All_RNA_218_samples_metadata_immune_data_freeze_v8.0.tsv',
              help="path to scrublet metadata",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
#input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay
cancer.type <- opt$cancer.type

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 100 * 1024^3) # for 100 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = '04_matched snRNA individual rds', trim_ws = T)
samples$Keep <- samples$`Include in snRNA analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`snRNAseq cellranger ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Disease Type` == cancer.type)
samples <- samples %>% dplyr::select(Sample,Piece_ID_RNA, `Data Type RNA`, `Seurat Object folder`)

samples.id <- samples$Sample %>% as.character()
piece.id <- samples$Piece_ID_RNA %>% as.character()
samples.type <- samples$`Data Type RNA` %>% as.character()
samples.folders <- samples$`Seurat Object folder`; names(samples.folders) = samples.id

cat (paste("Samples found:" ,length(samples.id), '\n'))

doublet.table <- fread(opt$scrublet.table, header = TRUE) %>%
  data.frame(row.names = 1, stringsAsFactors = FALSE)
doublet.table$Doublet_final <- as.logical(toupper(doublet.table$Doublet_final))
str(doublet.table)
table(doublet.table$Doublet_final)

#########################################################################
paths <- map2(samples.id, samples.folders, function(s, d) {
  print(s)
  p <- list.files(path = d, full.names = T, pattern = paste0(str_split_fixed(s, '_',2)[2],'.*rds'), all.files = T, recursive = T)
  print(length(p))
  return(p)
})

#stop if not all samples have RDS object
stopifnot(length(samples.id)==length(paths))

# make the list of atac objects
registerDoParallel(cores=10)
cat ('Reading in objects\n')
#atac=vector(mode = "list", length = length(samples.id))
rna <- foreach (i=1:length(samples.id), p = paths, .combine=c) %dopar% {
  print(samples.id[i])
  obj=readRDS(p)
  DefaultAssay(obj) <- 'RNA'
  obj<- DietSeurat(obj, assay = 'RNA')
  obj$dataset = samples.id[i]
  obj$Data.type.rna = samples.type[i]
  return(obj)
}
stopImplicitCluster()

cat ('Merging\n')
combined <- merge(x = rna[[1]], y = rna[-1], add.cell.ids = paste0(cancer.type, '_', piece.id))

cat ('Removing doublets\n')
combined <- AddMetaData(combined, doublet.table)
print(head(combined@meta.data, n=3))
print(table(combined$Doublet_final))

cat('before doublet removal\n')
print(dim(combined))
combined <- subset(combined, Doublet_final, invert = TRUE)
cat('after doublet removal\n')
print(dim(combined))
#remove individual objects
rm(rna)
gc()

combined <- combined %>% runAllNormalization (dims=50)

cat('saving the object...\n')
saveRDS(combined, paste0(cancer.type,'_',length(samples.id),"_snRNA_Merged_normalized_",add_filename,".rds"))

