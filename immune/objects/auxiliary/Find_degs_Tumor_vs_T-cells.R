## Alla Karpova
### compute DEGs between cancer and T-cells 

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

filter <- dplyr::filter
select <- dplyr::select


option_list = list(
  make_option(c("-i", "--input.folder"),
              type="character",
              default=NULL,
              help="path to folder with cancer level merged rds objects",
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
input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


plan("multicore", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM


cancers <- c('BRCA', 'ccRCC', 'GBM', 'CRC', 'HNSCC', 'MM', 'CESC', 'OV', 'UCEC', "PDAC", "SKCM")
paths <- map(cancers, function(c) {
  p <- list.files(path = input.path, 
                  full.names = T, 
                  pattern = paste0(c,'.*rds'), all.files = F, recursive = T)
  print(length(p))
  return(p)
})
paths <- unlist(paths)
paths <- paths[!grepl('_snRNA_', paths)]
names(paths) <- cancers
print(paths)

my.metadata <- fread(meta.path, data.table = F) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)

cat('opening cancer object...\n')
total.degs <- map(cancers, function(c) {
  p <- paths[c]
  print(p)
  obj=readRDS(p)
  DefaultAssay(obj) <- 'SCT'
  
  obj <- AddMetaData(obj, my.metadata)

  ct <- obj@meta.data %>% pull(cell_column) %>% unique %>% sort
  print(ct)
  stopifnot(!is.na(ct))
  
  
  Idents(obj) <- 'cell_type.harmonized.cancer'
  degs <- FindMarkers(obj, ident.1 = 'Tumor', ident.2 = 'T-cells', max.cells.per.ident = 5000, only.pos =TRUE) %>%
    mutate(gene = rownames(.),
           Cancer = c)
  
  fwrite(degs, glue::glue('Degs_tumor_vs_Tcells_{c}.tsv'), sep='\t', row.names = F)
  return(degs)
})  %>% rbindlist()


total.degs %>% fwrite('Degs_tumor_vs_Tcells_all.cancers.tsv', sep='\t', row.names = F)







