#run ROGUE

suppressMessages(library(ROGUE))
library(future)

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
library(glue)

option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object",
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
  make_option(c("-c","--cluster.file"),
              type="character",
              default=NULL,
              help = "path to file with clusters at different resolutions",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra


dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

cat('opening object \n')
rna.obj <- readRDS(input.path)
clusters.all <- fread(opt$cluster.file) %>%
  data.frame(row.names = 1) %>% 
  select(starts_with('integrated_snn_res'))

expr <- GetAssayData(rna.obj, assay = 'RNA', slot = 'counts')


ent.res <- SE_fun(expr)
head(ent.res)

fwrite(ent.res, glue('Entropy_from_ROGUE_{add_filename}.tsv'), sep='\t')

pdf(glue::glue('SEplot.pdf'), width=7, height = 6)
SEplot(ent.res)
dev.off()

registerDoParallel(cores=25)
foreach(column=clusters.all) %dopar% {
  
  #rna.obj <- AddMetaData(rna.obj, clusters.all[[column]], col.name = 'ct')
  meta <- rna.obj@meta.data %>% 
    bind_cols((clusters.all %>% select(all_of(column)))) %>%
    select(Piece_ID_RNA, ct) %>% 
    rename(Patient=Piece_ID_RNA)
  print(head(meta))
  
  rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)
  print(head(rogue.res))
  
  fwrite(rogue.res, glue('ROGUE_{column}_{add_filename}.tsv'), sep='\t', row.names = T)
  #saveRDS(rogue.res, glue('ROGUE_{column}_{add_filename}.rds'))
  
  rogue.boxplot(rogue.res)
  ggsave(glue::glue('Boxplot_{column}.pdf'), width=18, height = 5)
}
stopImplicitCluster()
  

