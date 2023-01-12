##################
library(future)
#plan("multicore", workers = 20)
#options(future.globals.maxSize = 100 * 1024 ^ 3)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

library(cicero)
library(monocle3)
library(SeuratWrappers)
library(BSgenome.Hsapiens.UCSC.hg38)

################################


################################

###options###
######################
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
  make_option(c("-m", "--metadata"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/All_159_samples_metadata_data_freeze_v5.0.tsv', 
              help="path to metadata",
              metavar="character"),
  make_option(c("-c", "--cancer"),
              type="character",
              default=NULL, 
              help="cancer type",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('conns', showWarnings = F)
dir.create('ccans', showWarnings = F)

cat('opening object \n')
int <- readRDS(input.path)

DefaultAssay(int) <- 'pancan'

# convert to CellDataSet format and make the cicero object
int.cds <- SeuratWrappers::as.cell_data_set(x = int)

int.cicero <- make_cicero_cds(int.cds, reduced_coordinates = reducedDims(int.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[!grepl('_', names(genome))]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(int.cicero, genomic_coords = genome.df, sample_num = 100)
fwrite(conns, paste0('conns/', add_filename, '_snATAC_CICERO_conns.tsv'), sep = '\t', row.names = F)

conns <- conns %>% filter(coaccess>0.25)
fwrite(conns, paste0('conns/', add_filename, '_snATAC_CICERO_conns.0.25_cutoff.tsv'), sep = '\t', row.names = F)

#Find cis-co-accessible networks (CCANs)
ccans <- generate_ccans(conns, coaccess_cutoff_override = 0.25)
fwrite(ccans, paste0('ccans/',add_filename, '_snATAC_CICERO.CCAN.tsv'), sep = '\t', row.names = F)



