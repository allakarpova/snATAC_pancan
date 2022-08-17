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
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
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

cat('opening object \n')
int <- readRDS(input.path)
DefaultAssay(int) <- 'ATAC_immune'

# convert to CellDataSet format and make the cicero object
int.cds <- SeuratWrappers::as.cell_data_set(x = int)

int.cicero <- make_cicero_cds(int.cds, reduced_coordinates = reducedDims(int.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[!grepl('_', names(genome))]

# use chromosome 1 to save some time
# omit this step to run on the whole genome
#genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(int.cicero, genomic_coords = genome.df, sample_num = 100)
saveRDS(conns, paste0('cicero.res.conns.rds'))

#Find cis-co-accessible networks (CCANs)
ccans <- generate_ccans(conns, coaccess_cutoff_override = 0.25)
saveRDS(ccans, paste0('cicero.res.ccans.rds'))

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(int) <- links

saveRDS(int, paste0(add_filename, '.cicero.rds'))


