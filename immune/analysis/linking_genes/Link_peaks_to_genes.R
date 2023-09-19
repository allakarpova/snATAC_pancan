# Alla Karpova run LinkPeaks on a multiome object immune cells

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =10)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))


###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="input object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')


out.obj <- str_replace(input_path, pattern = 'rds', replacement = 'linked.rds')
out.obj <- str_split(out.obj, '[/]')[[1]][str_count(out.obj, '/')+1]

cat('opening object \n')
obj <- readRDS(input_path)
DefaultAssay(obj) <- "ATAC_immune"
cat('done \n')

annot <- readRDS('/diskmnt/Projects/snATAC_analysis/immune/conda_env_files/Annotations.EnsDb.Hsapiens.v100.rds')
genome(annot) <- "NA"
seqlevelsStyle(annot) <- 'UCSC'
genome(annot) <- "hg38"
Annotation(obj) <- annot


min.cells.num <- 0.00005*ncol(obj)
print(min.cells.num)

# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "ATAC_immune",
  distance = 5e+05,
  n_sample = 1000,
  #min.cells = min.cells.num,
  expression.assay = "SCT"
)

toreturn <- Links(obj)
print(out.obj)
saveRDS(toreturn, paste0('out/', out.obj))
