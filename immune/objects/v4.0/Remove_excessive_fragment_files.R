suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(googlesheets4))
set.seed(1234)
suppressMessages(library(optparse))


option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to where you want to update fragment paths",
              metavar="character"),
  make_option(c("-o", "--output.object"),
              type="character",
              default=NULL, 
              help="output object name",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

###################################

# read in initial arguments
input.path <- opt$input.object
output.path <- opt$output.object

panc <- readRDS(input.path)

all.fragment.obj <- Fragments(panc)
all.fragment.obj.cell.count <- map_chr(all.fragment.obj, function(x) length(x@cells))
all.fragment.obj.upd <- all.fragment.obj[all.fragment.obj.cell.count > 0]

Fragments(panc) <- NULL
Fragments(panc) <- all.fragment.obj.upd


saveRDS(panc, output.path)



