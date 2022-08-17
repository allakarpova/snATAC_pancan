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

gs4_deauth()
fragments.sheet <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 18, trim_ws = T)

all.fragment.obj <- Fragments(panc)
all.fragment.paths <- map_chr(all.fragment.obj, function(x) x@path)

new.paths <- fragments.sheet$`Fragment file compute 1`[match(all.fragment.paths, fragments.sheet$`Fragment file katmai`)]

#Update paths
new.paths <- as.list(new.paths)
for (i in seq_along(all.fragment.obj)) {
  all.fragment.obj[[i]] <- UpdatePath(all.fragment.obj[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(panc) <- NULL 
Fragments(panc) <- all.fragment.obj # assign updated list back to the object
Fragments(panc)

saveRDS(panc, output.path)



