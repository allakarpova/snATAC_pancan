
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))


option_list = list(
  make_option(c("-i", "--input.folder"),
              type="character",
              default="/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v3.1", 
              help="output folder path",
              metavar="character"),
  make_option(c("-v", "--version"),
              type="character",
              default="v3.0", 
              help="version of the metadata",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input.folder <- opt$input.folder
ver <- opt$version

meta.files <- list.files(path = input.folder, pattern = '*tsv', full.names = T)
meta.files <- meta.files[grepl( ver, meta.files)]
meta.files <- meta.files[!grepl('All_', meta.files)]
meta.files <- meta.files[!grepl('archive', meta.files)]

# this script will take metadata tables from regular ATAC and combo objects and create a single metadats file
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples$`Disease Type`[samples$`Disease Type`=='PKD'] <- 'ccRCC'

samples.id <- samples$Sample %>% as.character()
length(samples.id)

total.metadata <- lapply(meta.files, function (x) {
  fread(x)})

total.colnames <- lapply(meta.files, function (x) {
  colnames(fread(x))})
common.columns <- Reduce('intersect', total.colnames)
print(common.columns)

total.metadata <- lapply(total.metadata, FUN = function(x) {
  x[,common.columns, with = F]
})

total.metadata <- rbindlist(total.metadata)
total.metadata[Cancer=='PKD', 'Cancer'] <- 'ccRCC'

fwrite(total.metadata, paste0('All_',length(samples.id),'_samples_metadata_data_freeze_', ver,'.tsv'), sep = '\t')
