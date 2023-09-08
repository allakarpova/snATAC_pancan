#run several rounds of clustering with different resolutions
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(optparse))


###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to object",
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
  make_option(c("--multiome"),
              type="logical",
              default=FALSE,
              help = "is it a multiome object or not?",
              metavar="logical")
  
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

select <- dplyr::select
filter <- dplyr::filter

# read in initial arguments
input.path<- opt$input.object
out_path <- opt$output
add_filename <- opt$extra

dir.create(out_path, showWarnings = F)
setwd(out_path)

panc<- readRDS(input.path)

for (resol in c(0.1, 0.3, 0.5, 0.7, seq(1, 1.8, 0.1), 2)) {
  if(opt$multiome) {
    panc <- FindClusters(panc, resolution = resol,graph.name = "wsnn", algorithm = 4, method = "igraph")
  } else {
    panc <- FindClusters(panc, resolution = resol, algorithm = 4, method = "igraph")
  }
  
  print(head(panc@meta.data))
}

cluster.tb <- panc@meta.data %>% select(dplyr::contains('res.'))
fwrite(cluster.tb, paste0('Clusters_res0.1_to_2_alg4_', add_filename, '.txt'), sep='\t', row.names = TRUE)


for (resol in c(0.1, 0.3, 0.5, 0.7, seq(1, 1.8, 0.1), 2)) {
  if(opt$multiome) {
    panc <- FindClusters(panc, resolution = resol,graph.name = "wsnn", algorithm = 1)
  } else {
    panc <- FindClusters(panc, resolution = resol, algorithm = 1)
  }
  
  print(head(panc@meta.data))
}

cluster.tb <- panc@meta.data %>% select(dplyr::contains('res.'))
fwrite(cluster.tb, paste0('Clusters_res0.1_to_2_alg1_', add_filename, '.txt'), sep='\t', row.names = TRUE)

for (resol in c(0.1, 0.3, 0.5, 0.7, seq(1, 1.8, 0.1), 2)) {
  if(opt$multiome) {
    panc <- FindClusters(panc, resolution = resol,graph.name = "wsnn", algorithm = 3)
  } else {
    panc <- FindClusters(panc, resolution = resol, algorithm = 3)
  }
  
  print(head(panc@meta.data))
}
cluster.tb <- panc@meta.data %>% select(dplyr::contains('res.'))
fwrite(cluster.tb, paste0('Clusters_res0.1_to_2_alg3_', add_filename, '.txt'), sep='\t', row.names = TRUE)





