suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)
suppressMessages(library(future))
plan("multicore", workers =10)
options(future.globals.maxSize = 50 * 1024^3)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(future))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(AUCell))
library(msigdbr)
library(cowplot)


Run_AUCell <- function(obj, geneset.list, name) {
  
  DefaultAssay(obj) <- 'SCT'
  exprMatrix <- GetAssayData(obj, slot = 'data')
  
  cells_rankings <- AUCell_buildRankings(exprMatrix)
  cells_AUC <- AUCell_calcAUC(geneset.list, cells_rankings,  nCores=20)
  saveRDS(cells_AUC, file=glue::glue("{name}_cells_AUC.rds"))
  
  pdf(glue::glue('{name}_AUCell_histograms.pdf'))
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
  print(cells_assignment)
  dev.off()
  
  saveRDS(cells_assignment, file=glue::glue("{name}_cells_assignment.rds"))
  
  obj[["AUCell"]] <- CreateAssayObject(getAUC(cells_AUC))
  
  pdf(glue::glue('{name}_AUCell_featureplots.pdf'), onefile = TRUE)
  #print(DimPlot(obj, group.by = 'cell_type_merged', cols = 'Dark2', label = TRUE))
  names(geneSets) %>% map(function(f) {
    tryCatch({
      print(suppressWarnings(FeaturePlot(obj, features = f)))
    }, error = function(e) {
      message(e)
      NULL
    })
  })
  dev.off()
  dev.off()
  
}




###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object to score",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="foo",
              help="add unique string identifier for your data",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
#meta.path <- opt$metadata
add_filename <- opt$extra


dir.create(out_path, showWarnings = F, recursive = TRUE)
setwd(out_path)

obj <- readRDS(input.path)
# meta <- fread(meta.path, header = TRUE) %>%
#   data.frame(row.names = 1)
#obj <- AddMetaData(obj, meta)
DefaultAssay(obj) <- 'RNA'

cat('doing one way scoring')
message('do CP genes')
c2_gene_sets = msigdbr(species = "human", category="C2")

msig_df <- c2_gene_sets %>% dplyr::filter(grepl('METABOL|TCA_|GLYCO|OXIDA|FATTY|STRESS|FOLATE', gs_name)) %>% 
  dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()
geneSets.msig <- split(msig_df$gene_symbol, msig_df$gs_name)



Run_AUCell(obj, geneSets.msig, name = add_filename)

