# find DEGs one cancer type vs the rest
# Alla Karpova
###libraries
##################

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(future))
plan("multicore", workers = 5)
options(future.globals.maxSize = 100 * 1024^3)

################################

###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/10.snRNA/1.Merge_snRNA/Merge.v.20220307/Merged_132_snRNA.SCTrandformed.3KVarFeatures.v3.20220307.rds', 
              help="path to merged object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="blah", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("--metadata"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/04_celltyped_rds/cell_type_snRNA_merged/v5.0_data_freeze/All_RNA_cell_type.harmonized.cancer.meta.data', 
              help="path to metadata")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata

dir.create(out_path, showWarnings = F)
setwd(out_path)

meta <- fread(meta.path, header = TRUE) %>% data.frame(row.names = 1)
meta <- meta %>% mutate(Cancer.subtype = case_when(grepl('HT029B1|HT035B1|1408|HT141B1|HT271B1|HT268B1', Piece_ID) ~'Basal',
                                                   grepl('HT214|HT265', Piece_ID) ~ 'Her2-enriched',
                                                   Cancer == 'BRCA' ~ 'Luminal',
                                                   TRUE ~ 'NA'))
meta$Cancer.new <- case_when(meta$Cancer.subtype=='Basal' ~ 'BRCA_Basal',
                            TRUE ~ meta$Cancer)
meta$Barcodes.cancer <- rownames(meta)
#print(head(meta))
rownames(meta) <- case_when(meta$Cancer == 'CRC' ~ str_split_fixed(rownames(meta), '_', 2)[,2],
                            TRUE ~ rownames(meta))
print(head(meta %>% filter(Cancer == 'CRC')))

cat('opening object \n')
panc <- readRDS(input.path)
cat('done \n')

panc <- AddMetaData(panc, meta)

cat('subsetting object \n')
panc <- subset(panc, cell_type.harmonized.cancer=='Tumor')

Idents(panc) <- 'Cancer.new'
DefaultAssay(panc) <- 'SCT'
unique(Idents(panc))

deg <- FindAllMarkers(
  object = panc,
  only.pos = F,
  min.pct = 0.1,
  min.diff.pct = 0,
  assay = 'SCT',
  logfc.threshold = 0
)

fwrite(deg, paste0('degs_oneCances_vs_Others.minPct0.1_SCT_', add_filename,'.txt'), sep = '\t', row.names = T)

DefaultAssay(panc) <- 'RNA'
panc <- NormalizeData(panc)
deg <- FindAllMarkers(
  object = panc,
  only.pos = F,
  min.pct = 0.1,
  min.diff.pct = 0,
  assay = 'RNA',
  logfc.threshold = 0
)

fwrite(deg, paste0('degs_oneCances_vs_Others.minPct0.1_RNA_', add_filename,'.txt'), sep = '\t', row.names = T)











