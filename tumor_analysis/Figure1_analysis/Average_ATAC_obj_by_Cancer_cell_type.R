#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))


###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/Merge.vers.20220207/PanCan_object/CellsSampled.v.20220220/Pancan_700cellsPerCluster_200PCs.chromVar.20220221.rds.gz', 
              help="path to merged object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--metadata"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v5.0_data_freeze/All_159_samples_metadata_data_freeze_v5.0.tsv', 
              help="path to metadata")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
meta.path <- opt$metadata

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')

cat('opening object \n')
RNA=readRDS(input.path)
cat('done \n')

rna.meta <- fread(meta.path) %>% data.frame(row.names = 1)

rna.meta <- rna.meta %>% mutate(Cancer.subtype = case_when(grepl('HT029B1|HT035B1|1408|HT141B1|HT271B1|HT268B1', Piece_ID) ~'Basal',
                                                           grepl('HT214|HT265', Piece_ID) ~ 'Her2-enriched',
                                                           Cancer == 'BRCA' ~ 'Luminal',
                                                           TRUE ~ 'NA'),
                                Cancer = case_when(Cancer.subtype=='Basal' ~ 'BRCA_Basal',
                                                   TRUE ~ Cancer))
RNA <- AddMetaData(RNA, rna.meta)
RNA$Cancer_cell_type <- paste(RNA$Cancer, RNA$cell_type.harmonized.cancer, sep = '_')

DefaultAssay(RNA)='pancan_s'
Idents(RNA) <- 'Cancer_cell_type'


rna.aver <- AverageExpression(RNA, slot = 'data',assays = c('pancan_s'), return.seurat = T)

saveRDS(rna.aver, paste0('out/Averaged_ATAC_by_Cancer_cell_type_Pancan_700cellsPerCluster_200PCs.chromVar.20220221.', format(Sys.Date(),  format="%Y%m%d"), '.rds'))

DefaultAssay(RNA) <- 'chromvar'

RNA@assays$chromvar@counts <- RNA@assays$chromvar@data
rna.aver <- AverageExpression(RNA,slot = 'counts',assays = c('chromvar'), return.seurat = T)

saveRDS(rna.aver, paste0('out/Averaged_chromvar_by_Cancer_cell_type_Pancan_700cellsPerCluster_200PCs.chromVar.20220221.', format(Sys.Date(),  format="%Y%m%d"), '.rds'))


