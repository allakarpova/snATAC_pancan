#Run for BRCA and BRCA_Basal separate tests:
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
  make_option(c("-t", "--type"),
              type="character",
              default='peaks', 
              help="peaks or motifs",
              metavar="numeric"),
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
ATAC=readRDS(input.path)
cat('done \n')

atac.meta <- fread(meta.path) %>% data.frame(row.names = 1) %>%
  dplyr::select(starts_with('cell_type'))

ATAC <- AddMetaData(ATAC, atac.meta)

ATAC$Disease=ifelse(ATAC$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
                                         "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3"), "BRCA_Basal",
                    ATAC$Disease)



###Try with FindMarkers:
peak.data <- GetAssayData(object = ATAC, assay = 'pancan_s', slot = "counts")
total_fragments_cell <- ATAC$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
ATAC <- AddMetaData(object = ATAC, metadata = frip, col.name = 'frip_pancan_s')
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_pancan_s')
DefaultAssay(ATAC)='pancan_s'


ATAC$totest <- case_when(ATAC$cell_type.harmonized.cancer == 'Tumor' ~ paste(ATAC$Disease, 'Tumor', sep ='-') ,
                        grepl('Plasma|B-cells|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune', ATAC$cell_type.harmonized.cancer) ~ 'Immune',
                        grepl('Endo|Fibr|Per', ATAC$cell_type.harmonized.cancer) ~ 'Stroma',
                        TRUE ~ 'Other')
print(ATAC$totest %>% unique)
Idents(ATAC)='totest'

if (opt$type=='peaks') {
  ###do peaks
  DefaultAssay(ATAC) <- 'pancan_s'
  
  
  all_da_peaks <- FindAllMarkers (ATAC, 
                                 only.pos = FALSE,
                                 min.pct = 0.1,
                                 min.diff.pct=0,
                                 logfc.threshold=0.25,
                                 test.use = 'LR',
                                 latent.vars = 'peak_RF_pancan_s')
  
  
  fwrite(all_da_peaks, paste0("out/da_peaks_cancerTypeTumor_vs_Others&ImmuneStroma_vs_Others_", format(Sys.Date(), format="%Y%B%d"), ".tsv"),
              sep="\t",row.names=FALSE)
}

if (opt$type=='motifs') {
  ###do motifs
  DefaultAssay(ATAC) <- 'chromvar'
  
  all_da_motifs <- FindAllMarkers (ATAC, 
                                  only.pos = FALSE,
                                  min.pct = 0.1,
                                  min.diff.pct=0,
                                  logfc.threshold=0.25,
                                  test.use = 'LR',
                                  latent.vars = 'peak_RF_pancan_s')
  
  
  fwrite(all_da_motifs, paste0("out/da_motifs_cancerTypeTumor_vs_Others_ImmuneStroma_vs_Others_", format(Sys.Date(), format="%Y%B%d"), ".tsv"),
         sep="\t",row.names=FALSE)
}
###ENDS HERE FOR NOW####



