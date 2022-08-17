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
              default='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/10.snRNA/1.Merge_snRNA/Merge.v.20220307/Merged_132_snRNA.SCTrandformed.3KVarFeatures.v3.20220307.rds', 
              help="path to merged object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
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

DefaultAssay(RNA)='SCT'

RNA$totest <- case_when(RNA$cell_type.harmonized.cancer == 'Tumor' ~ paste(RNA$Cancer, 'Tumor', sep ='-') ,
                         grepl('Plasma|B-cells|T-cell|DC|Macro|Mast|Microglia|NK|Treg|Immune', RNA$cell_type.harmonized.cancer) ~ 'Immune',
                         grepl('Endo|Fibr|Per', RNA$cell_type.harmonized.cancer) ~ 'Stroma',
                         TRUE ~ 'Other')
print(RNA$totest %>% unique)
Idents(RNA)='totest'
  
all_deg <- FindAllMarkers (RNA, 
                            only.pos = FALSE,
                            min.pct = 0.1,
                            min.diff.pct=0,
                            logfc.threshold=0.25,
                            test.use = 'LR',
                            latent.vars = 'nCount_RNA')


fwrite(all_deg, paste0("out/degs_cancerTypeTumor_vs_Others&ImmuneStroma_vs_Others_", format(Sys.Date(), format="%Y%B%d"), ".tsv"),
       sep="\t",row.names=FALSE)

