suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))
#plan("multicore", workers =4)
#options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
library(readxl)

###options###
######################
option_list = list(
  
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-d", "--driver"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

driver <- opt$driver
out_path <- opt$output
dir.create(out_path, showWarnings = F)
setwd(out_path)


obj <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/PanCan_merged_obj/Pancan_225Samples_1000cellsPerCluster_TumorNormal.200PCs.chromVar.20230118.rds.gz')
atac.meta <- fread('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/All_225_samples_metadata_data_freeze_v7.0_new_columns_added.tsv') %>%
  data.frame(row.names = 1)
obj <- AddMetaData(obj, atac.meta)

exclude <-  c('CE507-C1A2', 'CE354E1-S1', 'CE357E1-S1', 'CE336E1-S1', 'CE332E1-N1', 'PM565P1-T1N1')

# add driver annotation to the object
genomic.data <- read_excel('/diskmnt/Projects/snATAC_primary/Pancan_ATAC_bulk_data_freeze.v2.0/Pancan_ATAC_driver_mutation_metatable_v7_20230221.xlsx', 
                           sheet='Driver_mutation_metatable')
genomic.data <- genomic.data[,-1]
driver.piece.anno  <- genomic.data %>% dplyr::select('Piece_ID', all_of(driver))
obj@meta.data <- obj@meta.data  %>% cbind(driver.piece.anno[match(obj@meta.data$Piece_ID,driver.piece.anno$Piece_ID), -1])

set.seed(123)
Idents(obj) <- 'Cancer_Piece_cell_type.normal'
obj <- subset(obj, downsample = 200)

cells.wt <- obj@meta.data %>% 
  filter(cell_type.normal == 'Tumor') %>%
  filter(!(Piece_ID %in% exclude)) %>%
  filter(.data[[driver]] == 'WT') %>% 
  rownames

cells.trunc <- obj@meta.data %>% 
  filter(cell_type.normal == 'Tumor') %>%
  filter(!(Piece_ID %in% exclude)) %>%
  dplyr::filter(grepl('Nonsen|Shift', .data[[driver]])) %>% 
  rownames

obj$Cancer.for.test <- case_when(obj$Cancer.new == 'PDAC_Basal' ~ 'PDAC',
                                 TRUE ~ obj$Cancer.new)
print(table(obj$Cancer.for.test))

obj$Test.column <- case_when(colnames(obj) %in% cells.trunc ~ 'Truncation',
                             colnames(obj) %in% cells.wt ~ 'WT',
                             TRUE ~ 'Ignore')

print(table(obj$Test.column))


cell_t1 <- 'Truncation'
cell_t2 <- 'WT'
Idents(obj) <- 'Test.column'

obj <- subset(obj, downsample = 10000)
print(table(obj$Test.column))


da_peaks <- FindMarkers(
  object = obj,
  ident.1 = cell_t1,
  ident.2 = cell_t2,
  only.pos = FALSE,
  min.pct = 0.05,
  min.diff.pct=0,
  logfc.threshold=0.1,
  test.use = 'LR',
  latent.vars = c('peak_RF_pancan','Cancer.for.test')
)

da_peaks$cell_t1=cell_t1
da_peaks$cell_t2=cell_t2
da_peaks$Driver=driver
da_peaks$peak=rownames(da_peaks)

fwrite(da_peaks, glue::glue("{driver}_Mut_vs_WT.Pooled.pct.0.05.fCh.0.1.20230329_downsampled200.tsv"),
       sep="\t",row.names=FALSE)



