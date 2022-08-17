#Find DAPs and DEGs between tumor and normal in multiome objects
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
              default=NULL, 
              help="path to merged object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--cancer"),
              type="character",
              default=NULL, 
              help="cancer type")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
round.num <- as.numeric(opt$round)
cancer <- opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')

cat('opening object \n')
combo=readRDS(input.path)
cat('done \n')

combo$cell_type.totest <- case_when(combo$cell_type.harmonized.cancer.atac=='Tumor' &  combo$cell_type.harmonized.cancer.rna=='Tumor' ~ 'Tumor',
                                    combo$cell_type.harmonized.cancer.atac=='Normal epithelial cells' & combo$cell_type.harmonized.cancer.rna=='Normal epithelial cells' ~ 'Normal epithelial cells',
                                    combo$cell_type.harmonized.cancer.atac %in% c('Luminal mature', 'Luminal progenitor','Basal progenitor') & combo$cell_type.harmonized.cancer.rna %in% c('Luminal mature', 'Luminal progenitor','Basal progenitor') ~ 'Normal epithelial cells',
                                    combo$cell_type.harmonized.cancer.atac %in% c('Ciliated Endometrial epithelial cells', 'Secretory Endometrial epithelial cells') & combo$cell_type.harmonized.cancer.rna %in% c('Ciliated Endometrial epithelial cells', 'Secretory Endometrial epithelial cells')  ~ 'Normal epithelial cells',
                                    TRUE ~ 'Other')

if (!file.exists(paste("out/da_peaks_tumor_vs_normal_epithelial_",cancer,".tsv",sep=""))) {
  Idents(combo) <- 'cell_type.totest'
  DefaultAssay(combo) <- 'pancan'
  da_peaks <- FindMarkers(
    object = combo,
    ident.1 = 'Tumor',
    ident.2='Normal epithelial cells', 
    only.pos = FALSE,
    min.pct = 0.05,
    min.diff.pct=0,
    logfc.threshold=0,
    test.use = 'LR',
    latent.vars = 'peak_RF_pancan')
  da_peaks$peak <- rownames(da_peaks)
  write.table(da_peaks, paste("out/da_peaks_tumor_vs_normal_epithelial_",cancer,".tsv",sep=""),
              sep="\t",quote=FALSE,row.names=FALSE)
}


DefaultAssay(combo) <- 'RNA'
Idents(combo) <- 'cell_type.totest'
degs <- FindMarkers(
  object = combo,
  ident.1 = 'Tumor',
  ident.2='Normal epithelial cells', 
  only.pos = FALSE,
  min.pct = 0.05,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'nCount_RNA')
degs$gene <- rownames(degs)

write.table(degs, paste("out/degs_tumor_vs_normal_epithelial_",cancer,".tsv",sep=""),
            sep="\t",quote=FALSE,row.names=FALSE)
###ENDS HERE FOR NOW####

