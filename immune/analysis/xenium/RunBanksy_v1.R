
library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(gridExtra)
library(pals)
suppressMessages(library(tidyverse))
set.seed(1234)
suppressMessages(library(data.table))
library(cowplot)
library(optparse)


filter <- dplyr::filter
select <- dplyr::select

option_list = list(
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-l", "--lambda"),
              type="double",
              default=0.8, 
              help="lambda 0.2 to 0.8",
              metavar="numeric"),
  make_option(c("-k", "--k_geom"),
              type="double",
              default=15, 
              help="k_geom 15 or 30",
              metavar="numeric")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
out_path <- opt$output
lambda.f <- opt$lambda
k_geom.f <- opt$k_geom
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)



xenium.tb <- fread('/diskmnt/Projects/snATAC_analysis/immune/validation/xenium/Xenium.table.tsv', header = T, data.table=F) %>%
  select(Sample,Cancer, Path_full)
head(xenium.tb)


pwalk(unname(as.list(xenium.tb)), function( sample,cancer,pathf) {
  if(!file.exists(glue::glue('{sample}_Niches_banksy_lambda{lambda.f}_k_geom{k_geom.f}_res_0.1_0.2_0.5.tsv'))) {
    print(sample)
    obj <- readRDS(pathf)
    all.celltypes <- fread(glue::glue('/diskmnt/Projects/snATAC_analysis/immune/validation/xenium/{cancer}/{sample}/All_cells_cell_type_v1.tsv'), data.table=F, header = T) %>%
      column_to_rownames('V1') %>%
      mutate(cell_type = case_when(cell_type=='Myoepitelial/Normal ducts' ~ 'Myoepitelial_Normal ducts',
                                   TRUE ~ cell_type))
    obj <- AddMetaData(obj, all.celltypes)
    DefaultAssay(obj) <- 'Xenium'
    obj <- NormalizeData(obj, assay = 'Xenium')
    obj <- FindVariableFeatures(obj)
    VariableFeatures(obj) <- VariableFeatures(obj)[!grepl(":|-WT", VariableFeatures(obj))]
    obj <- ScaleData(obj)
    
    #For most modern high resolution technologies like Xenium, Visium HD, StereoSeq, MERFISH, 
    #STARmap PLUS, SeqFISH+, SlideSeq v2, and CosMx (and others), we recommend the usual defults for lambda: 
    #For cell typing, use lambda = 0.2 (as shown below, or in this vignette) and for domain segmentation, use lambda = 0.8. 
    
    obj <- RunBanksy(obj, lambda = lambda.f, verbose=TRUE, 
                     assay = 'Xenium', slot = 'data', features = 'variable',
                     k_geom = k_geom.f)
    
    obj <- RunPCA(obj, assay = 'BANKSY', features = rownames(obj), npcs = 30)
    #obj <- RunUMAP(obj, dims = 1:30)
    obj <- FindNeighbors(obj, dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.5)
    
    obj <- FindClusters(obj, resolution = 0.1)
    
    obj <- FindClusters(obj, resolution = 0.2)
    
    niche.tb <- obj@meta.data %>% 
      select(starts_with('BANKSY'))
    
    fwrite(niche.tb, glue::glue('{sample}_Niches_banksy_lambda{lambda.f}_k_geom{k_geom.f}_res_0.1_0.2_0.5.tsv'), sep='\t', col.names = T, row.names = T)
    rm(obj)
    gc()
  }
})
