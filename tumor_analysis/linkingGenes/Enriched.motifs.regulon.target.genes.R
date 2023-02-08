# Alla Karpova create pancan combo object with RNA and ATAC

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =10)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(JASPAR2020))
library(motifmatchr)
library(TFBSTools)

###options###
######################
option_list = list(
  # make_option(c("-i", "--input"),
  #             type="character",
  #             default=NULL, 
  #             help="input object",
  #             metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-c", "--cancer"),
              type="character",
              default=NULL, 
              help="which cancer type to test",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
#input_path <- opt$input
out_path <- opt$output
cancer.type <- opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)

all.links.filtered <- fread('/diskmnt/Projects/snATAC_analysis/tumor_Alla/linkingGenes_v7.0/annotate.links/tumor_only_multiome/Pancan_all.links.filtered_by_zscore_diffused_links_annotated_by_GeneHancer_scEnhancer_with_CNV.txt')

genes.gained <- all.links.filtered %>% 
  filter(cnv=='Gain' & cnv.pct > 0.25 & N.cells.with.cnv > 2000) %>%
  mutate(gene.cancer = paste(gene, Cancer, sep ='_')) %>% pull(gene.cancer) %>% unique

interesting.connections <- all.links.filtered %>% 
  filter(!(paste(gene, Cancer, sep ='_') %in% genes.gained)) %>% 
  ungroup() %>% 
  dplyr::select(-cnv, -N.cells.with.cnv, -cnv.pct) %>% distinct() #%>%
  #filter(Enhancer.scEn != "" | Enhancer.GH != "") # only include links from enhancers or promoter/enhancers


cat('opening object \n')
obj <- readRDS(glue::glue('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/Multiome/Merged_objects_PanCan_peaks/Cancer_level/out/{cancer.type}_atac.tumor_cells_multiome_obj.20230121.rds'))
DefaultAssay(obj) <- "pancan"
cat('done \n')

Idents(obj) <- 'Cancer'

obj <- RegionStats(obj, assay = 'pancan', genome = BSgenome.Hsapiens.UCSC.hg38)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)
obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

current.regulons <- fread('/diskmnt/Projects/snATAC_analysis/tumor_Alla/SCENIC_targets/v7.0/column_annot.RegulonsPlot.20230131.tsv')
current.regulons %>% head()
regulons.tf <- current.regulons$Regulon_1 %>% unique

scenic.target <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snRNA/SCENIC_results/UpdatedByFreq.0.8_regulons.20220124.rds')
names(scenic.target) <- str_replace(string = names(scenic.target), pattern = "[(][+][)]", replacement = '')


regulons.tf %>% walk(function(tf) {
  print(tf)
  tf.targets <- scenic.target[[tf]]
  length(tf.targets)
  
  peaks.to.test <- interesting.connections %>% 
    filter(gene %in% tf.targets & Cancer == cancer.type)  %>% pull(peak) %>% unique
  
  peaks.to.test.against <- interesting.connections %>% 
    filter(!(gene %in% tf.targets) & Cancer == cancer.type)  %>% pull(peak) %>% unique
  
  print(length(peaks.to.test))
  print(length(peaks.to.test.against))
  
  if(length(peaks.to.test) > 5) {
    print('do shit')
    #tryCatch({
      
      # match the overall GC content in the peak set
      meta.feature <- GetAssayData(obj, assay = "pancan", slot = "meta.features")
      peaks.matched <- MatchRegionStats(
        meta.feature = meta.feature[peaks.to.test.against, ],
        query.feature = meta.feature[peaks.to.test, ],
        n = 5000
      )
      
      enriched.motifs <- FindMotifs(
        #assay = 'pancan',
        object = obj,
        features = peaks.to.test,
        background=peaks.matched
      )
      fwrite(enriched.motifs, glue::glue('Enriched_motifs_in_links_with_{tf}_targets_in_{cancer.type}.tsv'), sep='\t')
    # }, error = function(e) {
    #   message(e)
    #   data.frame()
    # })
  }
  
})

