## Alla Karpova
### recluster T-cells use only protein coding genes, remove MT- genes

suppressMessages(library(Signac))
suppressMessages(library(Seurat))

suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))


suppressMessages(library(future))
suppressMessages(library(optparse))


suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
library("glmGamPoi")


runAllNormalization <- function(obj, dims=30) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      method = "glmGamPoi",
      #vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      conserve.memory = T,
      return.only.var.genes = T
    ) %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunUMAP(dims = 1:dims, assay = 'SCT')
  
  obj <- NormalizeData(obj, assay = 'RNA')
  obj <- FindNeighbors(object = obj,  dims = 1:dims)
  obj <- FindClusters(object = obj,resolution = 2,algorithm = 4,
                      method='igraph',  verbose = FALSE)
  return(obj)
  
}


filter <- dplyr::filter
select <- dplyr::select


option_list = list(
  make_option(c("-i", "--input.obj"),
              type="character",
              default=NULL,
              help="path to input object",
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
  make_option(c("-t","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type.harmonized',
              help = "cell_type_column",
              metavar="character"),
  make_option(c("--int_batch"),
              type="character",
              default="chemistry",
              help = "options include 'weird_Brca_Ov','sample', 'cancer','chemistry', 'cancer_chemistry'",
              metavar="character"),
  make_option(c("--do.reference"),
              type="logical",
              default=FALSE,
              help = "Do reference integartion against multiome samples or not",
              metavar="logical")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.obj
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


plan("multicore", workers = 20)
options(future.globals.maxSize = 50 * 1024^3) # for 250 Gb RAM


all.rna <- readRDS(input.path)
# add meta data if provided
if (!is.null(meta.path)) {
  my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F) %>% select(all_of(cell_column))
  all.rna <- AddMetaData(all.rna, metadata = my.metadata)
}
all.rna@meta.data$to_remove <- grepl('oublet', as.character(unlist(all.rna[[cell_column]])))
DefaultAssay(all.rna) <- 'RNA'

all.rna <- all.rna %>%  DietSeurat(assay = 'RNA', counts = TRUE, data = TRUE)
all.rna <- subset(x = all.rna, subset = to_remove, invert = TRUE)
dim(all.rna)


print(opt$int_batch)
if(opt$int_batch=='weird.gbm') {
  wierd.gbm <- c('C3N-02783', 'C3L-02705', 'C3N-01334', 'C3N-01798', 'C3L-03968')
  all.rna@meta.data$Batches <- case_when(all.rna$Piece_ID_RNA %in% wierd.gbm ~  paste(all.rna$Cancer, 'weird', sep = '__'),
                                         all.rna$Cancer %in% c('PBMC') ~ paste(all.rna$Cancer, all.rna$Chemistry, sep = '__'),
                                         all.rna$Cancer %in% c('MM') ~ all.rna$Cancer,
                                         TRUE ~ all.rna$Chemistry)
  
} else if(opt$int_batch=='weird_Brca_Ov') {
  message('doing weird_Brca_Ov')
  wierd.brca <- c('HT206B1-S1H4', 'HT378B1-S1H2')
  wierd.ov <- c('VF031V1-Tm1Y1', 'VF027V1-S1Y1', 'VF034V1-T1Y1')
  all.rna@meta.data$Batches <- case_when(all.rna$Piece_ID_RNA %in% wierd.brca ~  paste(all.rna$Cancer, 'weird', sep = '__'),
                                         all.rna$Piece_ID_RNA %in% wierd.ov ~  paste(all.rna$Cancer, 'weird', sep = '__'),
                                         all.rna$Cancer %in% c('PBMC') ~ paste(all.rna$Cancer, all.rna$Chemistry, sep = '__'),
                                         all.rna$Cancer %in% c('MM') ~ all.rna$Cancer,
                                         TRUE ~ all.rna$Chemistry)
}  else if(opt$int_batch=='sample') {
  all.rna@meta.data$Batches <- all.rna@meta.data$Piece_ID_RNA
  
} else if(opt$int_batch=='cancer') {
  all.rna@meta.data$Batches <- all.rna$Cancer
} else if(opt$int_batch=='chemistry') {
  all.rna@meta.data$Batches <- case_when(all.rna$Cancer %in% c('PBMC') ~ paste(all.rna$Cancer, all.rna$Chemistry, sep = '__'),
                                         all.rna$Cancer %in% c('MM') ~ all.rna$Cancer,
                                         TRUE ~ all.rna$Chemistry)
} else if(opt$int_batch=='cancer_chemistry') {
  all.rna@meta.data$Batches <- case_when(all.rna$Cancer %in% c('GBM') ~ all.rna$Cancer,
                                         TRUE ~ paste(all.rna$Cancer, all.rna$Chemistry, sep = '__'))
}

print(table(all.rna$Batches))

cat ('Run SCT on batches\n')
all.rna.list <- SplitObject(all.rna, split.by = 'Batches')

batches <- names(all.rna.list)
multiome.batches.n <- which(grepl('Multiome', batches))

tumor.genes <- fread('/diskmnt/Projects/snATAC_analysis/immune/obj/v8.0/auxiliary/DEGs_tumor_T-cells/Pval_aver_immune_cell_pct_vs_tumor_pct.tsv', data.table = F, header = TRUE)
trash.genes <- tumor.genes %>% 
  filter(p.vals<0.1, pct.tumor.fake>0, Cancer.rank < 6) %>% pull(features.plot) %>% unique()
interferon.genes <- fread('/diskmnt/Projects/snATAC_analysis/immune/cell_typing/v8.0/allRNA2/by_lineage/no_doublets/Int_chemistry_cancer_reference_multiome2/no_doublets2/Interferon.reponse.genes.tsv', data.table = F, header = T) %>%
  pull(V1)

all_coding_genes <- fread('/diskmnt/Projects/snATAC_analysis/immune/gene_sets/Coding_genes.tsv', header = TRUE)
head(all_coding_genes)
mito.genes <- rownames(all.rna)[grepl('^MT-', rownames(all.rna))]
genes.i.dont.want2 <- Reduce(union, list(mito.genes, trash.genes))
#genes.i.dont.want2 <- trash.genes
cellcycle.genes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- CellCycleScoring(x, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  
  print('boop')
  x <- x %>% SCTransform(
    assay = 'RNA',
    method = "glmGamPoi",
    min_cells=1,
    # vars.to.regress =  regress.me,
    conserve.memory = F,
    verbose = F,
    return.only.var.genes = T
  )
  
  good.variable.genes <- setdiff(rownames(x@assays$SCT@scale.data),genes.i.dont.want2)
  good.variable.genes <- intersect(all_coding_genes$hgnc_symbol, good.variable.genes) # retain only protein coding genes
  x@assays$SCT@scale.data <- x@assays$SCT@scale.data[good.variable.genes,]
  VariableFeatures(x) <- good.variable.genes
  #x <- GetResidual(x, features = cellcycle.genes)
  return(x)
})

message('Selecting integration features')
features <- SelectIntegrationFeatures(object.list = all.rna.list, nfeatures = 4500)
print(length(features))


all.rna.list <- PrepSCTIntegration(object.list = all.rna.list, anchor.features = features)
message('Run PCA on integration features')
all.rna.list <- lapply(X = all.rna.list, FUN = RunPCA, features = features, verbose=F)
message('Run FindIntegrationAnchors')
if(opt$do.reference) {
  rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, reference = multiome.batches.n,
                                        normalization.method = "SCT",
                                        anchor.features = features, dims = 1:50, reduction = "rpca")
  
} else {
  rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:50, reduction = "rpca")
  
}
message('Run IntegrateData')

int <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:50)
DefaultAssay(int) <- 'integrated'
#int <- CellCycleScoring(int, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
message('Run ScaleData on integrated assay')
int <- ScaleData(int, 
                 #vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
                 do.scale = FALSE,
                 do.center = TRUE)


int <- RunPCA(int, verbose = FALSE, npcs = 70)
int <- RunUMAP(int, reduction = "pca", dims = 1:60)
int <- FindNeighbors(int, reduction = "pca", dims = 1:60)
int <- FindClusters(int, resolution = 1.3, algorithm = 4,
                    method='igraph')

int <- PrepSCTFindMarkers(int)
VariableFeatures(int) <- features

cat('saving the object...\n')
saveRDS(int,   paste0("PanStroma_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".rds"))

dimplot=DimPlot(object = int, label = TRUE) + NoLegend()
pdf(paste0("Dimplot_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "Piece_ID_RNA")
pdf(paste0("Dimplot_Piece_ID_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".pdf"),height=12,width=25, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = F, group.by = "data.type.rna")
pdf(paste0("Dimplot_data.type_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".pdf"),height=12,width=14, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "cell_type.harmonized.cancer")
pdf(paste0("Dimplot_cell_type.harmonized.cancer_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = cell_column)
pdf(paste0("Dimplot_",cell_column,"_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

dimplot=DimPlot(object = int, label = TRUE, group.by = "Cancer")
pdf(paste0("Dimplot_Cancer_integrated_allRNA_by_",opt$int_batch,"_",add_filename,".pdf"),height=12,width=15, useDingbats = F)
print(dimplot)
dev.off()

