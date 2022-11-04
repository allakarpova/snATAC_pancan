#Transfer labels from integrated RNA to integrated ATAC object
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))



###options###
######################
option_list = list(
  make_option(c("-r", "--input.rna.object"),
              type="character",
              default=NULL, 
              help="path to integrated RNA object",
              metavar="character"),
  make_option(c("-a", "--input.atac.object"),
              type="character",
              default=NULL, 
              help="path to integrated ATAC object",
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
              help = "path to RNA metadata file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

# read in initial arguments
input.path.rna <- opt$input.rna.object
input.path.atac <- opt$input.atac.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

my.metadata <- fread(meta.path, data.table = F, header=TRUE) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v5.0/Colors_panatac_v2.0.rds')
panc.rna <- readRDS(input.path.rna)
panc.atac <- readRDS(input.path.atac)
panc.rna <- AddMetaData(panc.rna, my.metadata)

plan("multicore", workers = 20)
options(future.globals.maxSize = 250 * 1024^3) # for 250 Gb RAM

print(head(VariableFeatures(panc.rna)))
length(VariableFeatures(panc.rna))

if(!file.exists(paste0('Panimmune_all_ATAC_GENEACTIVITY', add_filename, '.rds'))) {
  cat('Run Gene ACTIVITY\n')
  gene.activities <- GeneActivity(panc.atac, features = VariableFeatures(panc.rna))
  
  saveRDS(gene.activities, paste0('Panimmune_all_ATAC_GENEACTIVITY', add_filename, '.rds'))
  
} else {
  gene.activities <- readRDS(paste0('Panimmune_all_ATAC_GENEACTIVITY', add_filename, '.rds'))
}

# add gene activities as a new assay
panc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
cat('Add gene activities to ATAC object\n')
DefaultAssay(panc.atac) <- "ACTIVITY"
panc.atac <- NormalizeData(panc.atac)
panc.atac <- ScaleData(panc.atac, features = rownames(panc.atac))

DefaultAssay(panc.atac) <- "ATAC_immune"
panc.atac <-  panc.atac %>% 
  FindTopFeatures(min.cutoff = 500) %>%
  RunTFIDF() %>%
  RunSVD()


# Identify anchors
cat('Run FindTransferAnchors\n')

DefaultAssay(panc.atac) <- "ACTIVITY"

panc.rna <-  panc.rna %>% SCTransform(
  assay = 'RNA',
  vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
  conserve.memory = T,
  return.only.var.genes = T
)

transfer.anchors <- FindTransferAnchors(reference = panc.rna, 
                                        query = panc.atac, features = VariableFeatures(object = panc.rna), 
                                        normalization.method = 'SCT',
                                        reference.assay = "SCT", query.assay = "ACTIVITY", reduction = "cca")

cat('Run TransferData\n')
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = as.character(panc.rna@meta.data[,cell_column]),
                                     weight.reduction = panc.atac[["lsi"]], dims = 2:50)

panc.atac <- AddMetaData(panc.atac, metadata = celltype.predictions)

fwrite(celltype.predictions, paste0('PanImmune_all_ATAC_', add_filename, '_predicted_labels_transfered_fromRNA.txt'),  sep='\t', row.names = TRUE)

p1 <- DimPlot(panc.rna,  group.by = cell_column, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference RNA")
p2 <- DimPlot(panc.atac,  group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query ATAC")

p1 | p2
ggsave(paste0('Dimplot_RNA_labels_transferred_to_ATAC_', add_filename, '_cell_type_predicted.pdf'), width = 12, height = 4.5)






