#Map any object on any object


suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))

###options###
######################
option_list = list(
  make_option(c("-r", "--input.reference"),
              type="character",
              default=NULL, 
              help="path to integrated RNA object",
              metavar="character"),
  make_option(c("-q", "--input.query.object"),
              type="character",
              default=NULL, 
              help="path to reference rna object",
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
  make_option(c("-t", "--meta"),
              type="character",
              default="./", 
              help="metadata path",
              metavar="character"),
  make_option(c("-c", "--cell.column"),
              type="character",
              default="cell_type_v8.5_multi", 
              help="column with cell type to transfer",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
query.path <- opt$input.query.object
ref.path <- opt$input.reference
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$meta
cell_column <- opt$cell.column

dir.create(out_path, showWarnings = F)
setwd(out_path)


colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/Colors_panatac_v3.0.rds')
set.seed(666)
query.obj <- readRDS(query.path)
ref.obj <- readRDS(ref.path)

meta <- fread(meta.path) %>%
  column_to_rownames(var = 'V1') %>%    
  select(all_of(cell_column))
ref.obj <- AddMetaData(ref.obj, meta)


DefaultAssay(ref.obj) <- 'integrated'
ref.obj <- RunUMAP(ref.obj, reduction = "pca", dims=1:50, return.model = TRUE)

DefaultAssay(query.obj) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = ref.obj,
  query = query.obj,
  project.query = TRUE,
  k.anchor = 20,
  k.filter = 5,
  reduction = 'rpca',
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)

options(repr.plot.width=15, repr.plot.height=7)
p1 <- DimPlot(ref.obj, group.by = cell_column,
        cells.highlight = colnames(ref.obj)[anchors@anchors[,'cell1']])

p2 <- DimPlot(query.obj, group.by = 'seurat_clusters', label = TRUE,
        cells.highlight = colnames(query.obj)[anchors@anchors[,'cell2']])
p1+p2
ggsave(glue::glue('Dimplot_anchors_on_query_ref_{add_filename}.pdf'), width = 12, height = 5)


query.obj <- MapQuery(
  anchorset = anchors,
  query = query.obj,
  reference = ref.obj,
  refdata = list(
    celltype.l1 = cell_column
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)


p1 = DimPlot(query.obj, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(ref.obj, reduction = "umap", group.by = cell_column, label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p1+p2
ggsave(glue::glue('Dimplot_predicted.celltype.l1_{add_filename}.pdf'), width = 12, height = 5)

saveRDS(query.obj, glue::glue('Mapped_object_{add_filename}.rds'))
query.obj@meta.data %>% fwrite( glue::glue('Metadata_mapped_object_{add_filename}.tsv'), sep = '\t', row.names = T)

genes.oi <- c(
  'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 
  'LEF1', 'CCR7', 'TCF7', 'SELL', 'CD27','ID3', 'FAS', 'ITGA4', 'CXCR3', 'ZEB1', 'FOXP1', 'FHIT',# naive or central memory
  'KLRD1', 'KLRF1', 'KIR2DL3', 'NCAM1', 'FCGR3A', 'KLRC2', 'IKZF3', #different NKs
  'GZMB', 'GZMA' ,'GZMK', 'GNLY','PRF1','IFNG', 'FASLG', # cytotoxic stuff
  'KLRG1', 'CX3CR1', 'TNF', 'BHLHE40', 'ZEB2', 'ID2','RORC', 'RORA', 'STAT1', # effector CD8 and CD4
  'ITGAE',  'ITGA1', 'ENTPD1', 'CXCR6',  'PRDM1','RUNX3', 'ZNF683', # tissue residency
  'PDCD1', 'TIGIT','LAG3','TOX','TOX2', 'NFATC2', 'HAVCR2','CXCL13', 'CD200R1', 'NR4A2', # exhausted
  'TCF7', 'CXCR5','SLAMF6','CD200','GNG4','IGFBP4', 'TNFRSF4','TNFRSF14', 'ICOS', # pre-exhausted
  'BTLA', 'CXCL13',  'CD200','CD40LG','BATF','ICOS', # Tfh
  'FOXP3', 'CTLA4', 'IL2RA', 'IKZF2','TNFRSF18' , 'PTPN2', 'TNFRSF9','CCL4',  #Tregs
  'RORC','IL7R','IL1R1','BCL6','KIT','IL23R','KLRB1','SLC4A10','TRAV1-2','NCR3','DUSP2', #MAIT
  'LTB', 'MAL', 'CCR6','CCR4', 'RORC','IL17A','IL17F','IL1R1', #Th2
  'CCR3','CCR4','CXCR4','GATA3','IL2', 'IL4', #Th17
  'MKI67', 'TOP2A', 'MT-CO1', 'MT-ND1',
  'IFI6', 'IFI44L', 'IFI44', 'IFNG-AS1',
  'TRDV1', 'TRDV2', 'TRGV9', 'CST7', 'KLRG1', 'CCL5', #gdT
  'GZMK', 'FXYD2', 'NUCB2', 'CD27', # DN t-cells
  'KLRC2', 'KLRF1', 'CMC1', 'TYROBP', 'NKG7' # Temra CD8
) 
DotPlot(query.obj, features = unique(genes.oi), group.by = "predicted.celltype.l1") + 
  theme(axis.text.x = element_text(size = 13, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 16)) +
  scale_size_area(limits = c(0,50), oob=scales::squish) +
  scale_color_distiller(palette = 'RdPu', direction = 1)
ggsave(glue::glue('Dotplot_marker_genes_predicted.celltype.l1_{add_filename}.pdf'), width = 25, height = 7, useDingbats = F)

DimPlot(query.obj, label = TRUE, ncol = 5,
        repel = TRUE, split.by = "predicted.celltype.l1")
ggsave(glue::glue('Dimplot_query_predicted.celltype.l1_split_by_cell_type_{add_filename}.pdf'), width = 17, height = 17)

DimPlot(query.obj, label = TRUE, 
        repel = TRUE, group.by = "predicted.celltype.l1")
ggsave(glue::glue('Dimplot_query_predicted.celltype.l1_split_by_cell_type_{add_filename}.pdf'), width = 12, height = 7)


# 
# ref.obj <- DietSeurat(ref.obj, counts = FALSE, dimreducs = "pca")
# query.obj <- DietSeurat(query.obj, counts = FALSE, dimreducs = "ref.pca")
# 
# 
# ref.obj$id <- 'reference'
# query.obj$id <- 'query'
# refquery <- merge(ref.obj, query.obj)
# refquery[["pca"]] <- merge(ref.obj[["pca"]], query.obj[["ref.pca"]])
# refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50, reduction.name = "ref.query.umap", reduction.key = "refqueryUMAP_")
# refquery@meta.data[[cell_column]][is.na(refquery@meta.data[[cell_column]])] <- refquery$predicted.celltype.l1[refquery@meta.data[[cell_column]]]
# saveRDS(refquery, glue::glue('Reference_query_new_umap_object_{add_filename}.rds'))
# 
# 
# p1 <- DimPlot(refquery,reduction = "ref.query.umap", group.by = 'id', shuffle = TRUE)
# p2 <- DimPlot(refquery,reduction = cell_column, group.by = 'id', shuffle = TRUE)
# p1+p2
# ggsave(glue::glue('Dimplot_reference_query_merged_predicted.celltype.l1_{add_filename}.pdf'), width = 12, height = 5)




