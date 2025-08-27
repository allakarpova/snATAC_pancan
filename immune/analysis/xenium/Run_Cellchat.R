#run cellchat on Xenium object
###libraries
##################


suppressMessages(library(Seurat))
suppressMessages(library(CellChat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))


################################

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
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("--metadata"),
              type="character",
              default=NULL, 
              help="path to metadata"),
  make_option(c("-c", "--cell.column"),
              type="character",
              default="cell_type_sen", 
              help="column with cell type to transfer",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata
cell_column <- opt$cell.column

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
obj <- readRDS(input.path)


library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024 ^ 3)

meta <- fread(meta.path, header = TRUE) %>%
  column_to_rownames(var = 'V1') 

meta <- meta %>%    
  select(all_of(cell_column))

print (head(meta))
obj <- AddMetaData(obj, meta)
  
ct_col  <- cell_column
species <- "human"            # or "mouse"
assay_preference <- c("Xenium","RNA","SCT")  # try these in order for expression

Idents(obj) <- obj[[ct_col, drop=TRUE]]
DefaultAssay(obj) <- "Xenium"

#d = computeCellDistance(coords) # cannot allocate vector of size 107602.6 Gb
conversion.factor = 0.2125 # 0.2125 um per 1 pixel
#spot.size = min(d)*conversion.factor
spot.size=10
sf <- data.frame(ratio = conversion.factor,   
                 tol = spot.size/2)  

coords <- GetTissueCoordinates(obj, image = 'fov')
coords <- as.matrix(coords[, 1:2, drop=FALSE])

cellchat <- createCellChat(
  object       = obj,
  group.by     = ct_col,
  assay        = DefaultAssay(obj),
  datatype     = "spatial",
  coordinates  = coords,
  spatial.factors=sf
)

if (tolower(species) == "mouse") {
  CellChatDB <- CellChatDB.mouse
} else {
  CellChatDB <- CellChatDB.human
}
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Use spatial distances
cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean", trim = 0.1,
  distance.use = TRUE,
  interaction.range = 100,# µm; try 50–150 for Xenium
  scale.distance = 8, #8 works
  contact.range=10 #Typically, 'contact.range = 10', which is a typical human cell size.
)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)


write_rds(cellchat, file = paste0(add_filename, "_cellchat.rds"))

pdf('plots.pdf', width = 7, height = 7)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

dev.off()


