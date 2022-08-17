#####author: Nadezhda V. Terekhanova
#####2020-09-26: v.3.0 - this version of the script was adapted for the latest Signac v.1.0.0
###do prior to running the script `export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp` - otherwise it crushes
### Alla Karpova adapted this script. works on signac 1.1.0 too

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)
library(optparse)
library(data.table)

option_list = list(
  make_option(c("-i", "--input.folder"),
              type="character",
              default=NULL, 
              help="path to folder with rds objects",
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
  make_option(c("-a", "--assay"),
              type="character",
              default="peaks", 
              help="which assay should be used to merge objects? X500peaksMACS2, peaks",
              metavar="character"),
  make_option(c("-s", "--samples.file"),
              type="character",
              default=NULL, 
              help="path to file with a list of samples in one column, no header",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay
sample.path <- opt$samples.file

dir.create(out_path, showWarnings = F)
setwd(out_path)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2) # for 50 Gb RAM

###### CHANGE THIS!!!!!! ########
samples <- fread(sample.path, data.table = F, header = F)
samples <- samples$V1 %>% as.character()
#########
cat (paste("Samples found:" ,length(samples), '\n'))

# make the list of atac objects
atac=vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
  atac[[i]]=readRDS(list.files(path = input.path, full.names = T, pattern = paste0(samples[i],'.*rds'), all.files = T, recursive = T))
  DefaultAssay(atac[[i]]) <- assay.towork
}


#####To obtain the best results - use ALL peaks!

combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
peaks.use=combined.peaks

#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.

matrix.counts=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
  matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$peaks),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
  ) 
}


for (i in 1:length(samples)){
  atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],fragments=Fragments(atac[[i]]@assays$peaks))
  atac[[i]]$dataset=samples[i]
  DefaultAssay(atac[[i]])<-'peaksinters'
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)

DefaultAssay(combined) <- "peaksinters"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(
  combined,
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work = 400
)

combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
combined <- FindNeighbors(
  object = combined,
  reduction = 'lsi',
  dims = 2:50
)
combined <- FindClusters( 
  object = combined,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

cat('saving the object...\n')
saveRDS(combined, paste0(length(samples),"_snATAC_Merged_",add_filename,".rds"))

cat('plotting...\n')
p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1) + 
  ggplot2::ggtitle("Combined snATAC samples")+
  scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),"grey"))


p2 <- DimPlot(combined, pt.size = 0.1,label=T)

pdf(paste0(length(samples),"_snATAC_Merged_", add_filename, ".pdf"),height=6,width=16, useDingbats = F)
p1+p2
dev.off()

write.table(samples,paste0("Samples_snATAC_Merged_",add_filename,".txt"),sep="\t",quote=FALSE)

