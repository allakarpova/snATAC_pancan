library(Signac)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(optparse)

option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL,
              help="path to the integrated objected",
              metavar="character"),
  make_option(c("-o", "--outdir"),
              type="character",
              default=NULL,
              help = "output directory path",
              metavar="character"),
  make_option(c("-n","--name_out_file"),
              type="character",
              default=NULL,
              help = "name of the output object without rds part",
              metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.object
out.dir <- opt$outdir
name.file <- opt$name_out_file

##################################################################
#############Call merged object, all samples, all cells:###
##################################################################


setwd(out.dir)

atac=readRDS(input.path)
DefaultAssay='peaksinters'

cat('calling peaks\n')
pbmc=atac
peaks <- CallPeaks(
  object = pbmc,
 # group.by = "predicted.id",
  macs2.path='/diskmnt/Projects/Users/allakarpova/Tools/anaconda3/envs/signac/bin/macs2'
)
cat('done!\n')

cat('make feature matrix\n')
matrix.counts <- FeatureMatrix(
    fragments = Fragments(pbmc@assays$peaksinters),
    features = peaks,
    sep = c("-","-"),
    cells = colnames(pbmc)
)
cat('done!\n')

atac=pbmc
atac[['peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                      fragments=Fragments(atac@assays$peaksinters))
DefaultAssay(atac)<-'peaksMACS2'

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:50)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:50)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)

saveRDS(atac,paste0(name.file, '.rds'))
