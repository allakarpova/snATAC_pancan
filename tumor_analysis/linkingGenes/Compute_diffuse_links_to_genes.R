# Alla Karpova create pancan combo object with RNA and ATAC

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))


###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="input object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')


cat('opening object \n')
obj <- readRDS(input_path)
DefaultAssay(obj) <- "pancan"
cat('done \n')

annot <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
Annotation(obj) <- annot


chr.size=fread('/diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt', header = FALSE,col.names = c('seqnames','chr_length') ) %>% 
  data.frame()
chr.size <- chr.size %>% filter(seqnames %in% paste0('chr', 1:22) | seqnames == 'chrX')
chr.size.vector <- as.numeric(chr.size$chr_length)
names(chr.size.vector) <- chr.size$seqnames
bins <- GenomicRanges::tileGenome(seqlengths = chr.size.vector, tilewidth = 100000, cut.last.tile.in.chrom = TRUE)

frag <- Fragments(obj@assays$pancan)
#this will remove fragment objects with just 1 or 0 cells because they fail FeatureMatrix function
frag.filtered <- frag[do.call( 'c', lapply(frag, function(x) length(x@cells) > 1))]

# remove pancan assay
DefaultAssay(obj)<-'SCT'
obj[['pancan']] <- NULL

plan("multicore", workers = 30)
options(future.globals.maxSize = 300 * 1024^3) # for 500 Gb RAM

print(length(bins))
pro_n <-  round(length(bins)/30)


exlude.samples <- obj@meta.data %>% group_by(Sample) %>% tally() %>% filter(n<=1) %>% pull(Sample)
#now exclude these samples from the object
obj <- subset(obj, subset = Sample %in% exlude.samples, invert=TRUE)


matrix.counts <- FeatureMatrix(
  fragments = frag.filtered,
  features = bins,
  process_n = pro_n,
  sep = c("-","-"),
  cells = colnames(obj)
)

obj[['X100kb']] <- CreateChromatinAssay(counts = matrix.counts,
                                                 annotation = annot,
                                                 genome = 'hg38',
                                                 fragments = frag.filtered, 
                                                  min.features = -1)

DefaultAssay(obj)<-'X100kb'


# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "X100kb",
  distance = 5e+05,
  score_cutoff = 0,
  n_sample = 1000,
  expression.assay = "SCT"
)

toreturn <- Links(obj)
out.obj <- str_replace(input_path, pattern = 'rds', replacement = 'diffuse.links.only.rds')
out.obj <- str_split(out.obj, '[/]')[[1]][str_count(out.obj, '/')+1]

print(out.obj)
saveRDS(toreturn, paste0('out/', out.obj))
