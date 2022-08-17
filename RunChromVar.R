library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

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
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)


atac=readRDS(paste('/diskmnt/Projects/cptac_scratch_4/CPTAC3_GBM/snATAC_seq/signac_nadya/4.Merge_atac/',
'Merge.v.20210108_peaks100_1K/13_snATAC_GBM.v20210108.rds',sep=''))
DefaultAssay(atac)='peaksinters'
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(atac),
  pwm = pfm,
  genome = 'BSgenome.Hsapiens.UCSC.hg38',
  use.counts = FALSE
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
atac <- SetAssayData(
  object = atac,
  assay = 'peaksinters',
  slot = 'motifs',
  new.data = motif
)

atac[["peaksinters"]]

atac <- RegionStats(object = atac, genome = BSgenome.Hsapiens.UCSC.hg38)

atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(atac,"13_snATAC_GBM.chromvar.v20210108.rds")