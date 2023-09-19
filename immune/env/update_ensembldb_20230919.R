#### Author: Fernanda Martins Rodrigues @WashU (fernanda@wustl.edu)
#### Last modified: 02/19/2023

#### Get updated Ensembl annotation for seurat's multiome script
#### Updating from EnsDb.Hsapiens.v86 to v100 (correspondent to gencode 34)

## source code: https://support.bioconductor.org/p/132230/#132231

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("AnnotationHub", "ensembldb", "AnnotationForge", "biovizBase"))

library(AnnotationHub)
library(ensembldb)
library(AnnotationForge)
library(biovizBase)

# start an AnnotationHub instance/connection.
ah <- AnnotationHub()

# /home/fernanda/.cache/R/AnnotationHub
#   does not exist, create directory? (yes/no): no
# using temporary cache /diskmnt/Projects/Users/fernanda/Software/tmp/RtmpXBEvEL/BiocFileCache
#   |======================================================================| 100%

# snapshotDate(): 2022-10-31

# query for availabel Rat Ensembl databases
EnsDb.hs <- query(ah, c("EnsDb", "Homo sapiens"))
EnsDb.hs

# AnnotationHub with 23 records
# # snapshotDate(): 2022-10-31
# # $dataprovider: Ensembl
# # $species: Homo sapiens
# # $rdataclass: EnsDb
# # additional mcols(): taxonomyid, genome, description,
# #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
# #   rdatapath, sourceurl, sourcetype
# # retrieve records with, e.g., 'object[["AH53211"]]'

#              title
#   AH53211  | Ensembl 87 EnsDb for Homo Sapiens
#   AH53715  | Ensembl 88 EnsDb for Homo Sapiens
#   AH56681  | Ensembl 89 EnsDb for Homo Sapiens
#   AH57757  | Ensembl 90 EnsDb for Homo Sapiens
#   AH60773  | Ensembl 91 EnsDb for Homo Sapiens
#   ...        ...
#   AH95744  | Ensembl 104 EnsDb for Homo sapiens
#   AH98047  | Ensembl 105 EnsDb for Homo sapiens
#   AH100643 | Ensembl 106 EnsDb for Homo sapiens
#   AH104864 | Ensembl 107 EnsDb for Homo sapiens
#   AH109336 | Ensembl 108 EnsDb for Homo sapiens


# Fetch the v100 EnsDb and put it in the cache.
EnsDb.hs.v100 <- EnsDb.hs[["AH79689"]]

# check
columns(EnsDb.hs.v100)

#  [1] "DESCRIPTION"         "ENTREZID"            "EXONID"
#  [4] "EXONIDX"             "EXONSEQEND"          "EXONSEQSTART"
#  [7] "GCCONTENT"           "GENEBIOTYPE"         "GENEID"
# [10] "GENEIDVERSION"       "GENENAME"            "GENESEQEND"
# [13] "GENESEQSTART"        "INTERPROACCESSION"   "ISCIRCULAR"
# [16] "PROTDOMEND"          "PROTDOMSTART"        "PROTEINDOMAINID"
# [19] "PROTEINDOMAINSOURCE" "PROTEINID"           "PROTEINSEQUENCE"
# [22] "SEQCOORDSYSTEM"      "SEQLENGTH"           "SEQNAME"
# [25] "SEQSTRAND"           "SYMBOL"              "TXBIOTYPE"
# [28] "TXCDSSEQEND"         "TXCDSSEQSTART"       "TXID"
# [31] "TXIDVERSION"         "TXNAME"              "TXSEQEND"
# [34] "TXSEQSTART"          "TXSUPPORTLEVEL"      "UNIPROTDB"
# [37] "UNIPROTID"           "UNIPROTMAPPINGTYPE"


# Now copy, and install the database locally.
# By doing so, you can just quickly load the library next time.
# see: ?makeEnsembldbPackage

# set working dir to store the relatively large files - change this to a directory in /diskmnt/Projects of your choice
# do not do this in your home directory

setwd("/diskmnt/Projects/snATAC_analysis/immune/conda_env_files")

# copy databse from the cache to working dir

file.copy(AnnotationHub::cache(ah["AH79689"]), "./EnsDb.v100.hs.sqlite")
# [1] TRUE

# now make it a package. Change name and email accordingly
makeEnsembldbPackage("EnsDb.v100.hs.sqlite", version="0.0.1", 
	maintainer = "Fernanda Martins Rodrigues <fernanda@wustl.edu>", 
	author = "Fernanda Martins Rodrigues <fernanda@wustl.edu>", 
	destDir=".", license="Artistic-2.0")

# Creating package in ./EnsDb.Hsapiens.v100
# [1] TRUE

# Install package in R.  
# note modifed name of created directory (not EnsDb.rat... but EnsDb.Rnorvegicus.v100)
install.packages("./EnsDb.Hsapiens.v100", type = "source", repos = NULL)

# * installing *source* package ‘EnsDb.Hsapiens.v100’ ...
# ** using staged installation
# ** R
# ** inst
# ** byte-compile and prepare package for lazy loading
# ** help
# *** installing help indices
# ** building package indices
# ** testing if installed package can be loaded from temporary location
# ** testing if installed package can be loaded from final location
# ** testing if installed package keeps a record of temporary installation path
# * DONE (EnsDb.Hsapiens.v100)

# load and test package

library(EnsDb.Hsapiens.v100)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(doParallel)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v100)





