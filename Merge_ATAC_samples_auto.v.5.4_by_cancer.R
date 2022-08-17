## create cancer specific objects

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

require(magrittr)
require(readr)
suppressMessages(library(Matrix))
suppressMessages(library(tidyr))
set.seed(1234)

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg38)

library(googlesheets4)
library(stringr)
suppressMessages(library(doParallel))


############## FUNCTIONS #####################
getFeatureMatrix <- function (obj, peaks) {
  frag <- Fragments(obj@assays$peaksinters)
  cat('Making a large count matrix...\n')
  matrix.counts <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    process_n = 3000,
    sep = c("-","-"),
    cells = colnames(obj)
  )
  return(matrix.counts)
}

###################################
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
              default="foo", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-a", "--assay"),
              type="character",
              default="X500peaksMACS2", 
              help="which assay should be used to merge objects? X500peaksMACS2, peaks",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
assay.towork <- opt$assay

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 100 * 1024^3) # for 100 Gb RAM

###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::select(`Disease Type`, Sample, `Data Type`, `Data folder`)
samples$`Disease Type`[samples$`Disease Type`=='PKD'] <- 'ccRCC'
total.samples.id <- samples$Sample %>% as.character()
cancers <- unique(samples$`Disease Type`)

for (cancer in cancers) {
  print(cancer)
  samples.cancer <- base::subset(samples, `Disease Type`==cancer)
  samples.id <- samples.cancer$Sample %>% as.character()
  samples.type <- samples.cancer$`Data Type` %>% as.character()
  cat (paste("Samples found:" ,length(samples.id), '\n'))
  if (file.exists(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,".rds"))) {
    cat('opening normalized new peaks object...\n')
    combined <- readRDS(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,".rds"))
  } else {
    if (file.exists(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_not_normalized_",add_filename,".rds"))) {
      cat('opening not normalized new peaks object...\n')
      combined <- readRDS(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_not_normalized_",add_filename,".rds"))
    } else {
      if (file.exists(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))) {
        cat('opening not normalized old peaks object...\n')
        combined <- readRDS(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))
      } else {
        paths <- NULL
        for (i in 1:length(samples.id)){
          print(samples.id[i])
          p <- list.files(path = input.path, full.names = T, pattern = paste0(str_split_fixed(samples.id[i], '_',2)[2],'.*rds'), all.files = T, recursive = T)
          print(length(p))
          paths <- c(paths, p)
        }
        #stop if not all samples have RDS object
        print(length(samples.id))
        print(length(paths))
        stopifnot(length(samples.id)==length(paths))
        
        # make the list of atac objects
        registerDoParallel(cores=10)
        cat ('Reading in objects\n')
        #atac=vector(mode = "list", length = length(samples.id))
        atac <- foreach (i=1:length(samples.id), p = paths, .combine=c) %dopar% {
          print(samples.id[i])
          obj=readRDS(p)
          DefaultAssay(obj) <- assay.towork
          if (!file.exists(Fragments(obj)[[1]]@path)) stop("Urgh, this sample object can't locate fragments file")
          obj<- DietSeurat(obj, assay = assay.towork)
          obj$dataset = samples.id[i]
          obj$Data.type = samples.type[i]
          return(obj)
        }
        stopImplicitCluster()
        
        cat ('Reducing peaks\n')
        combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
        peakwidths <- width(combined.peaks)
        combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
        combined.peaks
        combined.peaks=combined.peaks
        combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
        combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
        #peaks.use <- combined.peaks
        peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)
        
        registerDoParallel(cores=10)
        cat ('creating matrix counts\n')
        #matrix.counts=vector(mode = "list", length = length(samples.id))
        matrix.counts <- foreach (obj = atac, .combine=c) %dopar% {
          FeatureMatrix(
            fragments = Fragments(obj@assays$X500peaksMACS2),
            features = peaks.use,
            sep = c("-","-"),
            cells = colnames(obj)
          ) 
        }
        stopImplicitCluster()
        
        checking.n.cells <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
          return(ncol(obj)==ncol(co))
          #stopifnot(ncol(obj)==ncol(co))
        }
        names(checking.n.cells) <- samples.id
        #print(checking.n.cells)
        print(checking.n.cells[!checking.n.cells])
        #stopImplicitCluster()
        #str(atac)
    
        registerDoParallel(cores=10)
        cat ('creating peaksinters and removing useless assays\n')
        atac <- foreach (obj = atac, co = matrix.counts, .combine=c) %dopar% {
          obj[['peaksinters']] <- CreateChromatinAssay(counts = co,
                                                       fragments=Fragments(obj@assays$X500peaksMACS2), 
                                                       min.cells = -1, min.features = -1)
          #obj$dataset=samples.id[i]
          DefaultAssay(obj)<-'peaksinters'
          ###remove other assay
          obj[['X500peaksMACS2']]<-NULL
          return(obj)
        }
        stopImplicitCluster()
    
        ####Merging and reduction on old peaks
        cat ('Merging\n')
        combined <- merge(x = atac[[1]], y = atac[2:length(samples.id)], add.cell.ids = samples.id)
        DefaultAssay(combined) <- "peaksinters"
        
        #remove individual objects
        rm(atac)
        gc()
        
        cat('saving the object...\n')
        saveRDS(combined, paste0(cancer,'_', length(samples.id),"_snATAC_Merged_not_normalized_",add_filename,".rds"))
      }
      
      #my peaks predone by Merge_ATAC_samples_auto_v.5.3 script
      all_peaks.filtered <- fread( paste0('peaks/',length(total.samples.id),'_',cancer, '_recentered_final.reproducible.filtered.',add_filename,'.tsv'), data.table = T)
      
      recentered_p=StringToGRanges(unique(all_peaks.filtered$new_peak), sep = c("-", "-"))
      plan("multicore", workers = 30)
      options(future.globals.maxSize = 500 * 1024^3) # for 500 Gb RAM
      
      matrix.counts <- getFeatureMatrix(combined, recentered_p)
      frag <- Fragments(combined@assays$peaksinters)
      
      # extract gene annotations from EnsDb - on congaree fails with database malfunciton WTH??? - didnt fail on lupe
      if (file.exists('Annotations.EnsDb.Hsapiens.v86.rds')) {
        annotations <- readRDS('Annotations.EnsDb.Hsapiens.v86.rds')
      } else {
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,  standard.chromosomes = TRUE)
        seqlevelsStyle(annotations) <- 'UCSC'
        genome(annotations) <- "hg38"
        saveRDS(annotations, 'Annotations.EnsDb.Hsapiens.v86.rds')
      }
      
      cat('Creating chromatin assay...\n')
      combined[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                                           annotation = annotations,
                                                           genome = 'hg38',
                                                           fragments = frag)
      # remove ATAC assay
      DefaultAssay(combined)<-'X500peaksMACS2'
      combined[['peaksinters']] <- NULL
      cat('saving the object...\n')
      saveRDS(combined, paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_not_normalized_",add_filename,".rds"))
    }
    
    #run normalization
    cat('Normalizing...\n')
    combined <- RunTFIDF(combined)
    combined <- FindTopFeatures(combined, min.cutoff = 20)
    combined <- RunSVD(
      combined,
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400
    )
    
    combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
    combined <- FindNeighbors(
      object = combined,
      reduction = 'lsi',
      dims = 2:30
    )
    
    combined <- FindClusters( 
      object = combined,
      algorithm = 3,
      resolution = 1,
      verbose = FALSE
    )
    
    cat('saving the object with updated metadata...\n')
    saveRDS(combined, paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,".rds"))
    
  }
  #update metadata
  metadata <- fread(paste0('/diskmnt/Projects/snATAC_primary/04_celltyped_rds/v4.0/',cancer, '_',length(samples.id) ,'_samples_metadata_data_freeze_v2.0.tsv')) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F) %>%
    dplyr::rename(seurat_clusters_indiv = seurat_clusters)
  combined@meta.data <- combined@meta.data[,c("orig.ident",'dataset', 'nCount_X500peaksMACS2','nFeature_X500peaksMACS2', 'seurat_clusters')]
  combined <- AddMetaData(object = combined, metadata = metadata)
  
  total_fragments_cell <- combined$passed_filters
  peak.counts <- colSums(x = GetAssayData(combined, slot='counts'))
  frip <- peak.counts *100 / total_fragments_cell
  combined <- AddMetaData(object = combined, metadata = frip, col.name = 'pct_read_in_peaks_500MACS2')
  combined <- AddMetaData(object = combined, metadata = peak.counts, col.name = 'peak_RF_500MACS2')
  
  
  print(sum(is.na(combined$Sample)))
  library(RColorBrewer)
  
  cat('plotting...\n')
  p1 <- DimPlot(combined, group.by = 'Sample', pt.size = 0.1)
  ggsave(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_sample_", add_filename, ".pdf"),plot = p1,height=12,width=17, useDingbats = F)
  
  p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
  ggsave(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_dataset_", add_filename, ".pdf"),plot = p1,height=12,width=17, useDingbats = F)
  
  p2 <- DimPlot(combined, pt.size = 0.1,label=T)
  ggsave(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_cluster_", add_filename, ".pdf"), plot = p2, height=12,width=14, useDingbats = F)
 

  n <- length(unique(combined$cell_type.harmonized))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  p3 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'cell_type.harmonized', cols = col_vector)
  ggsave(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_cell_type.harm_", add_filename, ".pdf"), plot = p3,height=12,width=15, useDingbats = F)
  
  p4 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'data.type', cols = 'Paired')
  ggsave(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_data.type_", add_filename, ".pdf"), plot = p4,height=12,width=13, useDingbats = F)
  
  p5 <- DimPlot(combined, pt.size = 0.1, label=T, group.by = 'Sample_type', cols = 'Paired')
  ggsave(paste0(cancer,'_', length(samples.id),"_snATAC_Merged_Sample_type_", add_filename, ".pdf"), plot = p5,height=12,width=13, useDingbats = F)
  
  write.table(combined@meta.data, paste0(cancer,'_', length(samples.id),"_snATAC_Merged_new_peaks_normalized_",add_filename,"_metaData.txt"),sep="\t",quote=FALSE, row.names = T)
  write.table(samples.id,paste0(cancer,'_', "Samples_snATAC_Merged_",add_filename,".txt"),sep="\t",quote=FALSE)
  
}


