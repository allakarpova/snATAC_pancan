####2021-05-04
#Modified by Alla Karpova to include RNA processing from the combo data

#Modified by Nadezhda Terekhanova (nvterekhanova@wustl.edu) using updated peak calling with MACS2

#Mofidified on 08/31/20 by Nadezhda Terekhanova using updates from Signac v.1.0.0, see here changes compared to Signac v.0.2.5: https://cran.r-project.org/web/packages/Signac/news/news.html

#Reyka Jayasinghe (reyka@wustl.edu)
#Modified on 08/27/20

#References based on 
#https://satijalab.org/signac/articles/pbmc_vignette.html

library(optparse)
set.seed(1234)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

option_list = list(
 make_option(c("-s", "--sample"),
	type="character", 
	default=NULL, 
	help = "sample_name", 
	metavar="character"),
 make_option(c("-a","--atac_data"),
	type="character",
	default=NULL,
	help = "path to ATAC data folder (e.g. cellranger output's raw matrices folder)",
	metavar="character"),
 make_option(c("-r","--rna_data"),
             type="character",
             default=NULL,
             help = "path to RNA data folder (e.g. cellranger output's raw matrices folder)",
             metavar="character"),
 make_option(c("-m","--macs2_path"),
	type="character",
	default=NULL,
	help = "path to installed MACS2",
	metavar="character"),
 make_option(c("-o","--output_folder"),
             type="character",
             default=NULL,
             help = "output folder where a sample subfolder will be created",
             metavar="character"),
#CellRanger ATAC QC metrics
 make_option(c("--prf_min"),
	type="integer", 
	default=3000, 
	help = "peak_region_fragments_minimum value for filtering", 
	metavar="integer"),
 make_option(c("--prf_max"),
	type="integer", 
	default=20000, 
	help = "peak_region_fragments_maximum value for filtering", 
	metavar="integer"),
 make_option(c("--pct_min"),
	type="integer", 
	default=15, 
	help = "pct_reads_in_peaks_minimum value for filtering", 
	metavar="integer"),
 make_option(c("--bl_ratio"),
	type="double", 
	default=0.05, 
	help = "blacklist_ratio_minimum value for filtering", 
	metavar="double"),
#Changed to default=4, based on the latest Signac-vignette
 make_option(c("--ns_max"),
	type="integer", 
	default=4, 
	help = "nucleosome_signal_maximum value for filtering", 
	metavar="integer"),
 make_option(c("--tss"),
	type="integer", 
	default=2, 
	help = "tss_enrichment_minimum value for filtering", 
	metavar="integer"),
 make_option(c("--pc_num"),
	type="integer", 
	default=30, 
	help = "number of principal components to use", 
	metavar="integer"),
 make_option(c("--pc_first"),
	type="integer", 
	default=1, 
	help = "first principal components to use (should be 1 or 2)", 
	metavar="integer"),
#### RNA QC metrics
make_option(c("--pre_filter"),
            type="integer",
            default=300,
            help="min number of reads per cell to prefilter",
            metavar="integer"),
make_option(c("--nfeature_min"),
            type="integer",
            default=200,
            help="nFeature_RNA min value for filtering",
            metavar="integer"),
make_option(c("--nfeature_max"),
            type="integer",
            default=10000,
            help="nFeature_RNA max value for filtering",
            metavar="integer"),
make_option(c("--ncount_min"),
            type="integer",
            default=1000,
            help="nCount_RNA min value for filtering",
            metavar="integer"),
make_option(c("--ncount_max"),
            type="integer",
            default=80000,
            help="nCount_RNA max value for filtering",
            metavar="integer"),
make_option(c("--mito_max"),
            type="double",
            default=0.1,
            help="maximum allowed mitochondrial fraction",
            metavar="double")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (sample_name,atac_data).n", call.=FALSE)
}

print("Input parameters")
print(paste("--peak_region_fragments_min:",opt$prf_min,sep=""))
print(paste("--peak_region_fragments_max:",opt$prf_max,sep=""))
print(paste("--pct_reads_in_peaks_minimum:",opt$pct_min,sep=""))
print(paste("--blacklist_ratio_minimum:",opt$bl_ratio,sep=""))
print(paste("--nucleosome_signal_maximum:",opt$ns_max,sep=""))
print(paste("--tss_enrichment_minimum:",opt$tss,sep=""))
print(paste("--pc_first:",opt$pc_first,sep=""))
print(paste("--pc_num:",opt$pc_num,sep=""))

##input data
sample=opt$sample
atac_data_folder=opt$a
rna_data=opt$r

print(paste("ATAC data:",atac_data_folder,sep=""))
print(paste("RNA data:",rna_data,sep=""))

##output data
print(samples)

#####LOAD REQUIRED PACKAGES##########
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)



	print(sample)

	###########################
	########LOAD IN DATA#######
	###########################
	
	dir.create("out")
	outputpath=paste("out/",sample,"/",sep="")
	dir.create(outputpath)

	counts <- Read10X_h5(paste(atac_data_folder,"/",sample,"/outs/raw_feature_bc_matrix.h5",sep=""))
#	metadata <-read.csv(file=paste(atac_data_folder,"/",sample,"/outs/singlecell.csv",sep=""), header = TRUE, row.names = 1)
	fragment.path <- paste(atac_data_folder,"/",sample,"/outs/atac_fragments.tsv.gz",sep="")

	#remove min.cells filter - returns the error, can filter peaks in the downstreaam analysis
	#Create Seurat Object
	chrom_assay <- CreateChromatinAssay(
 	  counts = counts$Peaks,
  	  sep = c(":", "-"),
  	  genome = 'hg38',
  	  fragments = fragment.path,
  	  min.cells = -1,
  	  min.features = 200
	)


	pbmc <- CreateSeuratObject(
	  counts = chrom_assay,
	  assay = "peaks",
	  project = 'ATAC',
#	  meta.data = metadata
	)

####2021-03-20: Change Peaks to MACS2
	###########################################################
	############MACS2 peak calling#############################
	###########################################################
	peaks <- CallPeaks(
 	 object = pbmc,
 	 # group.by = "predicted.id",
  	 macs2.path=opt$m
	 )
	 
	 p=as.data.frame(peaks)
	 write.table(p,paste(outputpath,'MACS2_peaks.',sample,'.tsv',sep=''),sep='\t',quote=FALSE,
row.names=FALSE)

	 p$peak_center=p$start+p$relative_summit_position
	 p$recentered_start=p$peak_center-250
	 p$recentered_end=p$peak_center+250

	 ####Now check that new start and end don't go beyond the chromosome boundaries
	 chr_size=read.table('hg38.chrom.sizes.txt',sep='\t',header=FALSE)
	 colnames(chr_size)=c('seqnames','chr_length')
	 p1=merge(p,chr_size,all.x=TRUE)

	 p1=p1[p1$recentered_end<=p1$chr_length && p1$recentered_start>=0,]
	 p1$length=p1$recentered_end-p1$recentered_start+1
	 p1$new_peak=paste(p1$seqnames,p1$recentered_start,p1$recentered_end,sep='-')

	 recentered_p=StringToGRanges(p1$new_peak, sep = c("-", "-"))

	 olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
	 olap1=olap[olap$queryHits!=olap$subjectHits,]

	 recentered_non_olap=p1[-olap1$queryHits,]
	 recentered_olap=p1[olap1$queryHits,]

	 pairs=cbind(p1[olap1$queryHits,c(1:3,7)],olap1$queryHits,p1[olap1$subjectHits,c(1:3,7)],olap1$subjectHits)
	 colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

	 pairs=pairs[pairs$score_1>=pairs$score_2,]
	 all_st=NULL
	 for (i in 1:nrow(pairs)){
	     if (nrow(pairs)>0){
	         pairs=pairs[order(-pairs$score_1),]
		 if (pairs[1,4]>=pairs[1,9]){
       		    all_st=rbind(all_st,p1[rownames(p1)==pairs[1,5],])
		    pairs=pairs[-1,]
		    pairs=pairs[pairs$row_1!=pairs[1,10],]
       		 }else{
       		    all_st=rbind(all_st,p1[rownames(p1)==pairs[1,10],])
       		    pairs=pairs[-1,]
       		    pairs=pairs[pairs$row_1!=pairs[1,5],]
       		 }
	     }
	}
	all_st=as.data.frame(all_st)
	all_st=all_st[!duplicated(all_st),]

	recentered_final=rbind(recentered_non_olap,all_st)
	write.table(recentered_final,paste(outputpath,'recentered_final.filtered',sample,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)

	recentered_p=StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
	matrix.counts <- FeatureMatrix(
    	   fragments = Fragments(pbmc@assays$peaks),
    	   features = recentered_p,
    	   sep = c("-","-"),
    	   cells = colnames(pbmc)
	)

	pbmc[['X500peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
	fragments=Fragments(pbmc@assays$peaks))
	DefaultAssay(pbmc)<-'X500peaksMACS2'

	peak.data <- GetAssayData(object = pbmc, assay = 'X500peaksMACS2', slot = "counts")
	total_fragments_cell <- pbmc$passed_filters
	peak.counts <- colSums(x = peak.data)
	frip <- peak.counts / total_fragments_cell
	pbmc <- AddMetaData(object = pbmc, metadata = frip, col.name = 'frip_500MACS2')
	pbmc <- AddMetaData(object = pbmc, metadata = peak.counts, col.name = 'peak_RF_500MACS2')

	###########################################################
	############Adding annotations for hg38 to the object######
	###########################################################

	# extract gene annotations from EnsDb
	annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

	# change to UCSC style since the data was mapped to hg19
	seqlevelsStyle(annotations) <- 'UCSC'
	genome(annotations) <- "hg38"

	# add the gene information to the object
	Annotation(pbmc) <- annotations	


	###########################################################
	############Quality Control OF SC-ATAC DATA################
	###########################################################
	#https://satijalab.org/signac/articles/pbmc_vignette.html

	# compute nucleosome signal score per cell
	pbmc <- NucleosomeSignal(object = pbmc)

	# compute TSS enrichment score per cell
	pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

	# add blacklist ratio and fraction of reads in peaks
	pbmc$pct_reads_in_peaks <- pbmc$peak_RF_500MACS2 / pbmc$passed_filters * 100
	pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_RF_500MACS2

	# inspecting TSS-enrichment scores
	pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
	tss_plot=TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

	# inspecting fragment length periodicity
	pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > opt$ns_max, 'NS > opt$ns_max', 'NS < opt$ns_max')
	fragment_period_plot=FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')	
	
	QC_plot=VlnPlot(
 	  object = pbmc,
	  features = c('pct_reads_in_peaks', 'peak_RF_500MACS2',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  	  pt.size = 0.1,
          ncol = 5
	)	

	pdf(paste(outputpath,"/",sample,"_0_QC.pdf",sep=""),height=6,width=12)
	print(tss_plot)
	print(fragment_period_plot)
	print(QC_plot)
	dev.off()

	#remove cells that are outliers for these QC metrics
	pbmc <- subset(
	x = pbmc,
	subset = peak_RF_500MACS2 > opt$prf_min &
	peak_RF_500MACS2 < opt$prf_max &
	pct_reads_in_peaks > opt$pct_min &
	blacklist_ratio < opt$bl_ratio &
	nucleosome_signal < opt$ns_max &
	TSS.enrichment > opt$tss
)

	##################################################
	##Normalization and linear dimensional reduction##
	##################################################

	pbmc <- RunTFIDF(pbmc)
	pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
	pbmc <- RunSVD(pbmc)

	#Check if first LSI-component correlated with the sequencibg depth. If it is, then re-run using LSI components starting from 2 (for exaample, 2:30 instead of 1:30)
	depth_corr_plot=DepthCor(pbmc)
	pdf(paste(outputpath,"/",sample,"_DepthCorrelation_1_QC.pdf",sep=""),height=6,width=12)
	print(depth_corr_plot)
	dev.off()
	
	##################################################
	##Non-linear dimension reduction and clustering###
	##################################################

	# perform graph-based clustering and non-linear dimension reduction for visualization
	pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = opt$pc_first:opt$pc_num)
	pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = opt$pc_first:opt$pc_num)
	pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
	dimplot=DimPlot(object = pbmc, label = TRUE) + NoLegend()

	pdf(paste(outputpath,"/",sample,"_2_Dimplot.pdf",sep=""),height=6,width=6)
	print(dimplot)
	dev.off()

	
	##################################################
	############Create a gene activity matrix#########
	##################################################

	#extract gene coordinates and extend them to include the 2 kb upstream region (as promoter accessibility is often correlated with gene expression)
	gene.activities <- GeneActivity(pbmc)
	
	# add the gene activity matrix to the Seurat object as a new assay and normalize it
	pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
	pbmc <- NormalizeData(
  	  object = pbmc,
	  assay = 'RNA',
	  normalization.method = 'LogNormalize',
	 scale.factor = median(pbmc$nCount_RNA)
	)

	DefaultAssay(pbmc) <- 'RNA'
	
	#Save object
	saveRDS(pbmc,file = paste(outputpath,"/",sample, "_processed_atac.rds", sep=""))


