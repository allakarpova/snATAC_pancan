suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Signac))
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(cowplot)

setwd('/Users/allakarpova/lab_Ding/work/single_cell/snATAC_pancan/test_peaks_ccrcc')

all.nadya.peaks <- unique(fread('peaks/28_snATACmerged_allPeaks.Annotated.20210712.tsv')[['peak']])
important.peaks <- unique(fread('peaks/da_peaks_vs_PT_NAT.minPct0.1.20210712.tsv')[['peak']])

non.overlap <- fread('peaks/139_ccRCC_recentered_final.v4_samples.tsv')
non.overlap.unique <- unique(non.overlap[, -c('Sample'), with = F])
non.over.reprod <- fread('peaks/139_ccRCC_recentered_final.reproducible.v4_samples.tsv')
all.pancan <- fread('peaks/139_All_recentered_final.v4_samples.tsv')

peakAnno.important <- annotatePeak(StringToGRanges(important.peaks), tssRegion=c(-1000, 100),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')
peakAnno.non.overlap <- annotatePeak(StringToGRanges(unique(non.overlap$new_peak)), tssRegion=c(-1000, 100),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')
peakAnno.all.nadya <- annotatePeak(StringToGRanges(all.nadya.peaks$peak), tssRegion=c(-1000, 100),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')
peakAnno.non.over.reprod <- annotatePeak(StringToGRanges(unique(non.over.reprod$new_peak)), tssRegion=c(-1000, 100),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')
peakAnno.all.pancan <- annotatePeak(StringToGRanges(unique(all.pancan$new_peak)), tssRegion=c(-1000, 100),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = 'org.Hs.eg.db')


intersect (important.peaks$peak, non.overlap$new_peak) %>% length

peakAnno.all.nadya.dt <- as.data.table(peakAnno.all.nadya@anno)
peakAnno.all.nadya.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]
p.all <- ggplot(data = peakAnno.all.nadya.dt, aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('All peaks ccRCC from Nadya')

peakAnno.important.signif.dt <- as.data.table(peakAnno.important.signif@anno)
peakAnno.important.signif.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]

ggplot(data = peakAnno.important.dt, aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count', show.legend = F) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important peaks in ccRCC from Nadya')
ggsave('Important_peaks_ccRCC.pdf', useDingbats = F, width = 4.5, height = 4)

peakAnno.non.overlap.dt <- as.data.table(peakAnno.non.overlap@anno)
peakAnno.non.overlap.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]
p.my.all <- ggplot(data = peakAnno.non.overlap.dt, aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('All non-overlapping peaks ccRCC\nfrom pan-atac project')


peakAnno.non.over.reprod.dt <- as.data.table(peakAnno.non.over.reprod@anno)
peakAnno.non.over.reprod.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]
p.my.rep <- ggplot(data = peakAnno.non.over.reprod.dt, aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('All non-overlapping+reproducible peaks ccRCC\nfrom pan-atac project')

p.all + p.my.all + p.my.rep

# now find how many important ccrcc peaks are in my set of peaks
overlapping.impo.all <- findOverlaps(StringToGRanges(important.peaks), StringToGRanges(non.overlap$new_peak), 
                                     minoverlap =490 )
overlapping.impo.all <-  as.data.table(overlapping.impo.all)

peakAnno.important.dt <- as.data.table(peakAnno.important@anno)
peakAnno.important.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]
peakAnno.important.dt[, peak:=paste(seqnames, start, end, sep = '-')]

p.found <- ggplot(data = peakAnno.important.dt[overlapping.impo.all$queryHits,], 
       aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important ccRCC peaks found in all\nnon-overlapping peaks\nfrom pan-atac project')


p.missed <- ggplot(data = peakAnno.important.dt[-overlapping.impo.all$queryHits,], 
       aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important ccRCC peaks NOT found in\nall non-overlapping peaks\nfrom pan-atac project')
p.found + p.missed
dim(peakAnno.important.dt[overlapping.impo.all$queryHits,])
dim(peakAnno.important.dt[-overlapping.impo.all$queryHits,])
find.vector <- as.vector(mode = 'character', x = rep('Found', nrow(peakAnno.important.dt))); find.vector[-overlapping.impo.all$queryHits] <- 'Missed'
peakAnno.important.dt$Find <- find.vector
ggplot(data = peakAnno.important.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'stack') +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("Important ccRCC peaks in ccRCC\nnon-overlapping peaks\nfrom pan-atac project")
ggsave('Bar_stack_important_nadya_in_my_ccRCC_nonover400.pdf' , useDingbats = F, width = 4.5, height = 6)

ggplot(data = peakAnno.important.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'fill') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("Important ccRCC peaks in ccRCC\nnon-overlapping peaks\nfrom pan-atac project")
ggsave('Bar_fill_important_nadya_in_my_ccRCC_nonover.pdf' , useDingbats = F, width = 4.5, height = 6)


#overlap with ccRCC non-overlapping+reproducible
unique(important.peaks) %>% length
overlapping.impo.all <- findOverlaps(StringToGRanges(important.peaks), StringToGRanges(unique(non.over.reprod$new_peak)), 
                                     minoverlap = 490)
overlapping.impo.all <-  as.data.table(overlapping.impo.all)

peakAnno.important.dt <- as.data.table(peakAnno.important@anno)
peakAnno.important.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]
peakAnno.important.dt[, peak:=paste(seqnames, start, end, sep = '-')]

p.found <- ggplot(data = peakAnno.important.dt[overlapping.impo.all$queryHits,], 
                  aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important ccRCC peaks found in all\nnon-overlapping+reproducible peaks\nfrom pan-atac project')


p.missed <- ggplot(data = peakAnno.important.dt[-overlapping.impo.all$queryHits,], 
                   aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important ccRCC peaks NOT found in\nall non-overlapping+reproducible peaks\nfrom pan-atac project')
p.found + p.missed
find.vector <- as.vector(mode = 'character', x = rep('Found', nrow(peakAnno.important.dt))); find.vector[-overlapping.impo.all$queryHits] <- 'Missed'
peakAnno.important.dt$Find <- find.vector
ggplot(data = peakAnno.important.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'stack') +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("Important ccRCC peaks in ccRCC\nnon-overlapping+reproducible\nfrom pan-atac project")
ggsave('Bar_stack_important_nadya_in_my_ccRCC_nonover_rep400.pdf' , useDingbats = F, width = 4.5, height = 6)

ggplot(data = peakAnno.important.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'fill') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("Important ccRCC peaks in ccRCC\nnon-overlapping+reproducible\nfrom pan-atac project")
ggsave('Bar_fill_important_nadya_in_my_ccRCC_nonover_rep400.pdf' , useDingbats = F, width = 4.5, height = 6)


dim(peakAnno.important.dt[overlapping.impo.all$queryHits,])
dim(peakAnno.important.dt[-overlapping.impo.all$queryHits,])

peakAnno.all.pancan.dt <- as.data.table(peakAnno.all.pancan@anno)
peakAnno.all.pancan.dt[, peak.position:= str_split_fixed(annotation, ' ', 2)[,1]]
peakAnno.all.pancan.dt[, peak:=paste(seqnames, start, end, sep = '-')]
ggplot(data = peakAnno.all.pancan.dt, 
       aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Pancan non-overlapping+reproducible peaks\nfrom pan-atac project')



# overlap with pancan peaks
overlapping.impo.all <- findOverlaps(StringToGRanges(important.peaks$peak), StringToGRanges(unique(all.pancan$new_peak)), 
                                     minoverlap = 400)
overlapping.impo.all <-  as.data.table(overlapping.impo.all)
p.found <- ggplot(data = peakAnno.important.dt[overlapping.impo.all$queryHits,], 
                  aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important ccRCC peaks found in pancan\nnon-overlapping+reproducible peaks\nfrom pan-atac project')
p.missed <- ggplot(data = peakAnno.important.dt[-overlapping.impo.all$queryHits,], 
                   aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('Important ccRCC peaks NOT found in pancan\nnon-overlapping+reproducible peaks\nfrom pan-atac project')
p.found + p.missed
dim(peakAnno.important.dt[overlapping.impo.all$queryHits,])
dim(peakAnno.important.dt[-overlapping.impo.all$queryHits,])

find.vector <- as.vector(mode = 'character', x = rep('Found', nrow(peakAnno.important.dt))); find.vector[-overlapping.impo.all$queryHits] <- 'Missed'
peakAnno.important.dt$Find <- find.vector
ggplot(data = peakAnno.important.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'stack') +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("Important ccRCC peaks in pancan\nnon-overlapping+reproducible\nfrom pan-atac project")
ggsave('Bar_stack_important_nadya_in_my_pancan_nonover_rep.pdf' , useDingbats = F, width = 4.5, height = 6)

ggplot(data = peakAnno.important.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'fill') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("Important ccRCC peaks in pancan\nnon-overlapping+reproducible\nfrom pan-atac project")
ggsave('Bar_fill_important_nadya_in_my_pancan_nonover_rep.pdf' , useDingbats = F, width = 4.5, height = 6)


# overlap all Nadya's peaks with ccRCC non-verlapping peaks
overlapping.impo.all <- findOverlaps(StringToGRanges(all.nadya.peaks$peak), StringToGRanges(unique(non.overlap$new_peak)), 
                                     minoverlap = 251)
overlapping.impo.all <-  as.data.table(overlapping.impo.all)
p.found <- ggplot(data = peakAnno.all.nadya.dt[overlapping.impo.all$queryHits,], 
                  aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("All Nadya's ccRCC peaks found in ccRCC\nnon-overlapping peaks\nfrom pan-atac project")

p.missed <- ggplot(data = peakAnno.all.nadya.dt[-overlapping.impo.all$queryHits,], 
                   aes (x=peak.position, fill = peak.position)) +
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("All Nadya's ccRCC peaks NOT found in ccRCC\nnon-overlapping peaks\nfrom pan-atac project")
p.found + p.missed

find.vector <- as.vector(mode = 'character', x = rep('Found', nrow(peakAnno.all.nadya.dt))); find.vector[-overlapping.impo.all$queryHits] <- 'Missed'
peakAnno.all.nadya.dt$Find <- find.vector
ggplot(data = peakAnno.all.nadya.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'stack') +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("All Nadya's ccRCC peaks in ccRCC\nnon-overlapping peaks\nfrom pan-atac project")
ggsave('Bar_stack_all_nadya_in_my_nonover.pdf' , useDingbats = F, width = 4.5, height = 6)

ggplot(data = peakAnno.all.nadya.dt, 
       aes (x=Find)) +
  geom_bar(aes(fill = peak.position),stat = 'count', position = 'fill') +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = 'grey')) +
  scale_fill_brewer(palette = 'Set1') +
  ggtitle("All Nadya's ccRCC peaks in ccRCC\nnon-overlapping peaks\nfrom pan-atac project")
ggsave('Bar_fill_all_nadya_in_my_nonover.pdf' , useDingbats = F, width = 4.5, height = 6)


dim(peakAnno.all.nadya.dt[overlapping.impo.all$queryHits,])
dim(peakAnno.all.nadya.dt[-overlapping.impo.all$queryHits,])


