suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))


option_list = list(
  make_option(c( "--input"),
              type="character",
              default=NULL,
              help="path to data folder ",
              metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments

matrix.path = opt$input

all.counts <- list.files(matrix.path, full.names = T) %>%
    map(function(p) {
      
        file.name <- basename(p)#
        
        clean.file.name <- str_remove(file.name, '_1')
        clean.file.name <- str_replace(clean.file.name, '_', '-')
        case <- str_split_fixed(clean.file.name, '-', 2)[,1]
        sample <- str_split_fixed(clean.file.name, '[.]', 2)[,1]
        print(sample)
        tb <- fread(p, data.table = F, header = TRUE) %>% select(symbol, read_count, fpkm, fpkm_uq) %>% mutate(Sample=sample)
        if(nrow(tb) < 5) {
            print(sample)
        }
        return(tb)
    }) %>% rbindlist()


count.table <- all.counts %>% dcast(symbol~Sample, value.var = 'read_count', fill = 0)

count.table %>%
    fwrite('Current_readcount_table.tsv',
           sep='\t', row.names = F, col.names = T)

fpkm.table <- all.counts %>% dcast(symbol~Sample, value.var = 'fpkm', fill = 0)
fpkm.table %>%
  fwrite('Current_FPKM_table.tsv',
         sep='\t', row.names = F, col.names = T)

fpkm.uq.table <- all.counts %>% dcast(symbol~Sample, value.var = 'fpkm_uq', fill = 0)
fpkm.uq.table %>%
  fwrite('Current_FPKM_UQ_table.tsv',
         sep='\t', row.names = F, col.names = T)
