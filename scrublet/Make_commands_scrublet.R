library(googlesheets4)
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape))
suppressMessages(library(data.table))
library(foreach)
library(purrr)
library(glue)

gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 2, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data Folder`)
samples <- samples %>% dplyr::filter(`Data Type` == '10x_SC_Multi_ATAC_SEQ')

samples.id <- samples$Sample %>% as.character()
samples.folder <- samples$`Data Folder` %>% as.character()


commands.run <- foreach (s = samples.id, p = samples.folder, .combine = 'c')  %do% {
  c <- paste0("tmux new -d -s ",s,"  'python3 /diskmnt/Projects/Users/allakarpova/Projects/gbm/scripts/run_scrublet_v2.py -s ", 
             s, " -c ", p, " -u 0.2'")
  return (c)
}

commands.run <- data.frame(commands.run)
fwrite(commands.run, '~/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_commands_panatac_comboRNA.sh', sep = '\t', quote = F, col.names = F)

#################
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 2, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data Folder`)
samples <- samples %>% dplyr::filter(`Data Type` == '10x_SC_Multi_ATAC_SEQ')

samples.id <- samples$Sample %>% as.character()
samples.folder <- samples$`Data Folder` %>% as.character()


commands.run <- foreach (s = samples.id, p = samples.folder, .combine = 'c')  %do% {
  c <- paste0("tmux new -d -s ",s,"  'python3 /diskmnt/Projects/Users/allakarpova/scripts/snATAC_my/run_scrublet_comboRNA.py -s ", 
              s, " -c ", p, " -u 0.2'")
  return (c)
}

commands.run <- data.frame(commands.run)
fwrite(commands.run, '~/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_commands_panatac_comboRNA_v2.sh', sep = '\t', quote = F, col.names = F)

###### comboATAC
#################
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 2, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data Folder`)
samples <- samples %>% dplyr::filter(`Data Type` == '10x_SC_Multi_ATAC_SEQ')

samples.id <- samples$Sample %>% as.character()
samples.folder <- samples$`Data Folder` %>% as.character()


commands.run <- foreach (s = samples.id, p = samples.folder, .combine = 'c')  %do% {
  c <- paste0("tmux new -d -s ",s,"  'python3 /diskmnt/Projects/Users/allakarpova/scripts/snATAC_my/run_scrublet_comboATAC.py -s ", 
              s, " -c ", p, " -u 0.2'")
  return (c)
}

commands.run <- data.frame(commands.run)
fwrite(commands.run, '~/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_commands_panatac_comboATAC.sh', sep = '\t', quote = F, col.names = F)


#now for atac samples
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 2, trim_ws = T)
samples$Keep <- samples$`Include in the downstream analysis` %>% unlist()
samples$Sample = paste(samples$`Disease Type`, samples$`Sample ID`, sep = '_')
samples <- samples %>% dplyr::filter(Keep == 'TRUE')
samples <- samples %>% dplyr::filter(`Cellranger version` == 'v2.0')
samples <- samples %>% dplyr::select(Sample, `Data Type`, `Data Folder`)
samples <- samples %>% dplyr::filter(`Data Type` == 'snATAC')

samples.id <- samples$Sample %>% as.character()
samples.folder <- samples$`Data Folder` %>% as.character()

commands.run <- map2(samples.id, samples.folder, function (s, path) {
  glue( "tmux new -d -s {s} 'source /diskmnt/Projects/Users/allakarpova/Tools/scrublet/bin/activate&&python3 /diskmnt/Projects/Users/allakarpova/scripts/snATAC_my/run_scrublet_regAtac.py -s {s} -c {path} -u 0.2 -r 0.15'")
})

fwrite(commands.run, '~/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_commands_panatac_regATAC_20210831.sh', sep = '\n', quote = F, col.names = F)


commands.run <- foreach (s = samples.id, p = samples.folder, .combine = 'c')  %do% {
  c <- paste0("python3 /diskmnt/Projects/Users/allakarpova/scripts/snATAC_my/run_scrublet_regAtac.py -s ", 
              s, " -c ", p, " -u 0.2")
  return (c)
}
commands.run <- data.frame(commands.run)
fwrite(commands.run, '~/lab_Ding/work/single_cell/snATAC_pancan/scrublet/Scrublet_commands_panatac_regATAC_notmux.sh', sep = '\t', quote = F, col.names = F)

