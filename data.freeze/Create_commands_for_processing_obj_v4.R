#

library(googlesheets4)
library(stringr)
library(dplyr)
library(data.table)
gs4_deauth()
rds.obj <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 3, trim_ws = T)
rds.obj$Sample = paste(rds.obj$`Disease Type`, rds.obj$`Sample ID`, sep = '_')
rds.obj <- rds.obj %>% mutate(Keep = case_when (`Data Type`=='snATAC' & `Cellranger version` == "v2.0" ~ 'TRUE', 
                                                 TRUE ~'FALSE'))
rds.obj <- rds.obj %>% dplyr::filter(Keep == 'TRUE')
rds.obj <- rds.obj %>% dplyr::select(Sample,`Disease Type`, `Data folder`) %>% arrange(Sample)

cellranger <- read_sheet("https://docs.google.com/spreadsheets/d/1lfPnSIweII4cUC5wWVfBIFjKNdwWUI_CUThE2M7NzOs/edit?usp=sharing", sheet = 2, trim_ws = T)
cellranger$Sample = paste(cellranger$`Disease Type`, cellranger$`Sample ID`, sep = '_')
cellranger <- cellranger %>% mutate(Keep = case_when (`Data Type`=='snATAC' & `Cellranger version` == "v2.0" ~ 'TRUE', 
                                                TRUE ~'FALSE'))
cellranger <- cellranger %>% dplyr::filter(Keep == 'TRUE')
cellranger <- cellranger %>% dplyr::select(Sample,`Disease Type`, `Data Folder`) %>% arrange(Sample)

samples.id <- rds.obj$Sample %>% as.character()
com.list <- NULL
for (i in 1:length(samples.id)) {
  s=samples.id[i]
  cr.path <- cellranger %>% filter (Sample==s) %>% pull (`Data Folder`)
  tmux.session <- paste0(tolower(cellranger$`Disease Type`[i]), i)
  command <- str_glue("tmux new -d -s {tmux.session} \ntmux send-keys -t {tmux.session} 'conda activate signac && Rscript /diskmnt/Projects/Users/allakarpova/scripts/snATAC_my/atac_pipeline_v4.0.R \\\
                      -s {s} \\\
                      -a {cr.path} \\\
                      -o /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v4.0 \\\
                      --prf_min 1000 --pct_min 15 --ns_max 5 --pc_first 2 \\\
                      --chrom_size /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt \\\
                      --macs2_path /diskmnt/Projects/Users/allakarpova/Tools/anaconda3/envs/signac/bin/macs2' ENTER\n")
  
  com.list <- c(com.list, command)
}

fwrite(data.frame(com.list), '~/R_working_dir/scripts/snATAC/data.freeze/Create_commands_for_processing_obj_v4.sh', sep = '\t', quote = F, col.names = F)

