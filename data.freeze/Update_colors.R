#make colors for panATAC project
cell.type <- c('Doublet' = 'grey', 
                "Tumor" = '#e7298a', 
                "Fibroblasts" =  '#d95f02',
                "Pericytes" = '#ff7f00',
                "Endothelial" = '#66a61e',
                "Normal epithelial cells"= '#1b9e77', 
                "Endothelial" = "#e6ab02",
                "Macrophages" = '#7570b3', 
                "Mast" = '#6a3d9a',
                "T-cells" = '#1f78b4',
                "NK" = '#a6cee3',
                "Tregs" = '#b3cde3',
                "Plasma" = '#fb9a99',
                "B-cells" = '#e31a1c',
                "B-cells/Plasma" = '#cb181d',
                "DC" = '#cab2d6',
                "Oligodendrocytes" = '#33a02c',
                "Neurons" = '#b2df8a',
                "OPC" = '#fccde5',
                "Erythrocytes" = '#fbb4ae', 
                "Hepatocytes" = '#a65628',
                "Epithelial"= '#1b9e77',
                "Alveolar" = '#b3cde3',
                "Goblet" = '#fed9a6',
                "Acinar" = '#e78ac3',
                "Islets" = '#7fc97f')

cancer.col <- c("BRCA"= "#c70039", 
                "CESC"="#ff5733", 
                "CRC"="#ff8d1a", 
                "GBM"="#ffc300",
                "HNSCC"="#eddd53", 
                "MM"="#add45c", 
                "OV"="#57c785", 
                "PDAC"="#00baad",
                "UCEC"="#2a7b9b", 
                "ccRCC"="#3d3d6b")

colors <- list(cell_type = cell.type,
               Cancer = cancer.col)
saveRDS(colors, '~/R_working_dir/scripts/snATAC/Colors_panatac_v1.0.rds')

