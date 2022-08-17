library(Seurat)
library(data.table)
library(dplyr)
library(scatterpie)
library(paletteer)
library(ggplot2)

##### Spatial data #######
input.path <- '/Users/allakarpova/lab_Ding/work/Spatial/Interactions/Pilot/merge_st/CRC_112C1_spatial_merged.rds'
st_merged<- readRDS(input.path)
rctd.output2 <- '/Users/allakarpova/lab_Ding/work/Spatial/Interactions/Pilot/RCTD_CRC_112C1_2_close_to_real/RCTDmulti_deconvolution_weights.txt'
rctd.output1 <- '/Users/allakarpova/lab_Ding/work/Spatial/Interactions/Pilot/RCTD_CRC_112C1_1_close_to_real/RCTDmulti_deconvolution_weights.txt'

wd <- '/Users/allakarpova/lab_Ding/work/Spatial/Interactions/Pilot/tumor_subclones_merged_st'
dir.create(wd)
setwd(wd)

colnames(st_merged)
rctd.weigts1 <- fread(rctd.output1, data.table = F) %>% mutate(Barcode = paste('CRC_112C1_1',Barcode, sep = '_')) %>%  data.frame(row.names = 1, check.rows = F, check.names = F)
rctd.weigts2 <- fread(rctd.output2, data.table = F) %>% mutate(Barcode = paste('CRC_112C1_2',Barcode, sep = '_')) %>%  data.frame(row.names = 1, check.rows = F, check.names = F)
# add missing columns of zeros to each matrix
rctd.weigts1 <- cbind(rctd.weigts1, 
                      setNames(data.frame(matrix(ncol = length(setdiff(colnames(rctd.weigts2), colnames(rctd.weigts1))), nrow = nrow(rctd.weigts1), data = 0)),
                              setdiff(colnames(rctd.weigts2), colnames(rctd.weigts1))))
rctd.weigts2 <- cbind(rctd.weigts2, 
                      setNames(data.frame(matrix(ncol = length(setdiff(colnames(rctd.weigts1), colnames(rctd.weigts2))), nrow = nrow(rctd.weigts2), data = 0)),
                               setdiff(colnames(rctd.weigts1), colnames(rctd.weigts2))))

rctd.weigts <- rbind(rctd.weigts1, rctd.weigts2[colnames(rctd.weigts1)])

st_merged <- AddMetaData(st_merged, rctd.weigts)

st_merged@images$CRC_112C1_2@coordinates
toplot <- rctd.weigts %>%
  tibble::rownames_to_column(var = 'Barcode') %>% 
  left_join(rbind(st_merged@images$CRC_112C1_1@coordinates[c('imagerow', 'imagecol')],
                  st_merged@images$CRC_112C1_2@coordinates[c('imagerow', 'imagecol')]) %>%
    tibble::rownames_to_column(var = 'Barcode') %>%
    as_tibble(),
  by = 'Barcode'
) %>%
  left_join(st_merged$orig.ident %>%
            data.frame(check.rows = F, check.names = F) %>%
            setNames(nm = 'orig.ident') %>%
            tibble::rownames_to_column(var = 'Barcode') %>%
            as_tibble(),
            by = 'Barcode'
  )
str(toplot)
p = ggplot() +
  geom_scatterpie(
    aes(y = imagerow, x = imagecol, group = Barcode),
    cols = colnames(rctd.weigts),
    data = toplot,
    pie_scale = 0.3, size = 1,
    color = NA
  ) + 
  coord_fixed() + 
  facet_grid(.~orig.ident) +
  cowplot::theme_cowplot() + scale_y_reverse() +
  scale_fill_paletteer_d("ggthemes::Tableau_20")

ggsave('Spatial_piechart_RCTDmulti_several_cell_types_facet_wrap_by_slice.pdf',plot = p, width =15, height =7, useDingbats = F)



samples <- c('CRC_112C1_1', 'CRC_112C1_2')
coords <- lapply(samples, function (x) { GetTissueCoordinates(
  st_merged[[x]],
  scale = "lowres",
  cols = c( "imagecol","imagerow")) })

st_merged[['RCTDdmulti']] <- CreateAssayObject(counts = t(as.matrix(rctd.weigts)))
DefaultAssay(st_merged) <- 'RCTDdmulti'


SpatialDimPlot(st_merged, label = T) /DimPlot(st_merged)
ggsave('All_cells_dimplot.pdf',useDingbats = F, width = 5, height = 10)


st_merged_tumor <- subset(st_merged, subset = Tumor_c2 >= 0.9)

SpatialDimPlot(st_merged_tumor, label = F, pt.size.factor = 2.5)/DimPlot(st_merged_tumor) & coord_fixed()
ggsave('Pure_tumor_dimplot.pdf',useDingbats = F, width = 5, height = 10)

subpopulations <- FindMarkers(st_merged_tumor, ident.1 = '0', ident.2 = '2', assay = 'SCT', test.use = 'LR')
fwrite(subpopulations, 'Pure_tumor_spots_cluster0_vs2_SCT.txt', sep = '\t', row.names = T)
pos.markers <- subpopulations %>%
  subset(avg_log2FC > 0 & p_val_adj < 0.1)
neg.markers <- subpopulations %>%
  subset(avg_log2FC < 0 & p_val_adj < 0.1)

DefaultAssay(st_merged) <- 'SCT'
p<- SpatialFeaturePlot(object = st_merged, features = rownames(pos.markers)[1:70], ncol = 10)
ggsave('Spatial_positive_pure_tumor_markers_0_vs2_SCT.pdf',plot = p, useDingbats = F, width = 10, height = 70,limitsize = F)

p<- SpatialFeaturePlot(object = st_merged, features = rownames(neg.markers)[1:70], ncol = 10)
ggsave('Spatial_negative_pure_tumor_markers_0_vs2_SCT.pdf',plot = p, useDingbats = F, width = 10, height = 70,limitsize = F)


center_vs_periphery <- FindMarkers(st_merged, ident.1 = '4', ident.2 = '7', assay = 'SCT',test.use = 'LR')
fwrite(subpopulations, 'Pure_tumor_spots_center_vs_periphery.txt', sep = '\t', row.names = T)
pos.markers <- center_vs_periphery %>%
  subset(avg_log2FC > 0 & p_val_adj < 0.1)
neg.markers <- center_vs_periphery %>%
  subset(avg_log2FC < 0 & p_val_adj < 0.1)
p<- SpatialFeaturePlot(object = st_merged, features = rownames(pos.markers)[1:70], ncol = 10)
ggsave('Spatial_positive_markers_center_vs_periphery.pdf',plot = p, useDingbats = F, width = 20, height = 20)

p<- SpatialFeaturePlot(object = st_merged, features = rownames(neg.markers)[1:70], ncol = 10)
ggsave('Spatial_negative_markers_center_vs_periphery.pdf',plot = p, useDingbats = F, width = 20, height = 20)

