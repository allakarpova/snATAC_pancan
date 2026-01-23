#!/usr/bin/env python
# coding: utf-8


### Liver_xenium_squidpy_auto4
### v4
### - ignore celltypes < 100 cells
### - co-ocur prob for 500um distance
### - ECs and cc for co-occur prob
### - macrophages as celltypeoi
### - send warnings to stdout
### - do not subsample
### uses scanpy to process and annotate data
### desc: run squidpy co-occurence probability and Neighbors enrichment analyses
### based on RCC_xenium_squidpy1.2



### Import packages
import argparse
import scanpy as sc
import squidpy as sq
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
import warnings
import sys
import shutil
import pickle

###

def centrality_analysis(adata):
    #calculate
    sq.gr.centrality_scores(adata, cluster_key="celltype_final")
    ## vizualize
    sq.pl.centrality_scores(adata, cluster_key="celltype_final", figsize=(20, 5))
    plt.savefig("centrality_plot.pdf")


def main():

    ## Redirect warnings to stdout
    warnings.showwarning = lambda msg, *args, **kwargs: print(msg, file=sys.stdout)


    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-i','--input_folder', help='Input Xenium folder', required=True)
    parser.add_argument('-s','--section_id', help='Xenium section_id', required=True)
    parser.add_argument('-m','--celltype_meta', help='Celltype metadata path', required=True)
    parser.add_argument('-o','--output_folder', help='Output folder', required=True)
    #parser.add_argument('-c','--clinical_df', help='Clinical table path', required=True)
    #parser.add_argument('--celltype_version', help='Clinical table path', default = "default")

    args = vars(parser.parse_args())

    print(args)


    ### Set options
    pd.set_option('display.max_colwidth', None)
    ###


    ## reassign arguments
    section_id = args["section_id"]
    out = args["output_folder"]
    input_data_folder = args["input_folder"]
    output_dir = out + "/" + section_id + "/"

    try:
        os.makedirs(output_dir)
    except OSError as error:
        print("")
    # set output directory
    os.chdir(output_dir)
    #
    os.getcwd()

    # In[643]:


    ### Load stuff
    ## load celltypes
    addmeta = pd.read_csv(args["celltype_meta"], sep='\t')
    print(addmeta.head())
    # change col name
    addmeta.columns = ["barcode","celltype_final"]
    addmeta["celltype_final"] = pd.Categorical(addmeta["celltype_final"])
    print(addmeta.dtypes)


    ## remove rare celltypes
    # for Type in list(addmeta["celltype_final"].cat.categories):
    addmeta["celltype_final"] = addmeta["celltype_final"].astype(str)


    celltype_counts = addmeta["celltype_final"].value_counts()
    # print(celltype_counts)
    # print(celltype_counts.index.values)
    # print(celltype_counts > 1000)
    common_celltypes =  celltype_counts.index.values[celltype_counts > 10]
    
 
    addmeta["celltype_final"][~addmeta["celltype_final"].isin(common_celltypes)] = "Other"
    # make categorical
    addmeta["celltype_final"] = pd.Categorical(addmeta["celltype_final"])

    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor="white")


    coordinates_mine_df = pd.read_csv(input_data_folder + "/cells.csv.gz")
    coordinates_mine_df.rename(columns={"cell_id": "barcode"}, inplace = True)
    coordinates_mine_df.head()
    

    coordinates_mine_filtered_df= coordinates_mine_df[coordinates_mine_df["barcode"].isin(addmeta["barcode"])]
    

    ## make sure coordinates and celltype dfs have the same order
    coordinates_mine_filtered_df.sort_values("barcode")
    addmeta.sort_values("barcode")


    coordinates_mine_array = coordinates_mine_filtered_df.loc[:,["x_centroid","y_centroid"]].to_numpy()
    np.shape(coordinates_mine_array)


    ## generate fake counts for 4 fake genes
    counts_mine = np.reshape((0,) * coordinates_mine_array.shape[0]*4, (coordinates_mine_array.shape[0],4))


    print(f"count matrix shape is {counts_mine.shape}")

    ## create AnnData object
    adata_mine = AnnData(counts_mine, obsm={"spatial": coordinates_mine_array})
    # add cell names
    adata_mine.obs_names = coordinates_mine_filtered_df["barcode"].to_list()



    ## filter object by barcodes in celltype table
    adata_mine=adata_mine[addmeta["barcode"].to_list()].copy()
    ## add celltypes
    addmeta.index = addmeta['barcode']
    adata_mine.obs = adata_mine.obs.merge(addmeta, left_index=True, right_index = True)




    plt.rcParams.update({
        "figure.facecolor":  (1.0, 0.0, 0.0, 0),  # red   with alpha = 30%
    #     "axes.facecolor":    (0.0, 1.0, 0.0, 0.5),  # green with alpha = 50%
        "savefig.facecolor": (0.0, 0.0, 1.0, 0),  # blue  with alpha = 20%
    })

    celltype_colors = {
        "Cancer cells": "#B11226",
        "Cancer cells DCIS": "#D1495B",
        "CCL28 Cancer cells": "#8C1D18",
        "Plasma cancer cells": "#E07A8C",
        "Pre-cancer cells PanIN": "#C44536",
        "Necrosis": "#5A1A1A",
        "Necrotic": "#3D0F0F",

        "T-cells": "#1F77B4",
        "NK cells": "#2C7FB8",
        "B-cells": "#6BAED6",
        "Plasma": "#9ECAE1",
        "Langerhans cells": "#4A90E2",
        "pDC": "#5DA5DA",

        "Macrophages": "#2CA02C",
        "Myeloid": "#31A354",
        "mregDC": "#74C476",
        "DC": "#41AB5D",
        "Granulocytes": "#66C2A4",
        "Neutrophil": "#006D2C",
        "Neutrophils": "#006D2C",
        "Myeloid or Neutrophil": "#99D8C9",
        "GMP": "#A1D99B",
        "Mast": "#78C679",
        "Erythroid": "#FB6A4A",
        "Megakaryocytes": "#CB181D",

        "Endothelial cells": "#17BECF",
        "Lymphatic endothelial cells": "#9EDAE5",
        "LSECs": "#6BAED6",
        "Pericytes": "#8DD3C7",
        "vSMCs": "#4EB3D3",
        "Smooth muscle": "#3690C0",

        "Fibroblasts": "#8C6D31",
        "MSCs": "#B5A642",
        "Adipocytes": "#D9BF77",
        "Skeletal muscle": "#A6761D",
        "Osteoblasts": "#C49C94",
        "Dermal layer": "#B15928",

        "Normal epithelial": "#9467BD",
        "Myoepithelial and Normal ducts": "#7B6FD6",
        "Epithelial collecting duct": "#8E63CE",
        "Cholangiocytes": "#A55194",
        "Goblet cells": "#CE6DBD",
        "Intestinal epithelium": "#DD1C77",
        "Islets": "#E377C2",
        "Normal pancreas": "#F781BF",
        "Keratinocytes": "#FF7F0E",
        "Basal keratinocytes": "#F28E2B",
        "Squamous cells": "#E15759",
        "Hepatocytes": "#BCBD22",

        "Neurons": "#EDC948",
        "Glial": "#F1CE63",
        "Astrocytes": "#EFC94C",
        "Oligodendrocytes": "#FFD92F",
        "OPC": "#FFED6F",

        "Immune": "#BDBDBD",
        "Unknown": "#969696",
        "Low quality": "#737373"
    }
    with plt.style.context('dark_background'):
        sq.pl.spatial_scatter( 
            adata_mine,
            library_id="spatial",
            shape=None,
            palette = celltype_colors,
            color=[
                "celltype_final",
            ],
            figsize = (12,8),
            size = 0.2,
            dpi = 200,
    #         save = "spatial_plot.pdf",
            wspace=0.4)
        plt.savefig("spatial_plot.pdf",facecolor="black", transparent=False)
    plt.show()

if __name__ == "__main__":
    main()
