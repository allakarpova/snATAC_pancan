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
    addmeta = pd.read_csv(args["celltype_meta"])
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


    with plt.style.context('dark_background'):
        sq.pl.spatial_scatter( 
            adata_mine,
            library_id="spatial",
            shape=None,
            palette = 'Paired',
            color=[
                "celltype_final",
            ],
            figsize = (12,8),
            size = 1,
            dpi = 100,
    #         save = "spatial_plot.pdf",
            wspace=0.4)
        plt.savefig("spatial_plot.pdf",facecolor="black", transparent=False)
    plt.show()

if __name__ == "__main__":
    main()
