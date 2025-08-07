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

    # In[642]:


    # create output directory
    # create output directory
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


    

    # ## modify celltypes
    # # Define conditions and corresponding values for the new column
    # conditions = [addmeta['celltype_final'].str.contains("Macrophages"), addmeta['celltype_final'] == 'cancer cell', addmeta['celltype_final'] == 'Endothelial']
    # values = [addmeta['celltype_final'], "cancer cell", "Endothelial"]

    # # Use numpy's select() function to apply the conditions and assign values
    # addmeta['celltype_final'] = np.select(conditions, values, default="Other")
    # # make categorical
    # addmeta["celltype_final"] = pd.Categorical(addmeta["celltype_final"])

    ## remove rare celltypes
    # for Type in list(addmeta["celltype_final"].cat.categories):
    addmeta["celltype_final"] = addmeta["celltype_final"].astype(str)


    ##QC
    # celltype_final_column = addmeta["celltype_final"].to_list()
    # addmeta["celltype_final"] = [el[0:2] for el in celltype_final_column]


    celltype_counts = addmeta["celltype_final"].value_counts()
    # print(celltype_counts)
    # print(celltype_counts.index.values)
    # print(celltype_counts > 1000)
    common_celltypes =  celltype_counts.index.values[celltype_counts > 10]
    
    ## QC
    # common_celltypes =  ["Podocyte"]

    #
    # print(common_celltypes)
    

    # addmeta["celltype_final"] = ['Other' if ct not in common_celltypes else ct
    #                              for ct in addmeta["celltype_final"]]

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


    # with open('adata_unfiltered.pkl', 'wb') as outp:
    #     pickled_list = [adata_mine]
    #     pickle.dump(pickled_list, outp, pickle.HIGHEST_PROTOCOL)

    



    ## filter object by barcodes in celltype table
    adata_mine=adata_mine[addmeta["barcode"].to_list()].copy()
    ## add celltypes
    addmeta.index = addmeta['barcode']
    adata_mine.obs = adata_mine.obs.merge(addmeta, left_index=True, right_index = True)


    # ## add celltypes
    # addmeta = addmeta.set_axis(adata_mine.obs.index, axis = "index")
    # # print(addmeta.head())
    # adata_mine.obs = pd.concat([adata_mine.obs, addmeta], axis = 1)
    # # adata_mine = AnnData(adata_mine.X, obs = pd.concat([adata_mine.obs, addmeta], axis=1),var = adata_mine.var)
    # adata_mine


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


    ## This dataset contains Leiden cluster groups’ annotations in anndata.AnnData.obs, 
    ## which are used for calculation of centrality scores.
    ## First, we need to compute a connectivity matrix from spatial coordinates to calculate the centrality scores.
    ## We can use squidpy.gr.spatial_neighbors for this purpose. We use the coord_type="generic" based on the data 
    ## and the neighbors are classified with Delaunay triangulation by specifying delaunay=True.

    sq.gr.spatial_neighbors(adata_mine, coord_type="generic", delaunay=True)

    # with open('adata_filtered.pkl', 'wb') as outp:
    #     pickled_list = [adata_mine]
    #     pickle.dump(pickled_list, outp, pickle.HIGHEST_PROTOCOL)


    ## subset for speed inmporevement
    adata_mine = sc.pp.subsample(adata_mine, fraction=1, copy=True)



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


    ## Calculate co-occurence probability
    sq.gr.co_occurrence(
        adata_mine,
        cluster_key="celltype_final",
        n_splits = 25,
        interval = np.arange(1, 500, 10),
        # interval = 50,
        n_jobs = 50,
        show_progress_bar=False
    )

    # with open('adata_with_co_occurrence.pkl', 'wb') as outp:
    #     pickled_list = [adata_mine]
    #     pickle.dump(pickled_list, outp, pickle.HIGHEST_PROTOCOL)


    ## save co_occurence calculation results
    celltypes = adata_mine.obs.celltype_final.unique()
    print(sorted(celltypes))
    sorted_celltypes = sorted(celltypes)

    distance_measure_point = list(adata_mine.uns["celltype_final_co_occurrence"]["interval"])

    occur_prob_tup = tuple(el for el in adata_mine.uns["celltype_final_co_occurrence"]["occ"])
    occur_prob_2d = np.row_stack(occur_prob_tup)
    print(occur_prob_2d.shape)
    # np.row_stack()

    ## convert to df
    df = pd.DataFrame(occur_prob_2d)
    df.columns = ["V" + str(el) for el in distance_measure_point][0:-1]
    print(df.columns)
    ## add column with "celltype of interest"
    x=[]
    [x.extend([str(el)]*int(df.shape[0]**0.5)) for el in sorted_celltypes]

    #
    df["celltype"] = x

    ## add column with "neighbor_type"
    y=sorted_celltypes*len(sorted_celltypes)
    ## add 
    df["neighbor_type"] = y

    ## save
    df.to_csv('co_occurence_2d_mat.tsv', sep="\t", index = False, header = True)

    ## celltypes
    celltypes = adata_mine.obs.celltype_final.unique()

    ## Plot co-occurence probability

    # for Type in ["Macrophages","T or NK","ccRCC cancer cell","B or Plasma cells", "Endothelial"]:
    for Type in [el for el in celltypes if el.find("Macrophages") != -1]:
        
        try:
            sq.pl.co_occurrence(
                adata_mine,
                cluster_key="celltype_final",
                clusters=Type,
                figsize=(8, 6),
            )
            #
            plt.savefig(Type + "_cooccur_prob.pdf")
        except ValueError as error:
            print(error)


    ## plot spatial celltype distrbution
    for Type in list(addmeta["celltype_final"].cat.categories):

        try:

            with plt.style.context('dark_background'):
                sq.pl.spatial_scatter(
                    adata_mine,
                    library_id="spatial",
                    shape=None,
                    color="celltype_final",
                    groups = Type,
                    figsize = (8,6),
                    size = 2,
                    dpi = 100,
            #         save = "spatial_plot.pdf",
                    wspace=0.4)
                plt.savefig(Type+"_spatial_plot.pdf")

        except ValueError as error:
            print(error)


    ## run Neighbors enrichment analysis
    sq.gr.nhood_enrichment(adata_mine, cluster_key="celltype_final")


    ## plot Neighbors enrichment analysis
    p=sq.pl.nhood_enrichment(
        adata_mine,
        
        cluster_key="celltype_final",
        figsize=(4, 4),
        annotate = False,
        title="Neighborhood enrichment adata",
    #     save = "enrichment_heatmap2.pdf"
    #     return_ax = True,
    #     ax=ax[0],
    )

    plt.savefig("enrichment_heatmap.pdf", bbox_inches='tight')


if __name__ == "__main__":
    main()

