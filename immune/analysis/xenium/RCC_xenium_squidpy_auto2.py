#!/usr/bin/env python
# coding: utf-8


### RCC_xenium_squidpy_auto1
###
###
###
###
###
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
###



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


# if (a < b)
#     print('a is less than b')


section_id = args["section_id"]
out = args["output_folder"]
input_data_folder = args["input_folder"]
# In[642]:


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
addmeta.dtypes


# In[644]:




sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor="white")


coordinates_mine_df = pd.read_csv(input_data_folder + "/cells.csv.gz")
coordinates_mine_df.rename(columns={"cell_id": "barcode"}, inplace = True)
coordinates_mine_df.head()


# In[649]:


coordinates_mine_filtered_df= coordinates_mine_df[coordinates_mine_df["barcode"].isin(addmeta["barcode"])]
print(coordinates_mine_df.shape)
print(addmeta.shape)
print(coordinates_mine_filtered_df.shape)


# In[650]:


## make sure coordinates and celltype dfs have the same order
coordinates_mine_filtered_df.sort_values("barcode")
addmeta.sort_values("barcode")


# In[651]:


coordinates_mine_array = coordinates_mine_filtered_df.loc[:,["x_centroid","y_centroid"]].to_numpy()
np.shape(coordinates_mine_array)


# In[652]:


## generate fake counts for 4 fake genes
counts_mine = np.reshape((0,) * coordinates_mine_array.shape[0]*4, (coordinates_mine_array.shape[0],4))
counts_mine.shape


# In[653]:


adata_mine = AnnData(counts_mine, obsm={"spatial": coordinates_mine_array})


# In[654]:


adata_mine


# In[655]:


## add celltypes
addmeta = addmeta.set_axis(adata_mine.obs.index, axis = "index")
# print(addmeta.head())
adata_mine.obs = pd.concat([adata_mine.obs, addmeta], axis = 1)
# adata_mine = AnnData(adata_mine.X, obs = pd.concat([adata_mine.obs, addmeta], axis=1),var = adata_mine.var)
adata_mine


# In[656]:


plt.rcParams.update({
    "figure.facecolor":  (1.0, 0.0, 0.0, 0),  # red   with alpha = 30%
#     "axes.facecolor":    (0.0, 1.0, 0.0, 0.5),  # green with alpha = 50%
    "savefig.facecolor": (0.0, 0.0, 1.0, 0),  # blue  with alpha = 20%
})

sq.pl.spatial_scatter(
    adata_mine,
    library_id="spatial",
    shape=None,
    color=[
        "celltype_final",
    ],
    figsize = (5,5),
    size = 1,
    wspace=0.4,
)


# In[329]:


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
    plt.savefig("spatial_plot.pdf")
plt.show()


# In[333]:


## This dataset contains Leiden cluster groups’ annotations in anndata.AnnData.obs, 
## which are used for calculation of centrality scores.
## First, we need to compute a connectivity matrix from spatial coordinates to calculate the centrality scores.
## We can use squidpy.gr.spatial_neighbors for this purpose. We use the coord_type="generic" based on the data 
## and the neighbors are classified with Delaunay triangulation by specifying delaunay=True.

sq.gr.spatial_neighbors(adata_mine, coord_type="generic", delaunay=True)


# In[334]:


## Compute centrality scores
sq.gr.centrality_scores(adata_mine, cluster_key="celltype_final")


# In[336]:


## vizualize
sq.pl.centrality_scores(adata_mine, cluster_key="celltype_final", figsize=(20, 5))
plt.savefig("centrality_plot.pdf")


# In[338]:

## subset for speed inmporevement
adata_subsample = sc.pp.subsample(adata_mine, fraction=0.5, copy=True)


# In[341]:


## Calculate co-occurence probability
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="celltype_final",
    n_jobs = 30
)



## Plot co-occurence probability

for Type in ["Macrophages","T or NK","ccRCC cancer cell","B or Plasma cells", "Endothelial"]:
    
    try:
        sq.pl.co_occurrence(
            adata_subsample,
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
                adata_subsample,
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

# In[358]:


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



