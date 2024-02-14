#!/usr/bin/env python
# coding: utf-8

# In[737]:


### RCC_xenium_nolan_neighborhood1_auto
###
###
### - optimized for single sample
### desc: based on Neighborhood Identification jupyter notebook fron nolan lab https://github.com/nolanlab/NeighborhoodCoordination.git
### based on RCC_xenium_nolan_neighborhood1


import argparse
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import time
import sys
import matplotlib.pyplot as plt
from sklearn.cluster import MiniBatchKMeans
import seaborn as sns
import os
import random
import warnings

## Redirect warnings to stdout
warnings.showwarning = lambda msg, *args, **kwargs: print(msg, file=sys.stdout)


## for parallel processing
os.environ["LOKY_MAX_CPU_COUNT"] = "4"

parser = argparse.ArgumentParser(description='based on Neighborhood Identification jupyter notebook fron nolan lab https://github.com/nolanlab/NeighborhoodCoordination.git')
parser.add_argument('-s','--section_id', help='Xenium section_id', required=True)
parser.add_argument('-K','--k_nearest_nbors', help='How many nbors in the nhood', required=True)
parser.add_argument('-n','--nhoods', help='How many nhoods', required=True)
parser.add_argument('-o','--output_dir', help='Output directory', required=True)
parser.add_argument('-i','--input_table_path', help='Table with cell coordinates and cell types', required=True)
parser.add_argument('--analysis_id', help='Analysis subfolder name', required=True)


#
args = vars(parser.parse_args())

print(args)

## assign arguments
section_id = args["section_id"]
output_dir = args["output_dir"]
path_to_data = args["input_table_path"]
analysis_type = args["analysis_id"]
K = [int(args["k_nearest_nbors"])]
n_nhoods = int(args["nhoods"])



# ### Fill In:
#     -K = individual window sizes to collect.  
#         -This accelerates the process of trying different window sizes by collecting them all in one step.
#     -path_to_data:  file path for csv file with x,y cluster information
#     -X:  column to use for X location
#     -Y:  column to use for Y location
#     -reg:  column to use for region (should be different for every tissue)
#     -file_type('csv','pickle'):  file type
#     -cluster_col : column name of cell type cluster to be used for neighborhoods
#     -cellhierpath:  path to cellhier libary (should only have to set once per computer
#     -keep_cols:  columns to transfer from original dataframe to dataframe with windows and neighborhoods


X = 'x_centroid'
Y = 'y_centroid'
reg = 'reg'
file_type = 'csv'
cluster_col = 'celltype_final'
keep_cols = [X,Y,reg,cluster_col]
save_path = '.'


def main():

    # create output directory
    try:
        os.makedirs(output_dir)
    except OSError as error:
        pass
    # set output directory
    os.chdir(output_dir)


    ## create analysis sub-directory
    output_subdir = "./" + analysis_type
    try:
        os.makedirs(output_subdir)
    except OSError as error:
        pass
    # set output directory
    os.chdir(output_subdir)



    ## read in data and do some quick data rearrangement
    n_neighbors = max(K)
    assert (file_type=='csv' or file_type =='pickle') #


    if file_type == 'pickle':
        cells_raw = pd.read_pickle(path_to_data)
    if file_type == 'csv':
        cells_raw = pd.read_csv(path_to_data, sep = "\t")

    cells_raw = pd.concat([cells_raw,pd.get_dummies(cells_raw[cluster_col])],axis = 1)
    cells_raw = cells_raw.assign(reg = "Whole")
    #cells_raw = cells_raw.reset_index() #Uncomment this line if you do any subsetting of dataframe such as removing dirt etc or will throw error at end of next next code block (cell 6)
    sum_cols = cells_raw[cluster_col].unique()
    values = cells_raw[sum_cols].values

    ## find windows for each cell in each tissue region
    tissue_group = cells_raw[[X,Y,reg]].groupby(reg)
    exps = list(cells_raw[reg].unique())
    tissue_chunks = [(time.time(),exps.index(t),t,a) for t,indices in tissue_group.groups.items() for a in np.array_split(indices,1)] 
    tissues = [get_windows(job,n_neighbors,exps, tissue_group, X = X, Y = Y) for job in tissue_chunks]


    ## for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
    out_dict = {}
    for k in K:
        for neighbors,job in zip(tissues,tissue_chunks):

            chunk = np.arange(len(neighbors))#indices
            tissue_name = job[2]
            indices = job[3]
            window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
            out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)

    ## concatenate the summed windows and combine into one dataframe for each window size tested.
    windows = {}
    for k in K:
       
        window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],axis = 0)
        window = window.loc[cells_raw.index.values]
        window = pd.concat([cells_raw[keep_cols],window],axis = 1)
        windows[k] = window


    ## This is not important
    cell_order = list(cells_raw[cluster_col].unique())
    cell_order


    ## calculate neighborhoods 
    cells, fc = nhoods_calculate(k=K[0], n_neighborhoods = n_nhoods,cells = cells_raw, windows = windows, sum_cols = sum_cols, values = values)

    ## make plots
    nhoods_plot(cells = cells, fc =fc, k=K[0], n_neighborhoods = n_nhoods)
    

def get_windows(job,n_neighbors,exps, tissue_group, X, Y):
    '''
    For each region and each individual cell in dataset, return the indices of the nearest neighbors.
    
    'job:  meta data containing the start time,index of region, region name, indices of region in original dataframe
    n_neighbors:  the number of neighbors to find for each cell
    '''
    start_time,idx,tissue_name,indices = job
    job_start = time.time()
    
    print ("Starting:", str(idx+1)+'/'+str(len(exps)),': ' + exps[idx])

    tissue = tissue_group.get_group(tissue_name)
    to_fit = tissue.loc[indices][[X,Y]].values

#     fit = NearestNeighbors(n_neighbors=n_neighbors+1).fit(tissue[[X,Y]].values)
    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[[X,Y]].values)
    m = fit.kneighbors(to_fit)
#     m = m[0][:,1:], m[1][:,1:]
    m = m[0], m[1]
    

    #sort_neighbors
    args = m[0].argsort(axis = 1)
    add = np.arange(m[1].shape[0])*m[1].shape[1]
    sorted_indices = m[1].flatten()[args+add[:,None]]

    neighbors = tissue.index.values[sorted_indices]
   
    end_time = time.time()
   
    print ("Finishing:", str(idx+1)+"/"+str(len(exps)),": "+ exps[idx],end_time-job_start,end_time-start_time)
    return neighbors.astype(np.int32)


def nhoods_calculate(k = 10, n_neighborhoods = 10, enrich_plot = True, spatial_plot = True,
                    cells = None, windows = None, sum_cols = None, values = None):
    
    ## set up
    neighborhood_name = "neighbors" + str(k) + "_hoods" + str(n_neighborhoods)
    k_centroids = {}

    ### calculate neighborhoods
    windows2 = windows[k]

    km = MiniBatchKMeans(n_clusters = n_neighborhoods,random_state=0)

    labelskm = km.fit_predict(windows2[sum_cols].values)
    k_centroids[k] = km.cluster_centers_
    cells[neighborhood_name] = labelskm
    cells[neighborhood_name] = cells[neighborhood_name].astype('category')
    ###
    
    ### this plot shows the types of cells (ClusterIDs) in the different niches (0-7)
    nhood_clusters = (k_centroids[k])
    tissue_avgs = values.mean(axis = 0)
    fc = np.log2(((nhood_clusters+tissue_avgs)/(nhood_clusters+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
    fc = pd.DataFrame(fc,columns = sum_cols)
    # s=sns.clustermap(fc.loc[[0,2,3,4,5,6,7,8,9],cell_order], vmin =-3,vmax = 3,cmap = 'bwr',row_cluster = False)
    # s.savefig("raw_figs/celltypes_perniche_10.pdf")
    
    return cells, fc


    ###




def nhoods_plot(cells, fc, k, n_neighborhoods):
    
    ##
    neighborhood_name = "neighbors" + str(k) + "_hoods" + str(n_neighborhoods)
    
    
    ### plot celltype enrichment matrix
    s2=sns.clustermap(fc.transpose(), vmin =-3,vmax = 3,cmap = 'bwr',
                      row_cluster = True, col_cluster = True,
                      yticklabels=True,
                     figsize = (fc.shape[0]*0.5+4,fc.shape[1]*0.3+3))
    s2.savefig("celltypes_per_" + neighborhood_name + ".pdf")
    ###

    
    
    ### Make nhood spatial plot
    
    ## create random palette
    nhood_levels = cells[neighborhood_name].unique()
    colors = random.sample(list(sns.color_palette("colorblind",n_colors = n_neighborhoods).as_hex()), k = n_neighborhoods)
    pal = dict(zip(nhood_levels, colors))

    ## plot
    with plt.rc_context({'axes.grid': False}):
        with plt.style.context('dark_background'):

            # Create a seaborn scatter plot
            p = sns.lmplot(data = cells,x = X,y=Y,hue = neighborhood_name,
                           palette = pal,
                        height = 10, aspect = 1, fit_reg = False,
                           scatter_kws={"s": 2, "marker" : "circle", "edgecolors" : "none"}
                          )
            # change marker size
            for lh in p._legend.legend_handles: 
                lh._sizes = [150]
            # change legend title
            p.legend.set_title("")
            # change legend position
            sns.move_legend(p, "upper left", bbox_to_anchor=(0.75, 1))

            # Save the figure with specific dimensions
            plt.savefig("spatial" + neighborhood_name + ".png", dpi=150, bbox_inches='tight')

        
    ## v2 new colors. create random palette
    nhood_levels = cells[neighborhood_name].unique()
    colors = random.sample(list(sns.color_palette("colorblind",n_colors = n_neighborhoods).as_hex()), k = n_neighborhoods)
    pal = dict(zip(nhood_levels, colors))

    ## plot
    with plt.rc_context({'axes.grid': False}):
        with plt.style.context('dark_background'):

            # Create a seaborn scatter plot
            p = sns.lmplot(data = cells,x = X,y=Y,hue = neighborhood_name,
                           palette = pal,
                        height = 10, aspect = 1, fit_reg = False,
                           scatter_kws={"s": 2, "marker" : "circle", "edgecolors" : "none"}
                          )
            # change marker size
            for lh in p._legend.legend_handles: 
                lh._sizes = [150]
            # change legend title
            p.legend.set_title("")
            # change legend position
            sns.move_legend(p, "upper left", bbox_to_anchor=(0.75, 1))


            # Save the figure with specific dimensions
            plt.savefig("spatial" + neighborhood_name + "_alt_colors.png", dpi=150, bbox_inches='tight')

    

    ###


if __name__ == "__main__":
    # pass
    # print(output_dir)
    main()

