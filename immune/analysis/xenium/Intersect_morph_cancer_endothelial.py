import Morph
import numpy as np 
from skimage.draw import polygon2mask
import csv
import skimage
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
import argparse
import os
import Morph.features


def compute_layers_distance(image, cells, nameme):
    layer = Morph.features.Layer()
    layers = layer.minimum(image) 
    layersmax = layer.maximum(image)

    layers_min_out = os.path.join(out_dir, nameme+'_layer_min.csv')
    layers_max_out = os.path.join(out_dir, nameme+'_layer_max.csv')

    Morph.writers.xenium(layers_min_out, layers, cells)
    Morph.writers.xenium(layers_max_out, layersmax, cells)

    distance = Morph.features.Distance()
    distances = distance.minimum(image) # or distance.maximum(image)
    distancesmax = distance.maximum(image)

    dist_min_out = os.path.join(out_dir, nameme+'_distances_min.csv')
    dist_max_out = os.path.join(out_dir, nameme+'_distances_max.csv')

    Morph.writers.xenium(dist_min_out, distances, cells)
    Morph.writers.xenium(dist_max_out, distancesmax, cells)

def main():

    parser = argparse.ArgumentParser(description="Run Morph pipeline on Xenium data." )
    parser.add_argument("--input", required=True, help="Input Xenium output folder containing cells.csv(.gz)")
    parser.add_argument("--input_cancer", required=True, help="Input npy saved image otuput of morph with cancer microregions")
    parser.add_argument("--input_endo", required=True, help="Input npy saved image otuput of morph with endothelial microregions")
    
    parser.add_argument("--output", required=True, help="Output folder to write CSVs (will be created if missing)" )


    args = parser.parse_args()
    in_dir = os.path.abspath(args.input)
    in_cancer = os.path.abspath(args.input_cancer)
    in_endo = os.path.abspath(args.input_endo)
    out_dir = os.path.abspath(args.output)
    os.makedirs(out_dir, exist_ok=True)

    image_cancer = np.load(in_cancer)
    image_endo = np.load(in_endo)
    cells_path = os.path.join(in_dir, "cells.csv.gz")

    image_inter = np.multiply(image_cancer>0, image_endo>0).astype(int)
    plt.imshow(image_inter.T, cmap=cm.magma_r)
    plt.savefig(os.path.join(out_dir,"output_endothelial_inside_cancer.png"), dpi=300, bbox_inches="tight")
    plt.close()

    image_out = (image_endo>0).astype(int) - np.multiply(image>0, image_endo>0).astype(int)
    plt.imshow(image_out.T, cmap=cm.magma_r)
    plt.savefig(os.path.join(out_dir,"output_endothelial_outside_cancer.png"), dpi=300, bbox_inches="tight")
    plt.close()

    cells = Morph.readers.cells(cells_path)
    compute_layers_distance(image_inter, cells, 'Endothelial_inside_cancer')
    compute_layers_distance(image_out, cells, 'Endothelial_outside_cancer')

    print("[Done] All outputs written to:", out_dir)


if __name__ == "__main__":
    main()
