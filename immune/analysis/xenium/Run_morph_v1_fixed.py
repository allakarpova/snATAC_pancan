import Morph
import numpy as np 
from skimage.draw import polygon2mask
import csv
import skimage
import matplotlib.pyplot as plt
import argparse
import os
import Morph.features

def crop(image, coord, S, d):
    mask = polygon2mask(image.shape, coord // d)
    imagecrop = np.multiply(image, mask)
    return skimage.morphology.closing(imagecrop, S)


def main():

    parser = argparse.ArgumentParser(description="Run Morph pipeline on Xenium data." )
    parser.add_argument("--input", required=True, help="Input folder containing transcripts.csv(.gz) and cells.csv(.gz)")
    parser.add_argument("--output", required=True, help="Output folder to write CSVs (will be created if missing)" )

    # Parameters from the #parameters cell in the notebook
    parser.add_argument("--d", type=int, default=10, help="Tile size / neighborhood parameter (default: 10)")
    parser.add_argument("--genes", nargs="*", default=["KRT14","KRT5","GPRC5A","PI3","SERPINB2"],
                        help="Gene set G used by the pipeline (default from notebook)")
    parser.add_argument("--structuring_size", type=int, default=3,
                        help="Size for structuring element S as np.ones((n,n)) (default: 3)" )
    parser.add_argument("--tau", type=int, default=10, help="Number of gene transcripts per tile threshold (default: 10)" )
    parser.add_argument("--lambda_", dest="lambda_", type=int, default=40,
                        help="Minimum region size for area_opening (default: 40)" )

    # Optional polygon file
    parser.add_argument("--polygon", default=None,
                        help="Optional CSV of polygon coordinates to restrict analysis region (e.g., *_tumor_coordinates.csv)." )

    args = parser.parse_args()

    in_dir = os.path.abspath(args.input)
    out_dir = os.path.abspath(args.output)
    os.makedirs(out_dir, exist_ok=True)

    transcripts_path = os.path.join(in_dir, "transcripts.csv.gz")
    print(in_dir)
    print(transcripts_path)
    cells_path = os.path.join(in_dir, "cells.csv.gz")
    polygon_path = args.polygon
    S = np.ones((args.structuring_size, args.structuring_size), dtype=np.uint8)
    G = set(args.genes)

    data = Morph.readers.transcripts(transcripts_path)

    if polygon_path:
        with open(polygon_path, 'r') as f:
            i=0
            c=[]
            reader = csv.reader(f)
            for line in reader:
                i += 1
                if i<4:
                    continue
                c.append((int(float(line[0])), int(float(line[1]))))
        c = np.array(c)
        image = Morph.backbone(data, ['xenium', args.d], ['total', G], ['maximum'], ['custom', crop, c, S, args.d],['binary', args.tau], ['area_opening', args.lambda_], ['blob', S])

    else:
        image = Morph.backbone(data, ['xenium', args.d], ['total', G], ['maximum'], ['closing', S], ['binary', args.tau], ['area_opening', args.lambda_], ['blob', S])


    from matplotlib import cm
    plt.imshow(image.T, cmap=cm.magma_r)
    plt.savefig(os.path.join(out_dir,"output.png"), dpi=300, bbox_inches="tight")
    plt.close()

    cells = Morph.readers.cells(cells_path)
    mapper = Morph.modules.Mapper()
    cells = mapper.xenium(cells, args.d)
    cells_out = os.path.join(out_dir, 'output.csv')
    Morph.writers.xenium(cells_out, image, cells)

    layer = Morph.features.Layer()
    layers = layer.minimum(image) # or layer.maximum(image)
    layersmax = layer.maximum(image)

    layers_min_out = os.path.join(out_dir, 'output_layer_min.csv')
    layers_max_out = os.path.join(out_dir, 'output_layer_max.csv')

    Morph.writers.xenium(layers_min_out, layers, cells)
    Morph.writers.xenium(layers_max_out, layersmax, cells)

    distance = Morph.features.Distance()
    distances = distance.minimum(image) # or distance.maximum(image)
    distancesmax = distance.maximum(image)

    dist_min_out = os.path.join(out_dir, 'output_distances_min.csv')
    dist_max_out = os.path.join(out_dir, 'output_distances_max.csv')

    Morph.writers.xenium(dist_min_out, distances, cells)
    Morph.writers.xenium(dist_max_out, distancesmax, cells)
    print("[Done] All outputs written to:", out_dir)


if __name__ == "__main__":
    main()




