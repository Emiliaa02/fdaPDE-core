import pandas as pd
import json
import os
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import argparse


# Function to plot the Voronoi
def plot_voronoi(shape, MESHES_PATH, VORONOI_PATH):
    points_path = os.path.join(MESHES_PATH, shape, "points.csv")
    edges_path = os.path.join(MESHES_PATH, shape, "edges.csv")
    df_points = pd.read_csv(f"test/data/mesh/{shape}/points.csv").rename(columns={'Unnamed: 0': 'id'})

    _ = plt.figure(figsize=(10, 8))
    _ = plt.scatter(df_points['V1'], df_points['V2'], color='black', s=1)
    id2TriaNodes = {}
    for _, row in df_points.iterrows():
        id2TriaNodes[row['id']] = (row['V1'], row['V2'])

    df_edges = pd.read_csv(f"test/data/mesh/{shape}/edges.csv")
    if 'Unnamed: 0' in df_edges.columns:
        df_edges = df_edges.drop(columns='Unnamed: 0')
    triaEdges = []
    for _, row in df_edges.iterrows():
        triaEdges.append([id2TriaNodes[row['V1']], id2TriaNodes[row['V2']]])

    lc = LineCollection(triaEdges, colors='black', linewidths=0.5)

    # Add to current axes
    ax = plt.gca()
    ax.add_collection(lc)

    voronoi_path = os.path.join(VORONOI_PATH, f"{shape}.json")
    with open(f"test/plots/data/{shape}.json", 'r') as f:
        data = json.load(f)

    nodes = data['nodes']

    xs = [item['coords'][0] for item in nodes]
    ys = [item['coords'][1] for item in nodes]

    _ = plt.scatter(xs, ys, s=5, color='red')

    id2Node = {d['id']: (d['coords'][0], d['coords'][1]) for d in nodes}
    data
    edges = [[id2Node[d['from']], id2Node[d['to']]] for d in data['edges']]
    edges

    lc = LineCollection(edges, colors='red', linewidths=1)

    # Add to current axes
    ax = plt.gca()
    ax.add_collection(lc)  
    
    plt.savefig(f"test/plots/voronoi_{shape}", dpi=300)


def main():
    MESHES_PATH = "test/data/mesh"
    VORONOI_PATH = "test/plots/data"
    parser = argparse.ArgumentParser()
    parser.add_argument("--shape", required=True)
    args = parser.parse_args()

    shape = args.shape

    plot_voronoi(shape, MESHES_PATH, VORONOI_PATH)
 
 
if __name__ == "__main__":
    main()