import pandas as pd
import networkx as nx
import sys

def read_file(file_path):
    """
    Reads the tab-separated file and returns a list of edges
    using only the first two columns.
    """
    data = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1])
    edges = list(data.itertuples(index=False, name=None))
    return edges

def find_connected_components(edges):
    """
    Finds and returns the connected components given a list of edges.
    """
    G = nx.Graph()
    G.add_edges_from(edges)
    connected_components = list(nx.connected_components(G))
    return connected_components

def main(file_path):
    edges = read_file(file_path)
    components = find_connected_components(edges)
    
    for i, component in enumerate(components, start=1):
        print(f"Component {i}: {component}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_file>")
    else:
        file_path = sys.argv[1]
        main(file_path)

