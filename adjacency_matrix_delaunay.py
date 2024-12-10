# adjacency_matrix_delaunay.py
import numpy as np
import networkx as nx
from libpysal.weights import Delaunay

def compute_delaunay_matrices(adata):
    """
    Perform Delaunay triangulation on spatial coordinates and return the distance matrix and adjacency matrix.

    Parameters:
    adata: AnnData object containing the spatial coordinates in adata.obsm['spatial'].

    Returns:
    dist_matrix: A distance matrix (numpy array).
    connectivity_matrix: An adjacency matrix (scipy sparse matrix).
    """
    # Extract spatial coordinates
    coordinates = adata.obsm['spatial']
    
    # Perform Delaunay triangulation using libpysal
    w = Delaunay(coordinates)
    
    # Convert libpysal weights to networkx graph
    G = w.to_networkx()
    
    # Create a distance matrix and adjacency matrix
    dist_matrix = np.full((len(coordinates), len(coordinates)), np.inf)
    for edge in G.edges():
        dist_matrix[edge[0], edge[1]] = np.linalg.norm(coordinates[edge[0]] - coordinates[edge[1]])
        dist_matrix[edge[1], edge[0]] = dist_matrix[edge[0], edge[1]]  # Symmetric

    # Create a connectivity matrix (adjacency matrix)
    connectivity_matrix = nx.adjacency_matrix(G)
    
    return dist_matrix, connectivity_matrix