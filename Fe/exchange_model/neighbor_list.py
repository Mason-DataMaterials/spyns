import numpy as np
import pandas as pd

from pymatgen.transformations.standard_transformations import SupercellTransformation

from neighbormodels.structure import from_file
from neighbormodels.neighbors import count_neighbors
from neighbormodels.interactions import build_model



def neighbor_list(
    cell_structure: Structure,
    N: int,
    r: float
)
    """
    Get data frame of pairwise neighbor distances for each atom in the unit cell,
    out to a distance ``r``.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param N: number of sites on one side of box for NxNxN total sites count
    :param r: Radius of sphere.
    :return: number of sites, number of n_th neighbors for each site depending on r and 
             numpy array of neighbors for each site sorted by index and distance.
             
    """
     


    structure_supercell = cell_structure.copy()
    structure_supercell.make_supercell([N,N,N])
    neighbor_data = count_neighbors(cell_structure=structure_supercell, r)
    a = neighbor_data.data_frame \
    .merge(neighbor_data.data_frame.rename(columns={"j": "i", "i": "j"}), how="outer") \
    .sort_values(["i","distance_bin", "j"]) \
    .loc[:, ["i", "j", "distance_bin"]] \
    .reset_index(drop=True)
    
    
    neighbors_per_site = a.groupby(["i"]).count().reset_index().rename(columns={"j": "count"}).loc[:, "count"].values
    sites_per_distance_group = a.groupby(["i", "distance_bin"]).count().reset_index() \
                            .rename(columns={"j": "count"}).loc[:, "count"].values
    neighbor_indices = a["j"].values
    
    #neighbor_count, n_indices = np.unique(sites_per_distance_group, return_index=True)
    neighbor_count = np.dstack ( np.unique(sites_per_distance_group, return_index=True) )    
    neighbor_count.dtype = np.dtype([('v', neighbor_count.dtype), ('i', neighbor_count.dtype)])
    neighbor_count.sort(order='i', axis=1)
    neighbor_count = neighbor_count.flatten()['v'].tolist()
    
    
    return N*N*N, neighbor_count, neighbor_indices
    
