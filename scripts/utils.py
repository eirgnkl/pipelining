import pandas as pd # type: ignore
import yaml # type: ignore
import hashlib
import os
import numpy as np
from scipy.spatial import cKDTree


# from https://github.com/HCA-integration/hca_integration_toolbox/blob/main/workflow/utils/misc.py#L129
def create_hash(string: str, digest_size: int = 5):
    string = string.encode('utf-8')
    return hashlib.blake2b(string, digest_size=digest_size).hexdigest()

def create_tasks_df(config, save=None):
    tasks_df = []
    with open(config, "r") as stream:
        params = yaml.safe_load(stream)
    for task in params['TASKS']: 
        task_dict = params['TASKS'][task]
        params_list = []
        method_dfs = []
        for method in task_dict['methods']:
            method_params = task_dict['methods'][method]
            df_params = pd.read_csv(method_params, sep='\t', index_col=0) #this reads the params.tsv file - \t indic tabsep
            params_list = [str(row) for row in df_params.to_dict(orient='records')]
            method_df = {}
            method_df['params'] = params_list
            method_df['hash'] = [create_hash(row + method + task) for row in params_list]
            method_df['method'] = method
            method_dfs.append(pd.DataFrame(method_df)) #dictionary turns into df
        method_dfs = pd.concat(method_dfs)
        method_dfs['task'] = task
        for key in task_dict:
            if key != 'methods':
                method_dfs[key] = task_dict[key] 
        
        tasks_df.append(method_dfs)
    tasks_df = pd.concat(tasks_df)
    if save is not None:
        tasks_df.to_csv(save, sep='\t')
    return tasks_df


def find_closest_points(adata_rna, adata_msi):
    """
    Find the closest points between two AnnData objects based on their spatial coordinates.

    Parameters:
    adata_rna (AnnData): AnnData object containing RNA data with spatial coordinates in 'spatial_warp'.
    adata_msi (AnnData): AnnData object containing MSI data with spatial coordinates in 'spatial_warp'.

    Returns:
    tuple: A tuple containing:
        - matching_df (DataFrame): DataFrame showing the matching pairs of points and their distances.
        - matching_rna (ndarray): Array of RNA data corresponding to the closest points.
        - matching_msi (ndarray): Array of MSI data.
    """
    # Use the 'spatial_warp' coordinates
    adata_rna_coords = adata_rna.obsm['spatial_warp']
    adata_msi_coords = adata_msi.obsm['spatial_warp']

    # Step 1: Build a spatial tree for `adata_rna`
    tree_adata_rna = cKDTree(adata_rna_coords)

    # Step 2: Query the tree to find the closest point in `adata_rna` for each point in `adata_msi`
    distances, indices = tree_adata_rna.query(adata_msi_coords)

    # Step 3: Create a DataFrame to show the matching pairs
    matching_df = pd.DataFrame({
        'adata_msi_index': np.arange(len(adata_msi_coords)),
        'adata_msi_x': adata_msi_coords[:, 0],
        'adata_msi_y': adata_msi_coords[:, 1],
        'closest_adata_rna_index': indices,
        'adata_rna_x': adata_rna_coords[indices, 0],
        'adata_rna_y': adata_rna_coords[indices, 1],
        'distance': distances
    })

    # Extract the matching RNA and MSI data
    matching_rna = adata_rna.X[indices]
    matching_msi = adata_msi.X

    return matching_df, matching_rna.toarray(), matching_msi
