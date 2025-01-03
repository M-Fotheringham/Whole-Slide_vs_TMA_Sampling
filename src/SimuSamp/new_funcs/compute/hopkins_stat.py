import numpy as np
from SimuSamp.new_funcs.compute.n_neighbours import neighbours


def hopkins_stat(gdf, k=0.1, prop_k=True):
    """
    Compute the Hopkins statistic for the input GeoDataFrame.
    The Hopkins statistic is a measure of the spatial clustering of a dataset.
    The Hopkins statistic is computed as follows:
    1. Randomly select k points from the dataset.
    2. Compute the distance from each point to n nearest neighbours.
    3. Compute the distance from random points to n nearest neighbours.
    4. Compute the ratio of the sum of the distances.
    Args:
        gdf (geopandas.GeoDataFrame): The input GeoDataFrame.
        k (int): The number of random points to select (min 5).
    Returns:
        float: The Hopkins statistic.
    """
    gdf = gdf.reset_index(drop=True)

    if len(gdf) <= 5:
        return np.nan

    # Compute the nearest neighbour distances for the input GeoDataFrame
    dataset_sum = np.nansum(neighbours(gdf, 1))

    if prop_k:
        k = int(len(gdf) * k)

    if k < 5:
        k = 5

    # Randomly select k points from the dataset
    random_indices = np.random.choice(len(gdf), k, replace=False)
    random_df = gdf.iloc[random_indices]

    # Compute the distance from each random point to its nearest neighbour
    random_sum = np.nansum(neighbours(random_df, 1))

    # Compute the ratio of the sum of the distances
    hopkins_statistic = dataset_sum / (dataset_sum + random_sum)

    return hopkins_statistic
