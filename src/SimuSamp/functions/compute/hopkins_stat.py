import numpy as np
import pointpats
from shapely.geometry import Point
import geopandas as gpd
from SimuSamp.functions.compute.n_neighbours import neighbours


def hopkins_stat(anno, gdf, k=0.1, prop_k=True):
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

    if prop_k:
        k = int(len(gdf) * k)

    if k < 5:
        k = 5

    # Randomly select k points from the dataset
    df_sample = gdf.sample(k)
    # Compute the distance from each sampled point to its nearest neighbour
    df_sum = np.nansum(neighbours(df_sample, 1))

    # Generate random points within the region
    random_points = pointpats.random.poisson(anno, size=len(gdf))
    random_points = [Point(p) for p in random_points]
    random_gdf = gpd.GeoDataFrame(geometry=random_points)
    # Subset k random points
    random_gdf = random_gdf.sample(k)
    # Compute the distance from each random point to its nearest neighbour
    random_sum = np.nansum(neighbours(random_gdf, 1))

    # Compute the ratio of the sum of the distances
    hopkins_statistic = random_sum / (df_sum + random_sum)

    return hopkins_statistic
