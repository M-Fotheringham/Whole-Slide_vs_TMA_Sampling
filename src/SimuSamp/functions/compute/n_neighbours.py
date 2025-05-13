import numpy as np
from sklearn.neighbors import NearestNeighbors


def neighbours(gdf, k, microns_per_pixel=0.22715):
    """
    Compute the nearest neighbour distances for a GeoDataFrame.
    Args:
        gdf (geopandas.GeoDataFrame): The input GeoDataFrame.
        k (int): The number of nearest neighbours to consider.
    """
    if len(gdf) < k + 1:
        return np.nan
    
    # Get the coordinates of the input GeoDataFrame
    xy = np.array(gdf.geometry.apply(lambda geom: (geom.x, geom.y)).tolist())

    # Compute the distance from each point to n nearest neighbours
    nn = NearestNeighbors(n_neighbors=k+1).fit(xy)
    distances, _ = nn.kneighbors(xy)

    dataset_distances = distances[:, 1:]

    dataset_distances = dataset_distances * microns_per_pixel

    return dataset_distances
