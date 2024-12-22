import numpy as np
from sklearn.neighbors import NearestNeighbors


def hopkins_stat(gdf, k):
    """
    Compute the Hopkins statistic for the input GeoDataFrame.
    The Hopkins statistic is a measure of the spatial clustering of a dataset.
    The Hopkins statistic is computed as follows:
    1. Randomly select k points from the dataset.
    2. Compute the distance from each point in the dataset to its nearest neighbour.
    3. Compute the distance from each randomly selected point to its nearest neighbour.
    4. Compute the ratio of the sum of the distances from the dataset points to their nearest neighbours to the sum of the distances from the randomly selected points to their nearest neighbours.
    5. The Hopkins statistic is the ratio of the sum of the distances from the dataset points to their nearest neighbours to the sum of the distances from the randomly selected points to their nearest neighbours.
    Args:
        gdf (geopandas.GeoDataFrame): The input GeoDataFrame.
        k (int): The number of random points to select.
    Returns:
        float: The Hopkins statistic.
    """
    # Get the coordinates of the input GeoDataFrame
    coords = np.array(gdf.geometry.apply(lambda geom: (geom.x, geom.y)).tolist())

    # Compute the distance from each point in the dataset to its nearest neighbour
    nn = NearestNeighbors(n_neighbors=2).fit(coords)
    distances, _ = nn.kneighbors(coords)
    dataset_distances = distances[:, 1]

    # Randomly select k points from the dataset
    random_indices = np.random.choice(len(coords), k, replace=False)
    random_coords = coords[random_indices]

    # Compute the distance from each randomly selected point to its nearest neighbour
    nn = NearestNeighbors(n_neighbors=2).fit(random_coords)
    distances, _ = nn.kneighbors(random_coords)
    random_distances = distances[:, 1]

    # Compute the ratio of the sum of the distances from the dataset points to their nearest neighbours
    # to the sum of the distances from the randomly selected points to their nearest neighbours
    dataset_sum = np.sum(dataset_distances)
    random_sum = np.sum(random_distances)
    hopkins_statistic = dataset_sum / (dataset_sum + random_sum)

    return hopkins_statistic
