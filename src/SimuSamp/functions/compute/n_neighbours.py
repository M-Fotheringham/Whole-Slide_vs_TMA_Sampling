import numpy as np
from sklearn.neighbors import NearestNeighbors


def edge_correction_weight(point, core, distance):
    # Get fraction of the buffer within the analysis area
    circle = point.buffer(distance)
    intersected_area = circle.intersection(core).area

    weight = circle.area / intersected_area

    return weight


def neighbours(
    gdf,
    k,
    microns_per_pixel=0.22715,
    edge_correction=False,
    core=None,
    distance=None,
):
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
    nn = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    distances, _ = nn.kneighbors(xy)

    dataset_distances = distances[:, 1:]

    if edge_correction:
        if core is None or distance is None:
            raise ValueError(
                "Core area and radius must be provided for edge correction."
            )

        weights = []
        for p in gdf.geometry:
            weights.append(edge_correction_weight(p, core, distance))

        weights = np.array(weights)

        scale = 1.0 / np.sqrt(weights)
        dataset_distances = dataset_distances * scale[:, None]

    dataset_distances = dataset_distances * microns_per_pixel

    return dataset_distances
