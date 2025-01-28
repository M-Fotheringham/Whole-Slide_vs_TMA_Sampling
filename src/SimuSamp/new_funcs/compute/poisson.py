import pointpats
import geopandas as gpd
from shapely.geometry import Point


def poisson_cells(cell_df, anno_df, annotation, n_cells=None):
    """
    Args:
        cell_df (geopandas.GeoDataFrame): The cell data.
        anno_df (geopandas.GeoDataFrame): The annotation data.
        annotation (str): The annotation.
        n_cells (int): The number of points to sample.
    """

    if n_cells is None:

        cells = gpd.sjoin(cell_df, anno_df, predicate="within")

        cells = cells[cells["layer"] == annotation].reset_index(drop=True)

        n_cells = len(cells)

    anno = anno_df[anno_df["layer"] == annotation]["geometry"].iloc[0]

    random_points = pointpats.random.poisson(anno, size=n_cells)

    random_points = [Point(p) for p in random_points]

    random_gdf = gpd.GeoDataFrame(geometry=random_points)

    return random_gdf
