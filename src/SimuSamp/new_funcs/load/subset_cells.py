import geopandas as gpd


def subset(cell_df, anno_df, annotation):

    cells = gpd.sjoin(cell_df, anno_df, predicate="within")

    cells = cells[cells["layer"] == annotation].reset_index(drop=True)

    return cells
