"""Load TMA data for a given SampleID."""

import pandas as pd
from shapely import Polygon, MultiPolygon, Point, unary_union
import geopandas as gpd


def load_tma_data(sampleid, parent_filepath):
    """Load TMA data for a given SampleID."""

    filepath = parent_filepath.rsplit("/", 1)[0] + "/TMA_data"

    # Get map of core locations and other metadata =======================
    core_map = pd.read_excel(f"{filepath}/Core_Map.xlsx")

    core_map = core_map[core_map["SampleID"] == sampleid]

    core_map = core_map[core_map["QC2_Included"] == 1]

    core_map = core_map.reset_index(drop=True)
    # ====================================================================

    # Get single-cell data for cores in core_map =========================
    tma_data = pd.read_excel(f"{filepath}/TMA_ObjectData.xlsx")

    tma_data = tma_data[tma_data["Coords"].isin(core_map["Coords"])]

    tma_data = tma_data.reset_index(drop=True)

    # Compute central coordinates and convert to GeoDataFrame
    cell_coords = []
    for idx, row in tma_data.iterrows():
        x = (row["XMin"] + row["XMax"]) / 2
        y = (row["YMin"] + row["YMax"]) / 2
        cell_coords.append(Point(x, y))

    tma_data.loc[:, "geometry"] = cell_coords

    tma_data = tma_data.drop(["XMin", "XMax", "YMin", "YMax"], axis=1)

    tma_data = gpd.GeoDataFrame(tma_data, geometry="geometry")
    # ====================================================================

    # Still editing...
    # Get Annotation Data from xml-style files ===========================
    tma_anno = []
    for idx, label in enumerate(core_map["Annotation_File"]):
        with open(f"{filepath}/Annotations/{label}.annotations", "r") as file:
            lines = file.readlines()

        positive_region = []
        layers = []
        polygons = []
        vertices = []
        layer_name = []
        region_status = []
        for line in lines:
            if "<Annotation Name=" in line:
                # Start a new polygon, pull name from line
                layer_name = line.split('Name="')[1].split('" Visible=')[0]

            elif "<Region Type=" in line:
                # Ascertain if polygon is positive or negative
                region_status = line.split('NegativeROA="')[1]
                region_status = region_status.split('">')[0] == "0"
                vertices = []

            elif "<V X=" in line:
                # Extract the x, y coords from the line
                x = int(line.split('X="')[1].split('" Y=')[0])
                y = int(line.split('Y="')[1].split('" />')[0])
                vertices.append((x, y))

            elif "</Region>" in line:
                # End current polygon and add it to the list
                if len(vertices) >= 4:
                    polygons.append(Polygon(vertices))
                    layers.append(layer_name)
                    positive_region.append(region_status)

        annos = pd.DataFrame(
            {
                "Layer": layers,
                "Positive Region": positive_region,
                "Polygon": polygons,
                "Area": [Polygon(i).area for i in polygons],
            }
        )

        annos = annos.sort_values(
            ["Layer", "Positive Region", "Area"],
            ascending=[False, False, False],
        ).reset_index(drop=True)

        # Positive polygons for plotting
        positive_polygons = unary_union(
            MultiPolygon(
                annos[annos["Positive Region"] == 1]["Polygon"].tolist()
            ).buffer(0)
        )

        # Negative polygons for plotting
        negative_polygons = unary_union(
            MultiPolygon(
                annos[annos["Positive Region"] == 0]["Polygon"].tolist()
            ).buffer(0)
        )

        # Net polygons for analysis
        analysis_area = positive_polygons.difference(negative_polygons)

        anno_df = gpd.GeoDataFrame(
            {
                "layer": [
                    "positive_regions",
                    "negative_regions",
                    "analysis_area",
                ],
                "geometry": [
                    positive_polygons,
                    negative_polygons,
                    analysis_area,
                ],
            },
            geometry="geometry",
        )

        anno_df.loc[:, "Coords"] = core_map.loc[idx, "Coords"]

        tma_anno.append(anno_df)

    if len(tma_anno) == 0:
        raise ValueError(f"No annotated cores available for {sampleid}.")
    else:
        tma_anno = pd.concat(tma_anno).reset_index(drop=True)
    # =================================================================

    return core_map, tma_data, tma_anno
