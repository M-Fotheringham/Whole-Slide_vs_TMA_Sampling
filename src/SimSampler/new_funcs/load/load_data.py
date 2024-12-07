import pandas as pd
from shapely import Polygon, MultiPolygon, Point, unary_union
import geopandas as gpd


def load_data(filepath, cell_name="CD8", microns_per_pixel=0.22715):
    """
    Loads data exported from HALO into Shapely objects
    """
    # Get Object Data ========================================================
    object_data = pd.read_csv(f"{filepath}/Object_Data.csv")

    object_data = object_data[
        ["Analysis Region", "XMin", "XMax", "YMin",
            "YMax", f"{cell_name}"]
            ]

    object_data = object_data[
        (object_data["Analysis Region"] == "Partition Zone")
        & (object_data[f"{cell_name}"] == 1)
        ]

    cell_coords = []
    for idx, row in object_data.iterrows():
        x = (row["XMin"] + row["XMax"]) / 2
        y = (row["YMin"] + row["YMax"]) / 2
        cell_coords.append(Point(x, y))

    object_data.loc[:, "geometry"] = cell_coords

    object_data = object_data.drop(
        ["XMin", "XMax", "YMin", "YMax"], axis=1)

    object_data = gpd.GeoDataFrame(object_data, geometry="geometry")

    object_data = object_data.reset_index(drop=True)
    # ========================================================================

    # Get Annotation Data from xml-style file ================================
    with open(f"{filepath}/Annotations.annotations", "r") as file:
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
            region_status = line.split('NegativeROA="')[1].split('">')[0] == '0'
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
        {"Layer": layers,
            "Positive Region": positive_region,
            "Polygon": polygons,
            "Area": [Polygon(i).area for i in polygons]
        }
        )

    layer_filt = ["Tumour", "Peritumoral Zone", "Partition Zone"]
    annos = annos[annos["Layer"].isin(layer_filt)]

    annos = annos.sort_values(
        ["Layer", "Positive Region", "Area"],
        ascending=[False, False, False]).reset_index(drop=True)
    # ========================================================================

    # Compile relevant annotation layers =====================================
    annos = annos[
        ((annos["Layer"] == "Partition Zone")
            & (annos["Positive Region"] == 0))
        |
        ((annos["Layer"] == "Tumour")
            & (annos["Positive Region"] == 1))
        ]

    annos = annos.reset_index(drop=True)

    # Positive polygons for plotting
    tumour_polygons = unary_union(
        MultiPolygon(
            annos[annos["Layer"] == "Tumour"]["Polygon"].tolist()
            ).buffer(0)
            )
    outer_im_polygons = unary_union(
        MultiPolygon(
            annos[(annos["Layer"] == "Tumour")]["Polygon"].tolist()
            ).buffer(500 / microns_per_pixel)
            )
    inner_im_polygons = unary_union(
        MultiPolygon(
            annos[annos["Layer"] == "Tumour"]["Polygon"].tolist()
            ).buffer(-500 / microns_per_pixel)
            )
    # Negative polygons for plotting
    negative_polygons = unary_union(
        MultiPolygon(
            annos[annos["Layer"] == "Partition Zone"]["Polygon"].tolist()
            ).buffer(0)
            )

    # Net polygons for analysis
    net_tumour = tumour_polygons.difference(negative_polygons)

    net_outer_im = outer_im_polygons.difference(tumour_polygons)
    net_outer_im = net_outer_im.difference(negative_polygons)

    net_inner_im = tumour_polygons.difference(inner_im_polygons)
    net_inner_im = net_inner_im.difference(negative_polygons)

    partition = outer_im_polygons.difference(negative_polygons)
    net_im = partition.difference(inner_im_polygons)

    # Polygon with additional buffer for simulations
    analysis_polygon = unary_union(
        MultiPolygon(
            annos[annos["Layer"] == "Tumour"]["Polygon"].tolist()
            ).buffer(1000 / microns_per_pixel)
            )
    analysis_polygon = analysis_polygon.difference(negative_polygons)

    annotation_data = gpd.GeoDataFrame(
        {"layer": [
            "tumour",
            "IM",
            "inner_IM",
            "outer_IM",
            "partition",
            "extended_partition",
            "exclusions",
            "tumour_hull",
            "inner_IM_hull",
            "outer_IM_hull"
            ],
         "geometry": [
             net_tumour,
             net_im,
             net_inner_im,
             net_outer_im,
             partition,
             analysis_polygon,
             negative_polygons,
             tumour_polygons,
             inner_im_polygons,
             outer_im_polygons]},
        geometry="geometry")
    # ========================================================================

    return object_data, annotation_data
