import numpy as np
from shapely.strtree import STRtree
import geopandas as gpd


def core_sample(random_points, radius, region, tum, outer_im):

    diameter = radius * 2
    core_max_area = np.pi * radius**2

    # Select random cores from random points ==================================
    # Shuffle df for randomization
    random_points = random_points.sample(frac=1).reset_index(drop=True)

    point_list = []
    tree = STRtree(point_list)
    for idx, point in random_points.iterrows():
        # If not within diameter of another core, add to list
        if not any(
            tree.query(point.geometry, predicate="dwithin", distance=diameter)
        ):

            core = point.geometry.buffer(radius)

            if region == "tumour":
                intersection = core.intersection(tum)
                proportion = intersection.area / core_max_area
                # Ensure core is at least 50% tumour
                if proportion >= 0.5:
                    point_list.append(point.geometry)
                    tree = STRtree(point_list)
                else:
                    continue

            if region == "IM":
                im_intersect = core.intersection(outer_im)
                im_prop = im_intersect.area / core_max_area
                tum_intersect = core.intersection(tum)
                tum_prop = tum_intersect.area / core_max_area
                # Ensure core is 80% >= tumour > 0%, and 10% >= stroma
                tum_filt = (tum_prop > 0.0) & (tum_prop <= 0.8)
                strom_filt = im_prop >= 0.10
                if tum_filt & strom_filt:
                    point_list.append(point.geometry)
                    tree = STRtree(point_list)
                else:
                    continue

    # Iteratively sample
    filtered_points = gpd.GeoDataFrame(geometry=point_list)
    filtered_points = filtered_points.drop(0, axis=0)  # Drop initial point

    return filtered_points
