import numpy as np
from shapely import Point, STRtree
from SimuSamp.new_funcs.compute.hopkins_stat import hopkins_stat
from SimuSamp.new_funcs.compute.n_neighbours import neighbours


def sample_n_cores(
        sampleid,
        cell_data,
        primary_anno,
        secondary_anno=None,
        region="tumour",
        core_radius=0.5,
        n_cores=3,
        outer_im_anno=None,
        extended_partition=None,
        n_neighbours=False,
        microns_per_pixel=0.22715):
    """
    Args:
        sampleid (str): The sample ID.
        cell_data (geopandas.GeoDataFrame): The cell data.
        primary_anno (shapely.geometry.Polygon): The primary annotation.
        secondary_anno (shapely.geometry.Polygon): The secondary annotation.
        region (str): The region to sample from.
        core_radius (float): The core radius in microns.
        n_cores (int): The number of cores to sample.
        outer_im_anno (shapely.geometry.Polygon): The outer IM annotation.
        extended_partition (shapely.geometry.Polygon): The extended partition.
        n_neighbours (bool): Whether to compute the nearest neighbour distance.
        microns_per_pixel (float): The microns per pixel.
    """

    # Establish constants
    radius = core_radius * 1000 / microns_per_pixel
    diameter = radius * 2
    mm2_per_pixels2 = (microns_per_pixel / 1000) ** 2
    core_max_area = np.pi * radius ** 2

    # Collect core density information
    points = []
    points_tree = STRtree(points)
    cell_counts = []
    core_areas = []
    cell_densities = []
    n_den_stdev = []
    n_den_sterr = []
    nearest_neighbour = []
    hopkins_statistic = []

    # Ensure sufficient tissue for sampling
    if primary_anno.area > n_cores * core_max_area:

        # Iteratively sample points from annotation bounds and extend radius
        # This can be further optimized ======================================
        min_x, min_y, max_x, max_y = primary_anno.bounds
        for n in range(n_cores):

            attempt_n = 0
            while attempt_n < 200:

                x = np.random.uniform(min_x, max_x)
                y = np.random.uniform(min_y, max_y)
                centre = Point(x, y)

                # Ensure that there are no overlapping cores
                overlapping = False
                # for prev_point in points:
                #     if prev_point.distance(centre) < diameter:
                #         overlapping = True
                #         attempt_n += 1
                #         continue
                if any(points_tree.query(centre.buffer(diameter))):
                    attempt_n += 1
                    continue
        # ====================================================================

                if not overlapping and primary_anno.contains(centre):
                    core = centre.buffer(radius)
                    intersection = core.intersection(primary_anno)
                    proportion = intersection.area / core_max_area
                    if region == "tumour":
                        # Ensure core is at least 50% tumour
                        if proportion >= 0.5:
                            points.append(centre)
                            points_tree = STRtree(points)
                        else:
                            attempt_n += 1
                            continue
                    if region == "IM":
                        # Calculate tumour and stroma proportions
                        tum_inter = core.intersection(secondary_anno)
                        tum_prop = tum_inter.area / core_max_area
                        im_inter = core.intersection(outer_im_anno)
                        im_prop = im_inter.area / core_max_area
                        # Ensure 80% >= tumour >= 10%, and 10% >= stroma
                        tum_prop_filt = (tum_prop >= 0.10) & (tum_prop <= 0.8)
                        strom_prop_filt = im_prop >= 0.10
                        if tum_prop_filt & strom_prop_filt:
                            points.append(centre)
                            points_tree = STRtree(points)
                        else:
                            attempt_n += 1
                            continue

                    # Count cells in each core
                    core = core.intersection(extended_partition)
                    intersecting_cells = cell_data[cell_data.intersects(core)]
                    cell_count = len(intersecting_cells)
                    core_area = core.area * mm2_per_pixels2
                    cell_density = cell_count / core_area
                    # Collect for each core
                    cell_counts.append(cell_count)
                    core_areas.append(core_area)
                    cell_densities.append(cell_density)

                    if n_neighbours:
                        # Calculate Hopkins statistic for each core
                        h_stat = hopkins_stat(intersecting_cells, 0.05)
                        hopkins_statistic.append(h_stat)
                        # Calculate mean nearest neighbour distance
                        nn = neighbours(intersecting_cells, 1)
                        if nn is not np.nan:
                            nn = np.nanmean(nn)
                        nearest_neighbour.append(nn)
                    break

                else:
                    attempt_n += 1
                    continue

    # Compute averages, stdev, and sterror
    cores_sampled = len(points)
    if cores_sampled == 0:
        cell_density_n_mean = np.nan
        cell_counts_n_mean = np.nan
        core_areas_n_mean = np.nan
        n_den_stdev.append(np.nan)
        n_den_sterr.append(np.nan)
        hopkins_mean = np.nan
        nn_mean = np.nan
    else:
        cell_density_n_mean = np.nanmean(cell_densities)
        cell_counts_n_mean = np.nanmean(cell_counts)
        core_areas_n_mean = np.nanmean(core_areas)
        den_stdev = np.nanstd(cell_densities)
        n_den_stdev.append(den_stdev)
        n_den_sterr.append(den_stdev / np.sqrt(len(cell_densities)))
        hopkins_mean = np.nanmean(hopkins_statistic)
        nn_mean = np.nanmean(nearest_neighbour)

    # Store results in dict to grow in iterative loop outside of function
    sampling_results = {
        "Density_n_mean": cell_density_n_mean,
        "Den_stdev": n_den_stdev[0],
        "Den_sterr": n_den_sterr[0],
        "Region": region,
        "Sampleid": sampleid,
        "Radius": radius * microns_per_pixel,
        "n_cores": n_cores,
        "Cores_actually_sampled": cores_sampled,
        "Counts_n_mean": cell_counts_n_mean,
        "Areas_n_mean": core_areas_n_mean,
        "Nearest_neighbour_mean": nn_mean,
        "Hopkins_mean": hopkins_mean
        }

    if cores_sampled > 0:
        sampling_results["Density_top_core"] = max(cell_densities)
        sampling_results["Density_bottom_core"] = min(cell_densities)

    return sampling_results
