import pandas as pd
import geopandas as gpd


def sample_n_cores():

    # Edit ==================================================================
    # Randomly select core points
    radius = circle_radius * 1000 / microns_per_pixel
    diameter = radius * 2
    points = []
    cell_density = []
    cell_counts = []
    core_areas = []
    n_den_stdev = []
    n_den_sterr = []
    mm2_per_pixels2 = (microns_per_pixel / 1000) ** 2
    if boundary_type == "Tumour":
        if net_tumour.area > n_cores * np.pi * radius ** 2:
            min_x, min_y, max_x, max_y = net_tumour.bounds
            for i in range(n_cores):
                attempt_n = 0
                while attempt_n < 200:
                    x = np.random.uniform(min_x, max_x)
                    y = np.random.uniform(min_y, max_y)
                    centre = Point(x, y)
                    # Ensure that there are no overlapping cores
                    overlapping = False
                    for prev_point in points:
                        if prev_point.distance(centre) < diameter:
                            overlapping = True
                            attempt_n += 1
                            break
                    if not overlapping and net_tumour.contains(centre):
                        # Ensure core is at least 50% tumour
                        circle = centre.buffer(radius)
                        intersection = circle.intersection(net_tumour)
                        proportion = intersection.area / circle.area
                        if proportion >= 0.5:
                            points.append(centre)
                            # Find the cells that intersect with poly_hpf
                            intersecting_cells = cell_data_gdf[cell_data_gdf.intersects(circle.intersection(net_tumour))]
                            # Count the intersecting cells
                            cell_count = len(intersecting_cells)
                            cell_density.append(cell_count / (circle.intersection(analysis_polygon).area * mm2_per_pixels2))
                            cell_counts.append(cell_count)
                            core_areas.append(circle.intersection(analysis_polygon).area * mm2_per_pixels2)
                            break
                        else: attempt_n += 1
                    else:
                        attempt_n += 1
    if boundary_type == "IM":
        if net_IM.area > n_cores * np.pi * radius ** 2:
            min_x, min_y, max_x, max_y = net_IM.bounds
            for i in range(n_cores):
                attempt_n = 0
                while attempt_n < 200:
                    x = np.random.uniform(min_x, max_x)
                    y = np.random.uniform(min_y, max_y)
                    centre = Point(x, y)
                    # Ensure that there are no overlapping cores
                    overlapping = False
                    for prev_point in points:
                        if prev_point.distance(centre) < diameter:
                            overlapping = True
                            attempt_n += 1
                            break
                    if not overlapping and net_IM.contains(centre):
                        # Ensure 80% >= tumour >= 10%, and 10% >= stroma
                        circle = centre.buffer(radius)
                        intersection = circle.intersection(net_tumour)
                        proportion = intersection.area / circle.area
                        proportion_stroma = circle.intersection(partition.difference(net_tumour)).area / circle.area
                        if (proportion >= 0.10) & (proportion <= 0.80) & (proportion_stroma > 0.10):
                            points.append(centre)
                            # Find the cells that intersect with poly_hpf
                            intersecting_cells = cell_data_gdf[cell_data_gdf.intersects(circle.intersection(net_IM))]
                            # Count the intersecting cells
                            cell_count = len(intersecting_cells)
                            cell_density.append(cell_count / (circle.intersection(analysis_polygon).area * mm2_per_pixels2))
                            cell_counts.append(cell_count)
                            core_areas.append(circle.intersection(analysis_polygon).area * mm2_per_pixels2)
                            break
                        else: attempt_n += 1
                    else:
                        attempt_n += 1
    # Compute averages and stdev and sterror
    cores_sampled = len(cell_density)
    if cores_sampled == 0:
        cell_density_n_mean = np.nan
        cell_counts_n_mean = np.nan
        core_areas_n_mean = np.nan
        n_den_stdev.append(np.nan)
        n_den_sterr.append(np.nan)
    else:
        cell_density_n_mean = np.nanmean(cell_density)
        cell_counts_n_mean = np.nanmean(cell_counts)
        core_areas_n_mean = np.nanmean(core_areas)
        n_den_stdev.append(np.nanstd(cell_density))
        n_den_sterr.append(np.nanstd(cell_density) / np.sqrt(len(cell_density)))

    sampling_results = {"Density_n_mean": cell_density_n_mean, "Den_stdev": n_den_stdev, "Den_sterr": n_den_sterr, "Boundary": boundary_type, "Block": block, "Radius": radius, "n_cores": n_cores, "Cores_actually_sampled": cores_sampled, "Counts_n_mean": cell_counts_n_mean, "Areas_n_mean": core_areas_n_mean}

    if cores_sampled > 0:
        sampling_results["Density_top_core"] = max(cell_density)

    return sampling_results

    # =======================================================================

