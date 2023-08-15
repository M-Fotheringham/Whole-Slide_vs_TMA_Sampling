###
def n_core_sampler(block, object_data, plot_annotations, number_of_tumour_regions, circle_radius, boundary_type, n_cores, microns_per_pixel=0.22715):
    # if microns_per_pixel is None or not isinstance(microns_per_pixel, (int, float)):
    #     raise ValueError("microns_per_pixel is either missing or not a number. Should be 0.22715.")
    # if number_of_tumour_regions is None or not isinstance(number_of_tumour_regions, (int, float)):
    #     raise ValueError("number_of_tumour_layers is either missing or not a number.")
    # if circle_radius is None or not isinstance(circle_radius, (int, float)):
    #     raise ValueError("circle_radius is either missing or not a number.")
    # if boundary_type is None or not isinstance(boundary_type, str):
    #     raise ValueError("boundary_type is either missing or not a string. Must be 'Tumour' or 'IM'.")
    # if boundary_type not in ["Tumour", "IM"]:
    #     raise ValueError("boundary_type must be 'Tumour' or 'IM'.")
    # if n_cores is None or not isinstance(n_cores, int):
    #     raise ValueError("n_cores is either missing or not a number.")
    geometry = [Point(x, y) for x, y in object_data["Cell Centre"]]
    cell_data_gdf = gpd.GeoDataFrame(object_data, geometry=geometry)
    # Create a spatial index for faster intersection checks
    # cell_data_sindex = cell_data_gdf.sindex
    # Combine polygons, add to lists to make iterable for MultiPolygon
    tumour_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"].tolist()).buffer(0))
    outer_IM_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"][0:number_of_tumour_regions].tolist()).buffer(500 / microns_per_pixel))
    inner_IM_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"][0:number_of_tumour_regions].tolist()).buffer(-500 / microns_per_pixel))
    negative_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Partition Zone"]["Polygon"].tolist()).buffer(0))
    net_tumour = tumour_polygons.difference(negative_polygons)
    partition = outer_IM_polygons.difference(negative_polygons)
    net_IM = partition.difference(inner_IM_polygons)
    analysis_polygon = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"].tolist()).buffer(1000 / microns_per_pixel))
    analysis_polygon = analysis_polygon.difference(negative_polygons)
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
