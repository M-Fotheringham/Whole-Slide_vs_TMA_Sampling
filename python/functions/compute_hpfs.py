### | These functions are used to digitally manipulate the WS and simulate TMA core sampling, and display the resulting data. |###
#     Last Updated: June 1, 2023 by Michael Fotheringham  
def compute_hpfs(block, object_data, plot_annotations, number_of_tumour_regions, boundary_type, hpf_dimension_microns, microns_per_pixel=0.22715):
    # Convert the cell coordinates to a GeoDataFrame
    geometry = [Point(x, y) for x, y in object_data["Cell Centre"]]
    cell_data_gdf = gpd.GeoDataFrame(object_data, geometry=geometry)
    # Create a spatial index for faster intersection checks
    # cell_data_sindex = cell_data_gdf.sindex
    # Assemble annotations
    tumour_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"].tolist()).buffer(0))
    outer_IM_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"][0:number_of_tumour_regions].tolist()).buffer(500 / microns_per_pixel))
    inner_IM_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"][0:number_of_tumour_regions].tolist()).buffer(-500 / microns_per_pixel))
    negative_polygons = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Partition Zone"]["Polygon"].tolist()).buffer(0))
    net_tumour = tumour_polygons.difference(negative_polygons)
    partition = outer_IM_polygons.difference(negative_polygons)
    net_IM = partition.difference(inner_IM_polygons)
    analysis_polygon = unary_union(MultiPolygon(plot_annotations[plot_annotations["Layer"] == "Tumour"]["Polygon"].tolist()).buffer(1000 / microns_per_pixel))
    analysis_polygon = analysis_polygon.difference(negative_polygons)
    # Compute HPFs and densities
    hpf_dimension_pixels = hpf_dimension_microns / microns_per_pixel
    mm2_per_pixels2 = (microns_per_pixel / 1000) ** 2
    computed_hpfs = []
    if boundary_type == "Tumour":
        start_time = datetime.datetime.now()
        print(f"Computing {block} {boundary_type} HPFs...")
        min_x, min_y, max_x, max_y = net_tumour.bounds
        partition_x = int((max_x - min_x) / hpf_dimension_pixels)
        partition_y = int((max_y - min_y) / hpf_dimension_pixels)
        px, py = np.linspace(min_x, max_x, partition_x), np.linspace(min_y, max_y, partition_y)
        total_iterations = partition_x * partition_y
        total_sampled = 0
        for x in range(len(px) - 1):
            for y in range(len(py) - 1):
                if ((total_sampled / 10) % 200 == 0) & (total_sampled > 0):
                    eta = eta_counter(start_time=start_time, total_sampled=list(range(total_sampled)), total_iterations=total_iterations)[0]
                    print(f"ETA: {eta}.")
                total_sampled += 1
                poly_hpf = shapelyPolygon([[px[x], py[y]], [px[x], py[y + 1]], [px[x + 1], py[y + 1]], [px[x + 1], py[y]]])
                hpf_area = poly_hpf.intersection(net_tumour).area
                if hpf_area / poly_hpf.area >= 0.50:
                    # Find the cells that intersect with poly_hpf
                    intersecting_cells = cell_data_gdf[cell_data_gdf.intersects(poly_hpf.intersection(net_tumour))]
                    # Count the intersecting cells
                    cell_count = len(intersecting_cells)
                    cell_den = cell_count / (hpf_area * mm2_per_pixels2)
                    hpf_area = hpf_area * mm2_per_pixels2
                    computed_hpfs.append({"block": block, "region": boundary_type, "hpf": poly_hpf, "density": cell_den, "cell_count": cell_count, "hpf_area": hpf_area})
    if boundary_type == "IM":
        print(f"Computing {block} {boundary_type} HPFs...")
        min_x, min_y, max_x, max_y = net_IM.bounds
        partition_x = int((max_x - min_x) / hpf_dimension_pixels)
        partition_y = int((max_y - min_y) / hpf_dimension_pixels)
        px, py = np.linspace(min_x, max_x, partition_x), np.linspace(min_y, max_y, partition_y)
        for x in range(len(px) - 1):
            for y in range(len(py) - 1):
                poly_hpf = shapelyPolygon([[px[x], py[y]], [px[x], py[y + 1]], [px[x + 1], py[y + 1]], [px[x + 1], py[y]]])
                hpf_area = poly_hpf.intersection(net_IM).area
                if hpf_area / poly_hpf.area >= 0.50:
                    # Find the cells that intersect with poly_hpf
                    intersecting_cells = cell_data_gdf[cell_data_gdf.intersects(poly_hpf.intersection(net_IM))]
                    # Count the intersecting cells
                    cell_count = len(intersecting_cells)
                    cell_den = cell_count / (hpf_area * mm2_per_pixels2)
                    hpf_area = hpf_area * mm2_per_pixels2
                    computed_hpfs.append({"block": block, "region": boundary_type, "hpf": poly_hpf, "density": cell_den, "cell_count": cell_count, "hpf_area": hpf_area})
    return computed_hpfs
