###
def sample_hpfs(block, plot_annotations, hpf_data, percent_tissue_list, boundary_type, n_iterations, number_of_tumour_regions, microns_per_pixel=0.22715):
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
    # Refine dataframe
    hpf_data = hpf_data[(hpf_data["block"] == block) & (hpf_data["region"] == boundary_type)]
    mm2_per_pixels2 = (microns_per_pixel / 1000) ** 2
    sampled_hpfs = []
    if boundary_type == "Tumour":
        # Sample n_iterations per area
        for percent in percent_tissue_list:
            # Determine area of tissue by percent
            x_area = (net_tumour.area * mm2_per_pixels2) * percent / 100
            for iteration in range(n_iterations):
                print(f"Sampling {percent}% of {block} {boundary_type}. On iteration {iteration+1}/{n_iterations}.")
                shuffled_data = hpf_data.sample(frac=1)
                sampled_area = 0
                # Randomly select rows until area meets percentage
                for idx, row in shuffled_data.iterrows():
                    if sampled_area + row['hpf_area'] <= x_area:
                        hpf_cell_count = row["cell_count"]
                        hpf_area = row["hpf_area"]
                        sampled_hpfs.append({"block": block, "region": boundary_type, "percent_sampled": percent, "cell_count": hpf_cell_count, "hpf_area": hpf_area})
                        sampled_area += row['hpf_area']
                    else:
                        break
    if boundary_type == "IM":
        # Sample n_iterations per area
        for percent in percent_tissue_list:
            # Determine area of tissue by percent
            x_area = (net_IM.area * mm2_per_pixels2) * percent / 100
            for iteration in range(n_iterations):
                print(f"Sampling {percent}% of {block} {boundary_type}. On iteration {iteration+1}/{n_iterations}.")
                shuffled_data = hpf_data.sample(frac=1)
                sampled_area = 0
                # Randomly select rows until area meets percentage
                for idx, row in shuffled_data.iterrows():
                    if sampled_area + row['hpf_area'] <= x_area:
                        hpf_cell_count = row["cell_count"]
                        hpf_area = row["hpf_area"]
                        sampled_hpfs.append({"block": block, "region": boundary_type, "percent_sampled": percent, "cell_count": hpf_cell_count, "hpf_area": hpf_area})
                        sampled_area += row['hpf_area']
                    else:
                        break
    return sampled_hpfs
