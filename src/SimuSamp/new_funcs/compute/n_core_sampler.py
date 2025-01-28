import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, MultiPoint
import geopandas as gpd
from SimuSamp.new_funcs.compute.hopkins_stat import hopkins_stat
from SimuSamp.new_funcs.compute.n_neighbours import neighbours
from SimuSamp.new_funcs.compute.core_sample import core_sample


def sample_n_cores(
    spatdat,
    region="tumour",
    core_radius=0.5,
    n_core_list=None,
    iterations=100,
    microns_per_pixel=0.22715,
):
    """
    Args:
    """

    # Establish constants =====================================================
    radius = core_radius * 1000 / microns_per_pixel
    mm2_per_pixels2 = (microns_per_pixel / 1000) ** 2

    if n_core_list is None:
        n_core_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    cells = spatdat.subset_cells("extended_partition")
    ext_partition = spatdat.subset_annotation("extended_partition")
    tum = spatdat.subset_annotation("tumour")
    outer_im = spatdat.subset_annotation("outer_IM")
    # =========================================================================

    # Generate random points within the region ================================
    random_points = spatdat.poisson_distribution(region)

    # Loop through conditions:

    den_mean = []
    den_stdev = []
    den_sterr = []
    region_list = []
    sampleid_list = []
    radius_list = []
    n_cores_list = []
    n_sampled_list = []
    iteration_list = []
    counts_n_mean = []
    areas_n_mean = []
    rand_n_den_mean = []
    rand_counts_n_mean = []
    rand_counts_n_stdev = []
    nearest_neighbour_n_mean = []
    nearest_neighbour_n_stdev = []
    hopkins_statistic_n_mean = []
    hopkins_statistic_n_stdev = []
    cores_list = []
    filtered_points_list = []

    for n_cores in n_core_list:
        for i in range(iterations):

            filtered_points = core_sample(
                random_points, radius, region, tum, outer_im
            )

            # Calculate density, spatial metrics per core =====================
            filtered_cores = filtered_points.buffer(radius)

            # Subsample n_cores
            n_samples = min([n_cores, len(filtered_cores)])

            sampled_gdf = filtered_cores.sample(n_samples)

            # Compile core metrics
            core_counts = []
            core_areas = []
            core_densities = []
            rand_counts = []
            rand_densities = []
            nearest_neighbour = []
            hopkins_statistic = []

            for geom in sampled_gdf:
                core = geom.intersection(ext_partition)
                core_area = core.area * mm2_per_pixels2

                # gdf because sjoin is faster than overlay
                core_df = gpd.GeoDataFrame(geometry=[core])

                core_cells = gpd.sjoin(cells, core_df, predicate="within")
                core_random_points = gpd.sjoin(
                    random_points, core_df, predicate="within"
                )

                cell_count = len(core_cells)
                core_counts.append(cell_count)

                core_density = cell_count / core_area
                core_areas.append(core_area)
                core_densities.append(core_density)

                rand_count = len(core_random_points)
                rand_counts.append(rand_count)
                rand_density = rand_count / core_area
                rand_densities.append(rand_density)

                core_nearest_neighbour = neighbours(core_cells, 1)
                if core_nearest_neighbour is not np.nan:
                    core_nearest_neighbour = np.nanmean(core_nearest_neighbour)
                nearest_neighbour.append(core_nearest_neighbour)

                core_hopkins_statistic = hopkins_stat(core, core_cells)
                hopkins_statistic.append(core_hopkins_statistic)

            den_mean.append(np.nanmean(core_densities))
            den_stdev.append(np.nanstd(core_densities))
            den_sterr.append(
                np.nanstd(core_densities) / np.sqrt(len(core_densities))
            )
            region_list.append(region)
            sampleid_list.append(spatdat.sampleid)
            radius_list.append(core_radius)
            n_cores_list.append(n_cores)
            n_sampled_list.append(n_samples)
            iteration_list.append(i + 1)
            counts_n_mean.append(np.nanmean(core_counts))
            areas_n_mean.append(np.nanmean(core_areas))
            if np.nansum(rand_counts) > 0:
                rand_n_den_mean.append(np.nanmean(rand_densities))
                rand_counts_n_mean.append(np.nanmean(rand_counts))
                rand_counts_n_stdev.append(np.nanstd(rand_counts))
            else:
                rand_n_den_mean.append(np.nan)
                rand_counts_n_mean.append(np.nan)
                rand_counts_n_stdev.append(np.nan)
            if np.nansum(core_counts) > 1:
                nearest_neighbour_n_mean.append(np.nanmean(nearest_neighbour))
                nearest_neighbour_n_stdev.append(np.nanstd(nearest_neighbour))
            else:
                nearest_neighbour_n_mean.append(np.nan)
                nearest_neighbour_n_stdev.append(np.nan)
            if np.nansum(core_counts) > 5:
                hopkins_statistic_n_mean.append(np.nanmean(hopkins_statistic))
                hopkins_statistic_n_stdev.append(np.nanstd(hopkins_statistic))
            else:
                hopkins_statistic_n_mean.append(np.nan)
                hopkins_statistic_n_stdev.append(np.nan)
            cores_list.append(MultiPolygon([x for x in sampled_gdf]))
            filtered_points_list.append(
                MultiPoint([x for x in filtered_points.geometry])
            )

    results_df = pd.DataFrame(
        {
            "density_mean": den_mean,
            "density_stdev": den_stdev,
            "density_sterr": den_sterr,
            "region": region_list,
            "sampleid": sampleid_list,
            "core_radius": radius_list,
            "n_cores": n_cores_list,
            "n_sampled": n_sampled_list,
            "iteration": iteration_list,
            "cell_counts_mean": counts_n_mean,
            "areas_mean": areas_n_mean,
            "random_density_mean": rand_n_den_mean,
            "random_counts_mean": rand_counts_n_mean,
            "random_counts_stdev": rand_counts_n_stdev,
            "nearest_neighbour_mean": nearest_neighbour_n_mean,
            "nearest_neighbour_stdev": nearest_neighbour_n_stdev,
            "hopkins_statistic_mean": hopkins_statistic_n_mean,
            "hopkins_statistic_stdev": hopkins_statistic_n_stdev,
            "cores": cores_list,
            "eligible_points": filtered_points_list,
        }
    )

    return results_df
