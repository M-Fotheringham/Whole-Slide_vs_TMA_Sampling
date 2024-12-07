from shapely import Polygon
import numpy as np
import geopandas as gpd


def compute_hpfs(
        sampleid,
        cells,
        layer,
        annos,
        width_microns,
        microns_per_pixel=0.22715
        ):
    """
    Generate a partition of high-power fields and compute cell counts.
    """

    # Establish constants
    width_pixels = width_microns / microns_per_pixel
    mm2_per_pixels2 = (microns_per_pixel / 1000) ** 2
    partition = annos[annos["layer"] == layer]["geometry"]
    tumour = annos[annos["layer"] == "net_tumour"]["geometry"]
    im = annos[annos["layer"] == "net_IM"]["geometry"]

    # Generate HPFs from annotation bounds
    min_x, min_y, max_x, max_y = partition.bounds
    partition_x = int((max_x - min_x) / width_pixels)
    partition_y = int((max_y - min_y) / width_pixels)
    px = np.linspace(min_x, max_x, partition_x)
    py = np.linspace(min_y, max_y, partition_y)

    computed_hpfs = []
    for x in range(len(px) - 1):
        for y in range(len(py) - 1):

            poly_hpf = Polygon(
                [[px[x], py[y]], [px[x], py[y + 1]],
                 [px[x + 1], py[y + 1]], [px[x + 1], py[y]]])

            hpf_inter = poly_hpf.intersection(partition)
            hpf_area = hpf_inter.area * mm2_per_pixels2

            if hpf_area / poly_hpf.area >= 0.5:
                # Compute cell count
                intersecting_cells = cells[cells.intersects(hpf_inter)]
                cell_count = len(intersecting_cells)
                cell_den = cell_count / hpf_area

                # Estimate HPF region from layer area
                tum_frac = poly_hpf.intersection(tumour).area / hpf_area
                im_frac = poly_hpf.intersection(im).area / hpf_area
                if tum_frac > im_frac:
                    region = "tumour"
                else:
                    region = "im"

                # Store HPFs >= 50% area
                computed_hpfs.append(
                    {
                        "sampleid": sampleid,
                        "region": region,
                        "geometry": poly_hpf,
                        "density": cell_den,
                        "cell_count": cell_count,
                        "area": hpf_area
                        }
                        )

    computed_hpfs = gpd.GeoDataFrame(computed_hpfs, geometry="geometry")

    return computed_hpfs