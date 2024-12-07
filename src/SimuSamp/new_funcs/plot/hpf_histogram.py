import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt


def hpf_histogram(hpfs, cumulative=True, percentile=True, y_type="cell_count"):

    # Cumulative Density

    # Sort by density or count
    hpfs_sorted = hpfs.sort_values(f"{y_type}", ascending=False).reset_index(drop=True)
    
    # Calculate bins
    if percentile:
        # Get percentile bins
        bins = 100 - (hpfs_sorted.index / hpfs_sorted.index.max()) * 100
        p_label = " Percentile"
    else:
        bins = hpfs_sorted.index
        p_label = ""
    
    hpfs_sorted["bins"] = bins

    if cumulative:
        # Get cumulative sum of counts and area to calculate cumulative density
        hpfs_sorted["cum_cell_count"] = hpfs_sorted["cell_count"].cumsum()
        hpfs_sorted["cum_area"] = hpfs_sorted["area"].cumsum()
        hpfs_sorted["cum_density"] = hpfs_sorted["cum_cell_count"] / hpfs_sorted["cum_area"]
        hpfs_sorted["y_vals"] = hpfs_sorted[f"cum_{y_type}"]

        cum_label = "Cumulative "
    else:
        hpfs_sorted["y_vals"] = hpfs_sorted[f"{y_type}"]
        
        cum_label = ""
    
    # Plot the results
    plt.plot(hpfs_sorted["bins"], hpfs_sorted["y_vals"])
    # plt.gca().set_aspect("equal")
    plt.ylabel(f"{cum_label}{y_type}")
    plt.xlabel(f"HPF{p_label}")
    if percentile:
        plt.gca().invert_xaxis()
