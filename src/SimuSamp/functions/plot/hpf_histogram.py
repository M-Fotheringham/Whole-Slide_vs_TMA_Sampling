import pandas as pd
import matplotlib.pyplot as plt

def hpf_histogram(hpfs, bins, colour):


    # Edit ===================================================================
    if normalize == True:
        grouped = hpf_data.groupby(["block", "region"])
        den_mean = grouped["density"].transform("mean")
        den_std = grouped["density"].transform("std")
        hpf_data["density"] = (hpf_data["density"] - den_mean) / den_std
        upper_x_limit = -round(-max(hpf_data["density"]) // 1) * 1
        lower_x_limit = -round(-min(hpf_data["density"]) / 1) * 1
        # range requires integer intervals, multiplied everything by 100, calculated bins, then divide by 100 to retain float intervals
        bins = [*range(lower_x_limit*100, upper_x_limit*100, 100*(upper_x_limit - lower_x_limit) // 100)]
        bins = [bin / 100 for bin in bins]
        # plot_data = (plot_data - np.nanmean(plot_data)) / np.nanstd(plot_data)
    else:
        upper_x_limit = -round(-max(hpf_data["density"]) // 250) * 250
        lower_x_limit = 0
        bins = [*range(0, upper_x_limit*100, upper_x_limit*100 / 100)]
        
    # Grab parameters for sample size estimations (repeat for each region)
    plot_data = hpf_data[hpf_data["region"] == "Tumour"]["density"]
    plot_data = plot_data.sort_values()
    print(f"Std: {np.std(plot_data)}, mean: {np.mean(plot_data)}, variance: {(np.std(plot_data))**2}, median: {np.median(plot_data)}.")
    # Use these parameters on the website Wei provided- copy the equation and perform the equation

    # Distribution per region, per block
    for block in sampling_key["block"].unique():
            for boundary_type in ["Tumour", "IM"]:
                hpf_density_histograms(block=block, hpf_data=hpf_data, boundary_type=boundary_type, upper_x_limit=upper_x_limit, bins=bins, normalize=True)

    # Total distribution per region
    for boundary_type in ["Tumour", "IM"]:
        hpf_density_histograms(hpf_data=hpf_data, boundary_type=boundary_type, upper_x_limit=upper_x_limit, bins=bins, normalize=False)
    # =======================================================================


    # Plot
    plt.hist(hpfs, bins=bins, edgecolor="black", color=colour)
             
    # Formattting
    plt.ylabel("HPF Frequency")
    # plt.xlabel("Density (cells/mm$^2$)")
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.gca().tick_params(width=2)
