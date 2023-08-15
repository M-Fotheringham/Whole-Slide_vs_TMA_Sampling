###
def hpf_density_histograms(hpf_data, boundary_type, upper_x_limit, bins, block="total", microns_per_pixel=0.22715, normalize=False):
    if block == "total":
        plot_data = hpf_data[hpf_data["region"] == boundary_type]["density"]
    else:
        plot_data = hpf_data[(hpf_data["block"] == block) & (hpf_data["region"] == boundary_type)]["density"]
    plt.rcParams["font.family"] = ["Arial"]
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    colour = "lightcoral"
    if boundary_type == "IM":
        colour = "lightgreen"
    plt.hist(plot_data, bins=bins, edgecolor="black", color=colour)
    plt.title(f"{block} WS HPF {boundary_type} CD8 Density Distribution")
    if block =="total":
        plt.ylim(0,3000)
        plt.legend(["Tumour", "IM"], frameon=False)
    if normalize == True:
        plt.xlabel("Normalized CD8 Density (Stdev from Mean)")
    else:
        plt.xlabel("Density (cells/mm$^2$)")
    plt.ylabel("Frequency")
    # buffer = (upper_x_limit - lower_x_limit) * 0.05
    # plt.xlim(-buffer, upper_x_limit + buffer)
    # Graphpadify
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.gca().tick_params(width=2)
    plt.savefig(f"C:/Users/labuser/Desktop/Simulated WS Sampling/HPF_Distributions/{block}_{boundary_type}_HPF_distribution.png", dpi=300)
    plt.cla()
    plt.clf()
