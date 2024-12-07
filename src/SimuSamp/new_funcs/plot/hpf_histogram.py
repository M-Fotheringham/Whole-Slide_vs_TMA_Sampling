import pandas as pd
import matplotlib.pyplot as plt

def hpf_histogram(hpfs, bins, colour):

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
