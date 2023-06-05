###| This code performs random sampling (without replacement, ie not bootstrapping) by simulating n_cores of circle_radius n_iterations times each.
###| Updated June 1, 2023 by Michael Fotheringham

import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
from shapely.geometry import Polygon as shapelyPolygon
from shapely.geometry import Point, MultiPolygon
from shapely.ops import unary_union
from _0_Simulated_Sampling_Functions import load_sim_data, n_core_sampler

# Define global data and functions
microns_per_pixel = 0.22715
filepath = "C:/Users/.../Desktop/Simulated WS Sampling/block_data"
sampling_key = pd.read_excel("C:/Users/.../Desktop/Simulated WS Sampling/SamplingKey.xlsx")[["block", "number_of_tumour_regions"]]

# Function inputs
radius_list = [0.3, 0.5, 1]
sample_size = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
n_iterations = 50
total_sampling_results = []
for idx, block in enumerate(sampling_key["block"][1:]):
    number_of_tumour_regions = sampling_key["number_of_tumour_regions"][idx]
    object_data, plot_annotations = load_sim_data(block=block, parent_filepath=filepath)
    for boundary_type in ["Tumour", "IM"]:
        for circle_radius in radius_list:
            for n_cores in sample_size:
                for iteration in range(n_iterations):
                    # Update statements
                    print(f"Starting to randomly sample {block} {boundary_type} with {n_cores} {circle_radius}-mm (radius) cores. On iteration {iteration+1}/{n_iterations}.")
                    sampling_results = n_core_sampler(block=block, object_data=object_data, plot_annotations=plot_annotations, number_of_tumour_regions=number_of_tumour_regions, circle_radius=circle_radius, boundary_type=boundary_type, n_cores=n_cores)
                    total_sampling_results.append(sampling_results)


# Create beefy dataframe from list of dicts
sampling_data = pd.DataFrame(total_sampling_results)
# Save just in case
sampling_data.to_excel(r"C:\Users\...\Desktop\Simulated WS Sampling\WS_sampling_results\total_sampling_data.xlsx")
# Load from desktop, if returning
# sampling_data = pd.read_excel(r"C:\Users\...\Desktop\Simulated WS Sampling\WS_sampling_results\total_sampling_data.xlsx")
# Reformat for plotting
from ast import literal_eval
test_sample = sampling_data[sampling_data["Block"].isin(["01_F", "01_G", "02_G"])]
test_sample["Den_stdev"] = [int(literal_eval(x)[0]) for x in test_sample["Den_stdev"]]
test_sample["Den_sterr"] = [int(literal_eval(x)[0]) for x in test_sample["Den_sterr"]]

# For normalization
grouped = test_sample.groupby(["Block", "n_cores"])
mean_den = grouped["Density_n_mean"].transform('mean')
std_den = grouped["Density_n_mean"].transform('std')
test_sample["norm_den"] = (test_sample["Density_n_mean"] - mean_den)# / std_den
mean_err = grouped["Den_sterr"].transform('mean')
std_err = grouped["Den_sterr"].transform('std')
test_sample["norm_err"] = (test_sample["Den_sterr"] - mean_err)# / std_err
test_sample.fillna(0, inplace=True)

#####Create the plots
#To change fonts: plt.rcParams["font.family"] = ["font"]
plt.rcParams["font.family"] = ["Times New Roman"]
snsdata = test_sample[["norm_err", "n_cores"]]
ax = sns.boxplot(data=snsdata, x="n_cores", y="norm_err", palette="Blues", saturation=0.5, width=0.6, fliersize=4, linewidth=0.7)
ax.set_title("")
ax.set_ylabel("Mean-Adjusted Standard Error per Sample (cells/mm$^2$)")
ax.set_xlabel("Number of Cores per Sample")
ax.set_xticks([0, 1, 2, 3, 4], ["1", "2", "3", "4", "5"])
#each box is 6 lines, the median is the 5th
for x in [4, 10, 16, 22, 28]:
    ax.lines[x].set_color("royalblue")
for x in [0, 1, 2, 3, 4]:
    box = ax.patches[x]
    box.set_facecolor("white")

plt.savefig(r"C:\Users\...\Desktop\Simulated WS Sampling\cores_per_sample_50i", dpi=300)

