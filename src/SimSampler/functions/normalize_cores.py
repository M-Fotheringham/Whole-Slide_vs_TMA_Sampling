###
import pandas as pd
import numpy as np

###Import TMA Data
TMA = pd.read_csv("/Users/michael/OneDrive - Queen's University/Team TMA vs Whole Slide/Data/CD8 IHC TMA1314/TMA_CD8_SummaryData_MPIHC.csv")
#Drop unnecessary columns
TMA = TMA.drop(TMA.columns[np.r_[0:4, 7, 8, 10:29, 30:44, 45:49]], axis=1)
#Create density column
TMA['Density_(/mm)'] = (TMA['CD8 Cells'] / TMA['Tissue Area (\u03BCm\u00b2)'])*10**6
#Retrive columns from TMA Map
from _7_CalculatingAreaFromDepth_IMTMACores import TMA_MAP
TMA["Block"] = TMA_MAP["Block"]
TMA["Distance from Tumour"] = TMA_MAP["Distance from Tumour to Estimated Centre of Core (mm)"]
TMA["Revised Location"] = TMA_MAP["Revised Location"]
TMA["Estimated T Percentage"] = TMA_MAP["Estimated Tumour Percentages"]
#Drop invalid cores
TMA = TMA[TMA["Spot Valid"] == 1]
#Drop validity and case column
TMA = TMA.drop(TMA.columns[[0, 1]], axis=1)
#Extract IM
TMA_IM = TMA[TMA["Revised Location"] == "IM"].sort_values("Block").reset_index(drop=True)

###Import WS Data
WS = pd.read_csv("/Users/michael/OneDrive - Queen's University/Team TMA vs Whole Slide/Data/CD8 IHC Whole-Slide/WS_CD8_SummaryData_MPIHC.csv")
#Drop unnecessary rows
WS = WS[WS["Algorithm Name"] == "WS_Infiltration_Analysis"]
#Chop file extension off of Image Tag
WS['Image Tag'] = WS['Image Tag'].replace('_CD8.ndpi', '', regex=True)

#For reference
WS_col_indeces = pd.Series(list(WS.columns))

#Define density subset
density_data = WS[WS.columns[np.r_[1,36:76]]]
#Chop density bins into readable titles
density_data.columns = density_data.columns.str.replace("CD8 per mm\u00b2 \\[", "", regex=True)
density_data.columns = density_data.columns.str.replace("] of interface", "", regex=True)

#Correct for 0 areas, low areas
#Create boolean dataframe of desired areas (Image Tag needs to be True)
area_data = WS[WS.columns[np.r_[1,116:156]]]
area_data["Image Tag"] = [*[True]*len(area_data)]
area_data = area_data > 0.0002
#Chop column names to match density columns
area_data.columns = area_data.columns.str.replace("Band area in mm\u00b2 \\[", "", regex=True)
area_data.columns = area_data.columns.str.replace("] of interface", "", regex=True)
#Apply boolean dataframe to density dataframe
density_data = density_data[area_data]

#Normalize TMA IM cores to WS T bins (-1000 to 0)
stdevs_from_T_mean = []
for idx, block in enumerate(TMA_IM["Block"]):
    T_mean = np.nanmean(density_data[density_data["Image Tag"] == block].drop("Image Tag", axis=1))
    T_stdev = np.nanstd(density_data[density_data["Image Tag"] == block].drop("Image Tag", axis=1))
    norm_TMA_IM_den = (TMA_IM["Density_(/mm)"][idx] - T_mean) / T_stdev
    stdevs_from_T_mean.append(norm_TMA_IM_den)

TMA_IM["T_Norm_Density"] = stdevs_from_T_mean


del(TMA, TMA_MAP, WS, T_mean, T_stdev, WS_col_indeces, area_data, block, density_data, idx, norm_TMA_IM_den, stdevs_from_T_mean)
