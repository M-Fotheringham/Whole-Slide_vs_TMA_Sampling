# <div align="center"> Whole-Slide vs TMA Sampling </div>
#### <div align="center"> The contained code was developed to analyze whole-slide IHC data exported from HALO (Indica Labs) <div>
#### <div align="center">Contact: 16mrf6@queensu.ca</div>
  
  ## 1. Description
This is the description. Included are snippets of Python (v3.8.2) code (I'll add the rest later) that were used in concert to process the data and construct some of the figures discussed in "Tissue Microarray vs. Whole-Slide Analysis of CD8 in Non-Small Cell Lung Carcinoma" (Fotheringham et al., *unpublished*).


![WhaleFig](documents/WhaleSlide.png)

  ### 1.1. Functions
  #### 1.1.1. load_sim_data
Tissue annotations (in .annotations format, essentially an xml of listed vertices) and cell object data (in the format of minimum and maximum x, y bounds) were exported from HALO. Annotations were digested and reassembled using Shapely (Polygon, MultiPolygon). CD8<sup>+</sup> cells were filtered and their centre points were calculated from the average of their bounds.
  
    Input:
  
 *Block*: This is the unique tissue block identifier;
  
 *parent_filepath*: Where the data is stored. The object data and annotations file are in folders for each tissue block, simply named *block*.
  
    Returns:
  
  *object_data*: dataframe of the coordinate data for CD8<sup>+</sup> cells;
  
  *plot_annotations*: Shapely Polygons representing the imported tissue regions.
  
  #### 1.1.2. n_core_sampler
 The reformatted annotations and cell object data are fed into n_core_sampler to simulate random TMA sampling of each tissue region. Coordinates within the range of the tissue bounds are randomly generated until a point is within the tissue polygon and the area of the simulated core generated from extending *circle_radius* from that coordinate point matches the tissue-specific criteria:
 Tumour cores must contain at least 50% tumour by area.
 IM cores must contain 80% => tumour => 10% and stroma => 10%.
  
    Input:
  
 *Block*: This is the unique tissue block identifier;
  
  *object_data*: CD8<sup>+</sup> cell coordinates derived from *load_sim_data*;
  
 *plot_annnotations*: Shapely polygons derived from *load_sim_data*;
  
 *number_of_tumour_regions*:
  
  *circle_radius*:
  
  *boundary_type*:
  
  *n_cores*:
  
  *microns_per_pixel*: Always 0.22715 µm/px.
 
    Returns:
 
  *sampling_results*: A dictionary containing the mean CD8<sup>+</sup> cell density, stdev, std error, tissue region, block, core radius, number of cores attempted, number of cores actually sampled, the mean CD8<sup>+</sup> cell count, and mean tissue area per sampling iteration.
  
  
  
![WhaleFig](documents/SimulatedSampling.png)
*A given sampling iteration using n_cores=10 of circle_radius=0.6 in the invasive margin (IM) and central tumour (CT) visualized with Matplotlib (external code).*
  
  
  #### 1.1.3. eta_counter
Ain't nobody got time for that.
  
  #### 1.1.4. compute_hpfs
  
  
  #### 1.1.5. sample_hpfs
  
  
  #### 1.1.6. hpf_density_histograms
  
  
## 2. Instructions
### 2.1. Python Instructions
#### 2.1.1. Environment setup
These are the instructions.
