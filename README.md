# <div align="center"> Whole-Slide vs TMA Sampling </div>
#### <div align="center"> The contained code was developed to analyze whole-slide IHC data exported from HALO (Indica Labs) <div>
#### <div align="center">Contact: 16mrf6@queensu.ca</div>
  
  ## 1. Description
This is the description. Included are snippets of Python (v3.8.2) code that were used in concert to process the data and construct some of the figures discussed in "Tissue Microarray vs. Whole-Slide Analysis of CD8 in Non-Small Cell Lung Carcinoma" (Fotheringham et al., *unpublished*).


![WhaleFig](documents/WhaleSlide.png)

  ### 1.1. Functions
  #### 1.1.1. load_sim_data
Tissue annotations (in .annotations format, essentially an xml of listed vertices) and cell object data (in the format of minimum and maximum x, y bounds) were exported from HALO. Annotations were digested and reassembled using Shapely (Polygon, MultiPolygon). CD8<sup>+</sup> cells were filtered and their centre points were calculated from the average of their bounds. This function returns the coordinate data for CD8<sup>+</sup> cells (object_data) and the processed annotations (plot_annotations).
  
 **Block**: This is the unique tissue block identifier;
  
 **parent_filepath**: Where the data is stored. The object data and annotations file are in folders for each tissue block, simply named *block*.
  
  #### 1.1.2. n_core_sampler
 The reformatted annotations and cell object data are fed into n_core_sampler to simulate random TMA sampling of each tissue region. 
  
 **Block**: This is the unique tissue block identifier;
  
  **object_data**: CD8<sup>+</sup> cell coordinates derived from *load_sim_data*;
  
 **plot_annnotations**: Shapely polygons derived from *load_sim_data*;
  
 **number_of_tumour_regions**:
  
  **circle_radius**:
  
  **boundary_type**:
  
  **n_cores**:
  
  **microns_per_pixel**: Always 0.22715 µm/px.
 
 
  
![WhaleFig](documents/SimulatedSampling.png)
  
  #### 1.1.3. eta_counter
Ain't nobody got time for that.
  
  #### 1.1.4. compute_hpfs
  
  
  #### 1.1.5. sample_hpfs
  
  
  #### 1.1.6. hpf_density_histograms
  
  
## 2. Instructions
### 2.1. Python Instructions
#### 2.1.1. Environment setup
These are the instructions.
