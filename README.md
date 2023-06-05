# <div align="center"> Whole-Slide vs TMA Sampling </div>
#### <div align="center"> The contained code was developed to analyze whole-slide IHC data exported from HALO (Indica Labs) <div>
#### <div align="center">Contact: 16mrf6@queensu.ca</div>
  
  ## 1. Description
This is the description. Included are snippets of Python (v3.8.2) code that were used in concert to process the data and construct some of the figures discussed in "TISSUE MICROARRAY VS. WHOLE-SLIDE ANALYSIS OF CD8 IN NON-SMALL CELL LUNG CARCINOMA (Fotheringham et al., unpublished)".


![WhaleFig](documents/WhaleSlide.png)

  ### 1.1. Functions
  #### 1.1.1. load_sim_data
Tissue annotations in .annotations format (essentially an xml of listed vertices) and cell object data in the format of minimum and maximum x, y bounds. Annotations were digested and reassembled using Shapely. CD8+ cells were filtered and their centre points were calculated from the average of their bounds. 
  #### 1.1.2. n_core_sampler

![WhaleFig](documents/SimulatedSampling.png)
  
  #### 1.1.3. eta_counter
  
  #### 1.1.4. compute_hpfs
  
  #### 1.1.5. sample_hpfs
  
  #### 1.1.6. hpf_density_histograms
  
## 2. Instructions
### 2.1. Python Instructions
#### 2.1.1. Environment setup
These are the instructions.
