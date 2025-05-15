# <div align="center"> Whole-Slide vs TMA Sampling </div>
#### <div align="center"> The contained code was developed to analyze whole-slide IHC data exported from HALO (Indica Labs) <div>
#### <div align="center">Contact: 16mrf6@queensu.ca</div>
  
  ## 1. Description
Included are snippets of Python code that were used in concert to process the data and construct some of the figures discussed in "Tissue Microarray vs. Whole-Slide Analysis of CD8 in Non-Small Cell Lung Carcinoma" (Fotheringham et al., *unpublished*).
The data required to perform the present analyses are available in Figshare at [Whole-Slide CD8 IHC Data and Tissue Annotation Data (HALO)](https://doi.org/10.6084/m9.figshare.28991807.v2).


![WhaleFig](docs/WhaleSlide.png)

  ### 1.1. Classes and Functions
  #### 1.1.1. SpatDat
The `SpatDat` `Class` is used to retrieve pre-computed annotation and cell data for a given tissue specimen. 
Tissue annotations (in .annotations format, essentially an xml of listed vertices) and cell object data (in the format of minimum and maximum x, y bounds) were exported from HALO. Annotations were digested and reassembled using Shapely (`Polygon`, `MultiPolygon`). CD8<sup>+</sup> cells were filtered and their centre points were calculated from the average of their bounds.
  
  - *sampleid*: This is the unique tissue block identifier;
  - *parent_filepath*: Where the data is stored. The object data and annotations file are in folders for each tissue block, simply named by the *sampleid*.
  
  **load_data**:
  Default method that imports the CD8 IHC and annotation data:
  
  - *object_data*: `DataFrame` of the coordinate data for CD8<sup>+</sup> cells;
  - *plot_annotations*: Shapely `Polygons` representing the imported tissue regions.

    Usage:

        spatdat = SpatDat("01_F", "../../../../Simulated WS Sampling/block_data")
        
        # E.g.
        01_f_object_data = spatdat.object_data
        01_f_annotation_data = spatdat.annotation_data
        01_f_sampleid = spatdat.sampleid
        01_f_parent_filepath = spatdat.parent_filepath
        01_f_filepath = spatdat.filepath

  <br><br>
  
  **subset_cells**:
  Method that retrieves a cell `GeoDataFrame` for cells in a given annotation.

  - *annotation*: the `str` name of the annotation of the desired cells.

    Usage

        spatdat = SpatDat("01_F", "../../../../Simulated WS Sampling/block_data")
        im_cells = spatdat.subset_cells("IM")

  <br><br>
  
  **poisson_distribution**:
  Method that creates random set of `n_cells` points within a given annotation using a Poisson point process. Access the points with `SpatDat.poisson_cells[annotation]`.

  - *annotation*: the `str` name of the annotation of the desired cells.
  - *n_cells*: the `int` number of cells to be simulated.

    Usage:

        spatdat = SpatDat("01_F", "../../../../Simulated WS Sampling/block_data")
        random_points = spatdat.poisson_distribution("IM", 500)
        # Can alternatively retrieve points from class object:
        random_IM_pints = spatdat.poisson_cells["IM"]

  <br><br>
  
  **subset_annotation**:
  Method that retrieves a `Polygon` `GeoDataFrame` for a given annotation.

  - *annotation*: the `str` name of the desired annotation.

    Usage:
  
        spatdat = SpatDat("01_F", "../../../../Simulated WS Sampling/block_data")
        im_annotation = spatdat.subset_annotation("IM")

  <br><br>
  
  **compute_fields**:
  Method that partitions the whole tissue into tiles `width_micron` wide and computes cell densities in IM and CT tiles. Access the tiles with `SpatDat.hpfs`.

  - *width_microns*: the `float` width of the square tiles that will partition the tissue.

    Usage:

        spatdat = SpatDat("01_F", "../../../../Simulated WS Sampling/block_data")
        spatdat.compute_fields("300.0")
  
        hpfs = spatdat.hpfs

  <br><br>
  
  #### 1.1.2. n_core_sampler
 The reformatted annotations and cell object data are fed into `n_core_sampler` to simulate random TMA sampling of each tissue region. Coordinates within the range of the tissue bounds are randomly generated until a point is within the tissue `Polygon` and the area of the simulated core generated from extending `core_radius` from that coordinate point matches the tissue-specific criteria:
 Tumour cores must contain at least 50% tumour by area.
 IM cores must contain 80% => tumour > 0% and stroma => 10%.
  
  Input:
  
 *spatdat*: the `SpatDat` `Class` object for a given specimen;
 
 *region*: the `str` region (either "IM" or "tumour") in which to simulate sampling.
  
 *core_radius*: the `float` radius of the desired simulated cores.

 *n_core_list*: a `list` of the number of cores to be sampled over which the function will iterate.

 *iterations*: the `int` number of replicate sampling iterations for each sample number in `n_core_list`.
  
  *microns_per_pixel*: The HALO resolution. Always 0.22715 µm/px.
 
  Returns:
 
  *sampling_results*: a `DataFrame` containing:
  
  - the mean CD8<sup>+</sup> cell density,
  - the std of CD8<sup>+</sup> cell densities,
  - std error of CD8<sup>+</sup> cell densities,
  - the tissue region,
  - the sampleid,
  - the core radius,
  - the number of cores attempted,
  - the number of cores actually sampled,
  - the given iteration,
  - the mean CD8<sup>+</sup> cell count,
  - the mean tissue area per sampling iteration,
  - the mean density of random whole-slide points,
  - the mean count of random whole-slide points,
  - the std of count of random whole-slide points,
  - the mean nearest-neighbour distance,
  - the std of nearest-neighbour distances,
  - a `list` of the core `Polygons`,
  - a `list` of `Points` from which cores could have been selected.
  
  
![WhaleFig](docs/sampling_example_horizontal_white.png)
*A given sampling iteration using `n_core_list=[5]` of `core_radius=0.3`, `core_radius=0.5`, and `core_radius=1.0` in the invasive margin (IM) and central tumour (CT) visualized with Matplotlib.*

<br>
  
## 2. Instructions
### 2.1. Download Whole-Slide and TMA CD8 IHC Data

The present whole-slide and TMA CD8 IHC data is available in Figshare at [Whole-Slide CD8 IHC Data and Tissue Annotation Data (HALO)](https://doi.org/10.6084/m9.figshare.28991807.v2). Create a parent directory to house the package and the required data. Download the block data and place it in a directory labelled Simulated WS Sampling, inside the parent directory.

### 2.2. Python Instructions
#### 2.2.1. Environment setup

In the command line, navigate to the directory containing the Simulate WS Sampling folder. 

Clone the repository with:

    git clone https://github.com/M-Fotheringham/Whole-Slide_vs_TMA_Sampling

Navigate to the SimuSamp folder in the command line and pip install with:

    cd SimuSamp
    pip install .

#### 2.2.2. Using the Package

The `Jupyter Notebooks` contained in SimuSamp/Src/SimuSamp/Scripts/ contain intended usages for the `SpatDat` `Class` and associated functions, methods described above.

Users can perform select operations without downloading the CD8 IHC data hosted in Figshare by running SimuSamp/Src/SimuSamp/Scripts/Test.ipynb.


