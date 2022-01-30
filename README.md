# BCI_Photogrammetry
This project quantifies canopy disturbances between 2015 and 2020 across Barro Colorado Island, Panama, using drone photogrammetry. Disturbance rates are compared to potential drivers of spatial variation (soils, topography, forest age), and to patterns of standing forest structure in 2009 lidar data.

Methods and results are described in the manuscript "Soils and topography control natural disturbance rates and thereby forest structure in a lowland tropical landscape", Ecology Letters (2020) by K.C Cushman, Matteo Detto, Milton Garcia, and Helene C. Muller-Landau. https://doi.org/10.1111/ele.13978

![EB-03-03519_0053_0122](https://user-images.githubusercontent.com/15330340/149802977-a1c868a7-c12e-4c1b-963c-48f641a674ab.JPG)

Data supporting analyses are archived through Smithsonian Institution Figshare upon publication (DOI will become active when published: 10.25573/data.17102600). A zipped copy of this repository at the time of publication will also be archived in the Figshare data repository. Here, code scripts reference all data files from the folders in which they are organized on Figshare:

- Code_AlignDroneData: R script for tiling raw point cloud data, and .bat files for aligning point cloud tiles in CloudCompare's command line tools.

- Code_GapSizeFrequency: R scripts for fitting and plotting gap size frequency data from gap rasters and shapefiles.

- Code_INLA: R scripts for configuring data for INLA models, running INLA models, and analyzing INLA results.

- Code_MakeFigures: R scripts for making main and supplemental figures.

- Code_ProcessHeightData: R scripts for making canopy height rasters from point cloud data, defining gap rasters/polygons from canopy height data, and smoothing digital elevation model (DEM) topography data from 2009 lidar.

- Code_QAQC: R scripts to make data quality masks based on cloud and photogrammetric reconstruction quality, and to find the optimal height correction for 2015 data with lower image overlap.
