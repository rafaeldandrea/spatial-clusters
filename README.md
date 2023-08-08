# spatial-clusters
Code and result files for "A spatial signal of niche differentiation in tropical forests"

<br>

*Citation for the code and data repository:

<br>

*Authors:

Mihir Umarani:   mihir.umarani@gmail.com

John Wang:       johnwang@g.harvard.edu

James O'Dwyer:   jodwyer@illinois.edu

Rafael D'Andrea: rafael.dandrea@stonybrook.edu

<br>

*Study summary


Here, we quantify the signatures of niche differentiation in tree species in tropical forest dynamic plots.
We demontrate a novel method that quantifies the strength of niche differentiation based on spatial 
clustering tendencies between pairs of species. We also characterize these species clusters in terms of their 
association with soil nutrients. We find a strong spatial clustering of species;
further, these clusters are strongly associated with the spatial distribution of soil nutrients, indicating 
that the different groups correspond to different niche strategies for local soil conditions.

<br>

*Authors responsible for writing code

Rafael D'Andrea

Mihir Umarani

John Wang

<br>

*Description of the contents of the repository


Folders: 


1. **codes**: This folder contains all the codes written in R (R version 4.0.5)
   
 	1.1. bci.R_ and laplanada.R_ run the data analysis (cluster analysis, soil nutrient analysis etc.) of the BCI and La Planada
   	   forest dynamic plot datasets, respectively. _bci.R_ also contains the recruitment analysis across censuses and species trait analysis.
   
   	1.2. *clustering_functions_rann.R* and *clustering_functions_rann_lap.R*: These file contain necessary supporting functions for the data analysis. 
	   These files are sourced into _bci.R_ and _laplanada.R_ to call the functions therein. There is no need to execute their codes independently.
   
	1.3. _Validation_analysis.R_: This file contains the code required for the analysis of the validation of our method of quantifying niche differentiation. 
	   This file requires no external input and can be executed independently.

2. **Final datasets**: This folder contains 2 sub-folders containing the results of the clustering analysis and kernel density estimation of spatial clusters for BCI and La Planada plots.
   
	2.1. Files named ___cluster_analysis.rds show the results of clustering analysis and have the following columns:
	   sp: Species code
	   group: Cluster (estimated through the algorithm) to which the species belongs. 
	   algorithm: Name of the alogorithm to calculate modularity of the clustering network.
	   weighted: (Yes/No) Whether the edges in the network were weighted.
	   number_of_groups: Total number of distinct clusters found.
	   d_cutoff (model parameter): Distance cutoff (in meters) used to identify the 'neighbor' trees.
	   seed: seed values used for randomization of the pairwise distance matrix for species. 0 indicates no changes.
   
	2.2. Files named ___kde_full.rds show the results of kernel density estimation of distinct spatial cluster of 
	   species (basically shows the density of each cluster in every cell of a map) and have the following columns:
	   census: (if there are multiple censuses) 
	   d_cutoff (model parameter): Distance cutoff (in meters) used to identify the 'neighbor' trees
	   group: Cluster group inferred from the clustering analysis
	   x and y: x and y coordinates of a cell (10m X 10m)
	   density: Estimated density of the individuals from a given species cluster from the kernel density estimation.
	   soiltype: Cluster from the census 1 (reference cluster) corresponding to the group.
	   fdp: Forest dynamic plot.

<br>

Coding software and packages:
All the code was executed in R, version (4.3.1).

Essential packages used (versions):
tidyverse(2.0.0)
openxlsx(4.2.5.2)
magrittr(2.0.3)
furrr(0.3.1)
readxl(1.4.3)
parallelDist(0.2.6)
igraph (1.5.0.1)
RANN (2.6.1)
FactoClass (1.2.7)
C50(0.1.8)
caret (6.0-94)
sparr(2.3-10)

pcaMethods(1.92.0)
Use the following code to install:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("pcaMethods")
		   

RandomFields (3.3.14)
This package will need to be installed from CRAN archive.

Step 1. Install the two dependencies first: 
	a)RandomFieldsUtils (version 1.2.5)
	This will need to be installed from the archive source as well:
	Use this code: 

	require("devtools")
	urls="https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.2.5.tar.gz"
	install.packages(urls,repos=NULL,type='source')

	b) sp
	install.packages("sp")

Step 2:Install the RandomFields package from CRAN archive
	
	require("devtools")
	urls="https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.2.5.tar.gz"
	install.packages(urls,repos=NULL,type='source')