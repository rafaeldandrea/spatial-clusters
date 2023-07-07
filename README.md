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


1. codes:This folder contains all the codes written in R (R version 4.0.5)
 	1. *clustering_functions_rann.R* and *clustering_functions_rann_lap.R*: This file contains necessary supporting functions for the data analysis. 
	   These files are sourced into the rest of the R files to call the functions in them. There is no need to execute their codes independently.
	2. _bci.R_ and _laplanada.R_ contain the data analysis (cluster analysis, soil nutrient analysis etc.) of the BCI and La Planada
   	   forest dynamic plot datasets, respectively. _bci.R_ also contains the recruitment analysis across censuses and species trait analysis.
	3. _Validation_analysis.R_: This file contains the code required for the analysis of the validation of our method of quantifying niche differentiation. 
	   This file requires no external input and can be executed independently.

2. Final datasets: This folder contains 2 sub-folders containing the results of the clustering analysis and kernel density estimation of 
		 spatial clusters for BCI and La Planada plots.

		1. Files named ___cluster_analysis.rds show the results of clustering analysis and have the following columns:
		   sp: Species code
		   group: Cluster to which the species belongs 
		   algorithm: Name of the alogorithm to calculate modularity of the clustering network
		   weighted: Whether the edges in the network were weighted
		   number_of_groups: Total number of distinct clusters found
		   d_cutoff (model parameter): Distance cutoff (in meters) used to identify the 'neighbor' trees
		   seed: seed values used for randomization of the pairwise distance matrix for species. 0 indicates no changes. 
		2. Files named ___kde_full.rds show the results of kernel density estimation of distinct spatial cluster of 
		   species (basically shows the density of each cluster in every cell of a map) and have the following columns:
		   census: (if there are multiple censuses) 
		   d_cutoff (model parameter): Distance cutoff (in meters) used to identify the 'neighbor' trees
		   group: Cluster group inferred from the clustering analysis
		   x and y: x and y coordinates of a cell (10m X 10m)
		   density: Estimated density of the individuals from a given species cluster from the kernel density estimation.
		   soiltype: Cluster from the census 1 (reference cluster) corresponding to the group.
		   fdp: Forest dynamic plot.

			
Coding software and packages:
All the code was executed in R, version (4.0.5).
Essential packages used (versions)-
RANN (2.6.1)
igraph (1.3.0)
sparr(2.2-15)
		   




