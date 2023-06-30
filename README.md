# spatial-clusters
Code and result files for "A spatial signal of niche differentiation in tropical forests"


*Citation for the code and data repository:


*Authors:

Mihir Umarani:   mihir.umarani@gmail.com
John Wang:       johnwang@g.harvard.edu
James O'Dwyer:   jodwyer@illinois.edu
Rafael D'Andrea: rafael.dandrea@stonybrook.edu


*Study summary
Here, we aim to quantify the signatures of niche differentiation in tree species in three tropical forest dynamic plots.
We demontrate a novel method that quantifies distict niches the strength of niche differentiation based on spatial 
clustering tendancies between pairs of species. We also characterize these clusters of species in terms of their 
association with soil nutrients. We show that there is strong spatial clustering of species across all three plots;
further, these clusters are also strongly associated with the spatial distribution of soil nutrients, indicating 
that the different groups correspond to different niche strategies for local soil conditions.


*Authors responsible for writing code

Rafael D'Andrea
Mihir Umarani
John Wang


*Description of the contents of the repository
Folders: 
1. codes:This folder contains all the codes written in R (R version 4.0.5)
 	1. clustering_functions_rann.R: This file contains necessary functions for the data analysis. This file should be 
           sourced into the rest of the R files to call functions.
	2. bci.R, laplanada.R and yasuni.R contain the data analysis (cluster analysis, soil nutrient analysis etc.) of the 
   	   three forest dynamic plot datasets respectively. bci.R also contains the recruitment analysis across censuses.

2. Final datasets:This folder contains 3 sub-folders containing the raw census data (either openly accessible or with permission 
		  from the PIs) and processed results for each of the forest dynamic plots.

		1. Raw files: i) bci_raw.rds- BCI census data compiled from dryad data repository:https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
			     ii) lap_census_data.rds- La Planada census data compiled from the repository:http://i2d.humboldt.org.co/ceiba/resource.do?r=planada_parcelapermanente_censo2&request_locale=ja
			     iii)  yasuni.fulldata.txt- ____________________

		2. Files named ___cluster_analysis.rds show the results of clustering analysis and have the following columns:
		   sp: species code; 
		   group:cluster to which the species belongs; 
		   algorithm: Name of the alogorithm to calculate modularity of the clustering network;
		   weighted:Whether the edges in the network were weighted;
		   number_of_groups: Total number of distinct clusters found;
		   d_cutoff (model parameter): Distance cutoff (in meters) used to identify the 'neighbor' trees
		   seed: seed values used for randomization of the pairwise distance matrix for species. 0 indicates no
			changes. 
		3. Files named ___kde_full.rds shows the results of kernel density estimation of distinct spatial cluster of 
		   species (basically shows the density of each cluster in every cell of a map) and have the following columns:
		   census:(if there are multiple censuses) 
		   d_cutoff (model parameter):Distance cutoff (in meters) used to identify the 'neighbor' trees
		   group: cluster group inferred from the clustering analysis
		   x and y: x and y coordinates of a cell (10m X 10m)
		   density: Estimated density of the individuals from a given species cluster from the kernel density estimation.
		   soiltype:Cluster from the census 1 (reference cluster) corresponding to the group.
		   fdp: Forest dynamic plot.

			
Coding software and packages:
All the code was executed in R, version (4.0.5).
Essential packages used (versions)-
RANN (2.6.1)
igraph (1.3.0)
sparr(2.2-15)
		   

There are 3 R code files under the sub-folder "codes" which 
contain the analyses for each of the forest dynamic plots. The file labelled
"clustering_functions_rann.R" contains supporting functions. Please refer to 
the annotations in the code files for further details on the analyses.

All the data files are in the folder "Final datasets".

Following are the raw input files for each of the FDPs:

BCI:
Census data:        "bci_raw.rds" 

La Planada:
Census data: "lap_census_data.rds"

Yasuni:
Census data: "yasuni.fulldata.txt"

