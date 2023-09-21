## DO NOT SKIP THIS!

## Instructions for downloading required datasets

###############################################################################################################
# Download census data
###############################################################################################################

#Source doi:https://doi.org/10.15472/ekg6vs

#Download the DwC file. Unzip the folder. There are two text files with tables for tree occurrence and dbh data. 
#Place them in a working directory.
#Filter the occurrence data by census==2 (DatasetID=="Censo: 2"). Filter the dbh data for for only primary branches (CodTallo=="A").
#Combine both the files and label the new file as lap_census_data.rds (newly created files is saved as .rds and .csv in the repo).
#Note that the column names are changed from the original spanish labels. 


###############################################################################################################
# Download soil nutrient data
###############################################################################################################

#Source:https://doi.org/10.13012/B2IDB-6140727_V1 
#Download the excel file labelled, "laplanada.data.deposited". 
#Create an independent CSV file from the sheet labelled "Block Estimates" of the 20X20m krieged data. 

#Label it as "lap_20X20_soil.csv"

###############################################################################################################
# Download custom functions
###############################################################################################################

# Navigate to https://github.com/rafaeldandrea/spatial-clusters/blob/main/codes/clustering_functions_rann_lap.R

# Download the raw file to the working directory.
###############################################################################################################

## Required libraries
library(tidyverse)
library(openxlsx)
library(magrittr)
library(furrr)
library(readxl)
library(parallelDist) ## for function parDist()
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
## a graph then find communities in the graph
library(RANN)  ## for neighbor-finding function nn2()
library(FactoClass)
library(C50)
library(caret)
library(sparr) # for function bivariate.density() in KDE()
library(pcaMethods)
library(rcompanion)
library(ggpubr)

filter = dplyr::filter

source('clustering_functions_rann_lap.R')

Lx=500
Ly=500

#Plotting settings


theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)

cbpalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#C5E772","#4F2CAA")



###############################################################################################################
#Load census data

df=read.table('occurrence.txt',sep='\t',header=TRUE)

df1=df%>%select(id,recordNumber,"verbatimLatitude","verbatimLongitude","scientificNameID")

trdat=read.table('measurementorfact.txt',sep='\t',header=TRUE)

trdat1=trdat%>%select(id,measurementType,measurementValue)%>%
  pivot_wider(names_from=measurementType,values_from = measurementValue)%>%
  filter(CodTallo=="A")%>%
  select(id,DAP)

lap=df1%>%
    inner_join(trdat1)%>%
    select(-id)%>%
    rename("treeID"=recordNumber,
    "gx"=verbatimLongitude,
    "gy"=verbatimLatitude,
    "sp"=scientificNameID,
    "dbh"=DAP)%>%
    mutate(dbh=as.numeric(dbh),
           sp=tolower(sp))%>%
    as_tibble()


# Read Soil nutrient data


nutrients=read.csv("lap_20x20_soil.csv")%>%as_tibble()%>%
          mutate(gx=gx+10,
                 gy=gy+10)%>%
          rename(x=gx,
                 y=gy)%>%
          mutate(
          across(
            Al:pH, 
            function(x) scale(x)[, 1]
          )
        )

###################################################################################################################
#1. Spatial clustering analysis:

# Step 1: Count the number of heterospecific trees in the neighborhood of each tree,
# where the neighborhood is a circle of radius = 10m, 20m and 30m. 
# Step 2: Build an adjacency matrix where each entry (Aij)  in a matrix is the number of neighbors
# of species i which belong to species j (significance is determined by the relative no. of conspecific neighbors)
# Step 3: Calculate maximum modularity of the matrix and determine the strength of clustering
# as well as no. of clusters. 



#Cutoffs to determine the adult individuals of each species

#Baldeck cutoffs
baldeck_cutoff = 
 lap %>%
  group_by(sp)%>%
  summarize(
    baldeck = quantile(dbh, .56), 
    .groups ='drop'
  )

cutoff_methods<-'baldeck'


parameters=
  expand.grid(
    cutoff=cutoff_methods,
    algorithm='louvain',
    d_cutoff = c(10,20,30),
    weighted=TRUE,
    seed=0:100
  )


###WARNING! Lengthy code! Each parameter combination requires 2-3 minutes. There
# are 300 total parameter combinations! Reduce the no. of combinations above (no.of seeds)
#if you just want to test the code.

fdp_analyzed = 
  parameters %>%
  future_pmap_dfr(
    .f = function(
    cutoff,
    algorithm,
    d_cutoff,
    weighted,
    seed
    ){
      
      if(cutoff=='raw'){
        data=
          lap%>%
          filter(dbh>100)
      }
      if(cutoff=='baldeck'){
        data = 
          lap %>%
          inner_join(baldeck_cutoff) %>%
          filter(dbh > baldeck)
      }
        
      if(seed > 0){
        data %<>%
          mutate(sp = sample(sp))
      }
      
      result = 
        data %>%
        adjacency_matrix(
          d_cutoff = d_cutoff, 
          Lx = Lx, 
          Ly = Ly
        ) %>%
        cluster_analysis(
          algorithm = algorithm,
          weighted = weighted
        ) %>%
        pluck('result') %>%
        mutate(
          d_cutoff = d_cutoff,
          seed = seed,
          cutoff=cutoff
        )%>%
        return()
    },
    .options = furrr_options(seed = TRUE)
  ) %>%
  rename(sp = name)

saveRDS(fdp_analyzed, file = "lap_clustering_analysis.rds")

write.csv(fdp_analyzed, file = "lap_clustering_analysis.csv")

##################################################################################################################
#2. Kernel density estimation: Kernel smoothing of geographic areas of each of the inferred clusters.

#Need the clustering results for this-

#If it is not already created/loaded
cluster_data = read.csv('lap_clustering_analysis.csv')%>%
                mutate(sp=tolower(sp))%>%
                as_tibble()

data = 
  lap %>% 
  select(sp,gx, gy) %>%
  full_join(
      cluster_data %>%
      filter(seed == 0) %>%
      select(sp, d_cutoff, group),
      relationship = "many-to-many")%>%
  drop_na()%>%
  select( d_cutoff, gx, gy, group)%>%
  mutate(census=1)


#Ignore the warnings from the KDE function

kde_full = 
  expand_grid(
    census=unique(data$census),
    d_cutoff = unique(data$d_cutoff)
  ) %>%
  future_pmap_dfr(
    .f = KDE,
    .data = data,
    .options = furrr_options(seed = NULL)
  )


saveRDS(kde_full,"lap_kde_full.rds")

write.csv(kde_full,"lap_kde_full.csv")


#Plot spatial clusters by census and d_cutoffs


lap_clustmap_all=read.csv('lap_kde_full.csv')%>%
  group_by(d_cutoff,census,x,y)%>%
  slice_max(density,with_ties=FALSE)%>%
  ungroup()%>%
  mutate(soiltype=as.factor(soiltype))%>%
  ggplot(aes(x,y,fill=soiltype))+
  geom_tile()+
  theme(aspect.ratio = 1,
        strip.text = element_text(size =15,face='bold'),
        legend.title= element_text(size = 10,face="bold"),
        axis.title= element_text(size = 15,face="bold")
  )+
  labs(fill='Species cluster')+
  scale_fill_manual(values=cbpalette[c(1,2,3,4,5)])+
  facet_grid(~d_cutoff)

lap_clustmap=read.csv('lap_kde_full.csv')%>%
  filter(d_cutoff==20)%>%
  group_by(x,y)%>%
  slice_max(density,with_ties=FALSE)%>%
  ungroup()%>%
  mutate(soiltype=as.factor(soiltype))%>%
  ggplot(aes(x,y,fill=soiltype))+
  geom_tile()+
  theme(aspect.ratio = 0.5)+
  labs(fill='Species cluster')+
  scale_fill_manual(values=cbpalette[c(1,2,3,4,5)])+
  theme(legend.position="top",
        legend.title= element_text(size = 15,face="bold"),
        axis.title= element_text(size = 15,face="bold"))


##################################################################################################################

#Perform PCA and kmeans clustering to calculate covariation and spatial distribution
# of different soil nutrients.


#load kde and clustering data if not already loaded
cluster_data=read.csv("lap_clustering_analysis.csv")%>%as_tibble()

kde_full=read.csv('lap_kde_full.csv')%>%as_tibble()


data=kde_full%>%
  filter(d_cutoff==20)%>%
  group_by(x,y)%>%
  slice_max(density)%>%
  ungroup()%>%
  select(x,y,group)%>%
  inner_join(nutrients)

set.seed(0)

pca = 
  nutrients %>% 
  select(Al:pH) %>% 
  pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)

df_scores = 
  pca@scores %>% 
  as_tibble() %>% 
  bind_cols(data)

#K-means clustering 

#Warning: Lengthy code 
k =
  expand_grid(
    groups = 2:10,
    seed = 0:10
  ) %>%
  future_pmap_dfr(
    .f = 
      function(.data, groups, seed) {
        foo = .data
        if (seed > 0) {
          set.seed(seed)
          foo = apply(foo, 2, sample)
        }
        wss = sum(kmeansW(foo, centers = groups, nstart = 100)$withinss)
        return(
          tibble(
            groups = groups,
            seed = seed,
            wss = wss,
            lwss = log(wss)
          )
        )
      },
    .data = df_scores %>% select(PC1:PC3),
    .options = furrr_options(seed = NULL)
  )



k = 
  future_pmap_dfr(
    expand_grid(
      groups = 2:10, 
      seed = 0:100
    ),
    function(groups, seed){
      foo = 
        df_scores %>%
        select(PC1:PC3)
      if(seed > 0){
        set.seed(seed)
        foo = apply(foo, 2, sample)
      }
      wss = sum(
        kmeansW(
          foo, 
          centers = groups, 
          nstart = 100
        )$withinss)
      return(tibble(groups = groups, seed = seed, wss = wss))
    },
    .options=furrr_options(seed=NULL) 
  )

k%<>%mutate(lwss=log(wss))

gap =
  k %>%
  filter(seed == 0) %>%
  inner_join(
    k %>%
      filter(seed > 0) %>%
      group_by(groups) %>%
      summarize(
        sd = sd(lwss),
        null = mean(lwss),
        .groups = 'drop'
      ),
    by = 'groups'
  ) %>%
  mutate(gap = null - lwss)

plot_gap =
  gap %>%
  ggplot(aes(groups, gap)) +
  geom_line() +
  geom_point() +
  xlab('number of groups') +
  ylab('Gap statistic (log)') +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 15, face = 'bold')
  )

gap_tib = 
  gap |>
  mutate(
    gap_diff = c(NA, diff(gap)) - c(NA, sd[-1]),
    test = gap_diff < 0
  )

number_of_clusters = with(gap_tib, which(test == TRUE)[1])


#3 clusters!

k3 = 
  kmeans(
    df_scores %>% 
      select(PC1:PC3),
    centers = 3,
    nstart = 100
  )

df_scores = 
  df_scores %>%
  mutate(kmeans = k3$cluster)


lap_soilclus=df_scores%>%
  mutate(kmeans=as.factor(kmeans))%>%
  ggplot(aes(x,y,fill=kmeans))+
  geom_tile()+
  labs(fill=gsub('\\s','\n',"Soil-Nutrient Cluster"))+
  scale_fill_manual(values=cbpalette[1:4])+
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.7, 'cm'),
        legend.title=element_text(size=10))


lap_cor.kmeans=df_scores%>%
  pivot_longer(Al:pH,names_to = 'nutrients',values_to = 'conc')%>%
  group_by(kmeans,nutrients)%>%
  summarize(mean=mean(conc))%>%
  ungroup()%>%
  mutate(kmeans=as.factor(kmeans))%>%
  ggplot(aes(nutrients,mean,fill=kmeans))+
  geom_col(position='dodge')+
  facet_grid(cols=vars(kmeans))+
  theme(axis.text.x = element_text(angle = 90,size=8),
        aspect.ratio=0.75)+
  scale_fill_manual(values=cbpalette[c(1,2,3,4)])

#Plot the correlations of spatial clusters and soil nutrients
cordat=kde_full%>%
  filter(d_cutoff==20)%>%
  inner_join(nutrients%>%
               pivot_longer(cols=Al:pH,names_to='nutrient',values_to='value')%>%
               group_by(nutrient)%>%
               mutate(values=(value - min(value)) / (max(value) - min(value)))%>%
               ungroup(),
             relationship =
               "many-to-many")%>%
  group_by(group,nutrient)%>%
  summarize(
    cor = cor(density, values, use = 'complete.obs'),
    p.value = cor.test(density, values)$p.value,
    p.adj = p.adjust(p.value, method = 'hochberg'),
    significant = (p.adj < .05),
    cor_sig = ifelse(significant == TRUE, cor, NA)
  )%>%
  ungroup()

#Ignore the warning message with this code. 
cordat%>%
  mutate(group=as.factor(group))%>%
  ggplot(aes(nutrient,cor_sig,fill=group))+
  geom_bar(stat='identity',position=position_dodge())+
  facet_wrap(vars(group))+
  theme(axis.text=element_text(size=10,angle=90),
        axis.title=element_text(size=12,face='bold'),
        legend.title=element_text(size=12,face='bold'))+
  xlab('Nutrients')+
  ylab("Pearson's correlation coefficient")+
  labs(fill='Spatial cluster')+
  scale_fill_manual(values=cbpalette[c(1,2,3,6)])

#############################################################################################

#Perform C5.0 analysis to find associations between spatial clusters and soil nutrient variables


cluster_data=read.csv("lap_clustering_analysis.csv")%>%as_tibble()

kde_full=read.csv('lap_kde_full.csv')%>%as_tibble()

nutrients1=nutrients %>%
  pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
  group_by(nutrient) %>%
  mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

dtf = 
  nutrients1 %>%
  select(-value) %>%
  pivot_wider(names_from = nutrient, values_from = standardized)%>%
  full_join(
    kde_full%>%
      filter(d_cutoff==20)%>%
      select(x,y,group)
  )

#Start the analysis
#Vary the parameter 'mincases' to determine the minimum no. of samples per split
#of the dataset.
#Record the results in terms of kappa

C5_res=C5.0(dtf%>%select(Al:pH),
            dtf$group,
            rule=TRUE,
            control=C5.0Control(bands=10,winnow=TRUE,minCases=10),
            trials=1)

cv <- trainControl(method = "repeatedcv", number =10, repeats =10)


C5_model=train(dtf%>%select(Al:water),dtf$group,method='C5.0',metric='Kappa',tunelength=10,trControl=cv)

