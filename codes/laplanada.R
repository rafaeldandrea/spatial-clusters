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

#Source doi:https://doi.org/10.15472/ekg6vs

#Download the DwC file. Unzip the folder. There are two text files with tables for tree occurrence and dbh data. PLace them in a working
#directory.
#Filter the occurrence data by census==2 (DatasetID=="Censo: 2"). Filter the dbh data for for only primary branches (CodTallo=="A").
#Combine both the files and label the new file as lap_census_data.rds (newly created files is saved as .rds and .csv in the repo).
#Note that the column names are changed from the original spanish labels. 

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
    "dbh"=DAP)%>%as_tibble()


#Soil nutrient data
#Source:https://doi.org/10.13012/B2IDB-6140727_V1 
#Download the excel file labelled, "laplanada.data.deposited". 
#Create an independent CSV file from the sheet called "Block Estimates" of the 20X20m krieged data. 
# Standardize the gx and gy values between 1 and 25 (from the range of 0-500).
#Label it as "lap_20X20_soil.csv"

nutrients=read.csv("lap_20x20_soil.csv")%>%as_tibble()%>%
          mutate(
          across(
            Al:pH, 
            function(x) scale(x)[, 1]
          )
        )

###################################################################################################################
#1. Spatial clustering analysis:

# Step 1:For each census, count the number of heterospecific trees in the neighborhood of each tree,
# where the neighborhood is a circle of radius = 10m, 20m and 30m. 
# Step 2: Build an adjacency matrix where each entry (Aij)  in a matrix is the number of neighbors
# of species i which belong to species j (significance is determined by the relative no. of conspecific neighbors)
# Step 3: Calculate maximum modularity of the matrix and determine the strength of clustering
# as well as no. of clusters. 
#Caution- lengthy execution time! One iteration of parameter combination takes ~5-6 minutes.

Lx=500
Ly=500

#Cutoffs to determine the adult individuals of each species

#Baldeck cutoffs
baldeck_cutoff = 
 lap %>%
  group_by(sp) %>%
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
        
      if(cutoff=='rf'){
        data = 
          lap %>%
          inner_join(rf_cutoff) %>%
          filter(dbh > rf) 
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
          census = thecensus,
          d_cutoff = d_cutoff,
          seed = seed,
          method=method,
          cutoff=cutoff
        )%>%
        return()
    },
    .options = furrr_options(seed = TRUE)
  ) %>%
  rename(sp = name)

saveRDS(fdp_analyzed, file = "lap_clustering_analysis.rds")


##################################################################################################################
#2. Kernel density estimation: Kernel smoothing of geographic areas of each of the inferred clusters.


#Need the clustering results for this-

#If it is not already created/loaded
cluster_data = readRDS('lap_clustering_analysis.rds')

data = 
  lap %>% 
  select(sp,gx, gy) %>%
  full_join(
      cluster_data %>%
      filter(seed == 0) %>%
      select(sp, d_cutoff, group)
  ) %>%drop_na()%>%
  select( d_cutoff, gx, gy, group)


kde_full = 
  expand_grid(census=1,
    d_cutoff = unique(data$d_cutoff)
  ) %>%
  future_pmap_dfr(
    .f = KDE,
    .data = data,
    .options = furrr_options(seed = NULL)
  )

kde_full=kde_full%>%
  inner_join(
    kde_full %>%
      group_by(d_cutoff, x, y) %>%
      slice_max(density, n = 1, with_ties = FALSE) %>%
      rename(soiltype = group) %>%
      ungroup() %>%
      select(-density))%>%
  mutate(soiltype=as.factor(soiltype))

saveRDS(kde_full,"lap_kde_full.RDS")


#Plot spatial clusters by census and d_cutoffs

lap_clustmap_all=readRDS('lap_kde_full.rds')%>%
  group_by(d_cutoff,census,x,y)%>%
  slice_max(density,with_ties=FALSE)%>%
  ungroup()%>%
  ggplot(aes(x,y,fill=soiltype))+
  geom_tile()+
  #facet_wrap(~census)+
  theme(aspect.ratio = 0.5,
        strip.text = element_text(size =15,face='bold'),
        legend.title= element_text(size = 10,face="bold"),
        axis.title= element_text(size = 15,face="bold")
  )+
  labs(fill='Species cluster')+
  scale_fill_manual(values=cbpalette[c(1,2,3,4,5)])+
  #theme(legend.position="top")+
  facet_grid(~d_cutoff)

lap_clustmap=readRDS('lap_kde_full.rds')%>%
  filter(d_cutoff==20)%>%
  group_by(x,y)%>%
  slice_max(density,with_ties=FALSE)%>%
  ungroup()%>%
  ggplot(aes(x,y,fill=soiltype))+
  geom_tile()+
  #facet_wrap(~census)+
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

cluster_data=readRDS("lap_clustering_analysis.rds")

kde_full=readRDS('lap_kde_full.rds')


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
  select(Al:H) %>% 
  pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)

df_scores = 
  pca@scores %>% 
  as_tibble() %>% 
  bind_cols(data)

k = 
  future_pmap_dfr(
    expand_grid(
      groups = 2:10, 
      seed = 0:50
    ),
    function(groups, seed){
      foo = 
        df_scores %>%
        select(Al:pH)
      if(seed > 0){
        set.seed(seed)
        foo = apply(foo, 2, sample)
      }
      wss = 
        kmeans(
          foo, 
          centers = groups, 
          nstart = 100
        )$tot.withinss
      return(tibble(groups = groups, seed = seed, wss = wss))
    } 
  )

k%<>%mutate(wss=-log(wss))

plot_gap = 
  k %>%
  filter(seed == 0) %>% 
  inner_join(
    k %>% 
      filter(seed > 0) %>% 
      group_by(groups) %>% 
      summarize(
        null = mean(wss), 
        .groups = 'drop')
  ) %>% 
  mutate(gap = null - wss) %>%
  ggplot(aes(groups, gap)) +
  geom_line() + 
  geom_point() +
  xlab('number of groups') +
  ylab('Gap statistic (log)')+
  theme(aspect.ratio = 1,
        axis.title=element_text(size=15,face='bold'))

k4 = 
  kmeans(
    df_scores %>% 
      select(Al:pH),
    centers = 4,
    nstart = 100
  )

df_scores = 
  df_scores %>%
  mutate(kmeans = k4$cluster)


lap_soilclus=df_scores%>%
  mutate(kmeans=as.factor(kmeans))%>%
  ggplot(aes(x,y,fill=kmeans))+
  geom_tile()+
  labs(fill=gsub('\\s','\n',"Soil-Nutrient Cluster"))+
  #scale_fill_brewer(palette='YlGnBu')+
  scale_fill_manual(values=cbpalette[1:4])+
  theme(aspect.ratio = 0.5,
        legend.key.size = unit(0.4, 'cm'),
        legend.title=element_text(size=8))


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
               pivot_longer(cols=Al:water,names_to='nutrient',values_to='value')%>%
               group_by(nutrient)%>%
               mutate(values=(value - min(value)) / (max(value) - min(value)))%>%
               ungroup())%>%
  group_by(group,nutrient)%>%
  summarize(
    cor = cor(density, values, use = 'complete.obs'),
    p.value = cor.test(density, values)$p.value,
    p.adj = p.adjust(p.value, method = 'hochberg'),
    significant = (p.adj < .05),
    cor_sig = ifelse(significant == TRUE, cor, NA)
  )%>%
  ungroup()

cordat%>%
  mutate(group=as.factor(group),
         nutrient=as.factor(nutrient))%>%
  ggplot(aes(nutrient,cor_sig,fill=group))+
  geom_bar(stat='identity')+
  facet_grid(~group)%>%
  scale_fill_manual(values=cbpalette)+
  theme(axis.text=element_text(size=12,face='bold'),
        legend.title=element_text(size=12,face='bold'))+
  xlab('Nutrients')+
  ylab("Pearson's correlation coefficient")

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


cluster_data=readRDS("lap_clustering_analysis.rds")

kde_full=readRDS('lap_kde_full.rds')

nutrients %<>%
  pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
  group_by(nutrient) %>%
  mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

dtf = 
  nutrients %>%
  select(-value) %>%
  pivot_wider(names_from = nutrient, values_from = standardized)%>%
  full_join(
    kde_full%>%
      filter(d_cutoff==20)%>%
      select(x,y,group)
  )

C5_res=C5.0(dtf%>%select(-group),
            dtf$group,
            rule=TRUE,
            control=C5.0Control(bands=10,winnow=TRUE),
            trials=1)


C5_model=train(dtf%>%select(-group),dtf$group,method='C5.0',metric='Kappa',tunelength=10,trControl=cv)

C5_model$results %>%
  as_tibble %>%
  bind_cols(parms[index, ]) %>%
  return()
