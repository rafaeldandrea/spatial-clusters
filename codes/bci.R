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


#Source important functions
source('C:/Users/mihir/Documents/spatial-clusters/codes/clustering_functions_rann.R')

#Plotting settings

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)

cbpalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#C5E772","#4F2CAA")

###############################################################################################################
#Load census data
###############################################################################################################

#Download the dryad repository https://doi.org/10.15146/5xcp-0d46
#Unzip the folder labelled "bci.tree.zip"

raw.files=list.files()[intersect(grep(".rdata",list.files()),grep("bci.tree",list.files()))]

#Load .rdata files
for(i in raw.files){
  load(i)
}

tree.files=ls()[grep("bci.tree",ls())]

mydat=lapply(tree.files,function(x){
  get(x)
})


all=do.call("rbind", mydat)%>%tibble()
bci_raw=all%>%mutate(census=rep(1:8,sapply(mydat, nrow)))
bci<-bci_raw%>%select(sp,gx,gy,dbh,census)%>%drop_na()

#Load BCI plant trait data
#Source for trait data: https://doi.org/10.6084/m9.figshare.c.3303654.v1
#Source for leaf stoichiometry data: https://doi.org/10.25573/data.23463254
#Convert the tab separated text files into csv files labeled as 'bci_trait_data.rds' and 'bci_st_tr.csv', respectively.

trait_data_raw = readRDS('bci_trait_data.rds')
trait_data_st=read.csv('bci_st_tr.csv')%>%as_tibble()

#BCI species information data
load('bci.spptable.rdata')

splist=bci%>%
        inner_join(
          bci %>%
            group_by(sp) %>%
            summarize(
              baldeck = quantile(dbh, .56), 
              .groups ='drop'
            ))%>%
        filter(dbh>baldeck)%>%
        select(sp)%>%
        unique()%>%
        pull()

bci.sp=bci.spptable%>%as_tibble()%>%
  filter(sp%in%splist)%>%
  select(sp,Genus,Species)%>%
  rename(genus=Genus,
         species=Species)%>%
  mutate(latin=paste0(genus," ",species))



#Data taken from a url "http://ctfs.si.edu/webatlas/datasets/bci/soilmaps/bci.block20.data.xls" 
#and save it in your working directory as "bci.block20.data.xls"

nutrients=read_xls("bci.block20.data.xls",sheet="means (20 X 20 m)")%>%
          as_tibble()


#Load soil water level data
#Source:https://doi.org/10.1038/s41597-019-0072-z
water =
  read.table('BCI_SWP_map_mid_dry_season_regular.txt',
    header = TRUE
  ) %>%
  as_tibble() %>%
  mutate(
    x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)],
    y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
  ) %>%
  replace_na(list(x = 10, y = 10)) %>%
  group_by(x, y) %>%
  summarize(water = mean(swp), .groups = 'drop')%>%
  mutate(water=scale(water)[,1])

nutrients=nutrients%>%inner_join(water)




##################################################################################################################
#1. Spatial clustering analysis:

# Step 1:For each census, count the number of heterospecific trees in the neighborhood of each tree,
# where the neighborhood is a circle of radius = 10m, 20m and 30m. 
# Step 2: Build an adjacency matrix where each entry (Aij)  in a matrix is the number of neighbors
# of species i which belong to species j (significance is determined by the relative no. of conspecific neighbors)
# Step 3: Calculate maximum modularity of the matrix and determine the strength of clustering
# as well as no. of clusters. 
#Caution- lengthy execution time! One iteration of parameter combination takes ~5-6 minutes.

Lx=1000
Ly=500



#Cutoffs to determine the adult individuals of each species

#Baldeck cutoffs
baldeck_cutoff = 
  bci %>%
  group_by(sp) %>%
  summarize(
    baldeck = quantile(dbh, .56), 
    .groups ='drop'
  )

cutoff_methods<-'baldeck'

thecensus<-bci%>%pull(census)%>%unique()

parameters=
  expand.grid(
    cutoff=cutoff_methods,
    thecensus=thecensus,
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
    thecensus,
    algorithm,
    d_cutoff,
    weighted,
    seed
    ){
      
      if(cutoff=='raw'){
        dat<-
          bci%>%
          filter(dbh>100)%>%
          mutate(fdp='bci')
      }
      if(cutoff=='baldeck'){
        dat = 
          bci %>%
          inner_join(baldeck_cutoff) %>%
          filter(dbh > baldeck) %>%
          mutate(fdp = 'bci')
      }
      if(cutoff=='rf'){
        dat = 
          bci %>%
          inner_join(rf_cutoff) %>%
          filter(dbh > rf) %>%
          mutate(fdp = 'bci')
      }
      
      
      data =
        dat %>%
        filter(census == thecensus)
      
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

saveRDS(fdp_analyzed, file = "clustering_analysis.rds")


#Make sure that the inferred clusters across censuses are consistent in terms of their species identity

ref=cluster_data%>%
  filter(seed==0)%>%
  select(sp,algorithm,census,d_cutoff,number_of_groups)%>%
  unique()%>%
  slice_min(number_of_groups,n=1,with_ties=FALSE)

ref=cluster_data%>%
  filter(census==8, d_cutoff==20)%>%unique()

x = 
  cluster_data %>% 
  filter(seed == 0) %>% 
  select(algorithm,census, d_cutoff, sp, group, number_of_groups)

x0 = 
  x %>%
  filter(census == 8, d_cutoff == 20)

x0=x0%>%
  inner_join(x0%>%
               count(group)%>%
               mutate(newgr=(length(n)+1)-rank(n))
  )%>%
  select(-group)%>%
  mutate(group=newgr)



par.cal=x%>%select(census,d_cutoff)%>%
  unique()%>%as.matrix()


res = NULL

for(i1 in 1:nrow(par.cal)){
  
  x1=x%>%filter(census==par.cal[i1,1],d_cutoff==par.cal[i1,2])
  res1=NULL
  
  for(i2 in 1:length(unique((x0$group)))){
    foo=x1%>%
      group_by(algorithm, census, d_cutoff,group)%>%
      summarize(
        int=length(intersect(sp,x0 %>% filter(group == i2) %>% pull(sp)))/(x0 %>% filter(group == i2) %>% pull(sp)%>%length()),
        .groups='drop')%>%
      mutate(reference_group=i2)
    
    res1=res1%>%bind_rows(foo)
  }
  res1<-find.consistency(res1)
  
  res=res%>%bind_rows(res1)
}

consistent_cluster_data = cluster_data%>%
  inner_join(res)%>%
  select(-group)%>%
  rename(group=newgroup)%>%select(-int)%>%select(-reference_group)

saveRDS(consistent_cluster_data, "bci_clustering_analysis_consistent.rds")



#2. Kernel density estimation: Kernel smoothing of geographic areas of each of the inferred clusters.
##################################################################################################################

#Need the clustering results for this-

#If it is not already created/loaded
consistent_cluster_data = readRDS('bci_clustering_analysis_consistent.rds')

data = 
  bci %>% 
  select(sp, census, gx, gy) %>%
  full_join(
    consistent_cluster_data %>%
      filter(seed == 0) %>%
      select(sp, census, d_cutoff, group)
  ) %>%drop_na()%>%
  select(census, d_cutoff, gx, gy, group)


kde_full = 
  expand_grid(
    census = unique(data$census),
    d_cutoff = unique(data$d_cutoff),
    Lx=Lx,
    Ly=Ly
  ) %>%
  future_pmap_dfr(
    .f = KDE,
    .data = data,
    .options = furrr_options(seed = NULL)
  )

kde_full=kde_full%>%
  inner_join(
    kde_full %>%
      group_by(census, d_cutoff, x, y) %>%
      slice_max(density, n = 1, with_ties = FALSE) %>%
      rename(soiltype = group) %>%
      ungroup() %>%
      select(-density))%>%
      mutate(soiltype=as.factor(soiltype))

saveRDS(kde_full,"bci_kde_full.RDS")


#Plot spatial clusters by census and d_cutoffs

bci_clustmap_all=readRDS('bci_kde_full.rds')%>%
  #filter(d_cutoff==20,census==7)%>%
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
  facet_grid(census~d_cutoff)

bci_clustmap_7=readRDS('bci_kde_full.rds')%>%
  filter(d_cutoff==20,census==7)%>%
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
  
  

#3. Perform the recruitment analysis to test the survival of the saplings, juveniles 
#and adult in and out of their inferred spatial niches (from the kde analysis).
# ## Definition of theta := P(recruit | match) / P(recruit | !match), where "match" means 
# finding an individual tree within its own spatial niche. 
##################################################################################################################

filter = dplyr::filter

cores = detectCores() 
plan(multisession, workers = cores)


#Data filters- 
#Species- Use only those species for which kde and clustering analysis was performed. 

bci_dat=bci_raw%>%
  drop_na()%>%
  select(census,sp,treeID,gx,gy,dbh)

kde_full=readRDS('bci_kde_full.rds')%>%
  select(census,d_cutoff,x,y,soiltype)%>%
  unique()

cluster_data=readRDS('bci_clustering_analysis_consistent.rds')%>%
  select(census,d_cutoff,sp,group)%>%
  unique()

#Baldeck cutoffs

lh_cutoff=bci_dat%>%
  group_by(sp)%>%
  summarize(juv=quantile(dbh,0.56),
            adult=quantile(dbh,0.88))%>%
  ungroup()

saplings=bci_dat %>%
  inner_join(lh_cutoff) %>%
  filter(dbh <= juv)%>%select(-c(juv,adult))

juveniles=bci_dat %>%
  inner_join(lh_cutoff) %>%
  filter(dbh > juv,dbh<=adult)%>%select(-c(juv,adult))

adults= bci_dat %>%
  inner_join(lh_cutoff) %>%
  filter(dbh > adult)%>%select(-c(juv,adult))


#When do the adult trees first appear (except for the ones from census 1)
recruits.saplings = 
  saplings %>%
  filter(census > 1) %>%
  group_by(treeID) %>%
  slice_min(census) %>%
  ungroup %>%
  mutate(recruit = TRUE)

recruits.juv = 
  juveniles %>%
  filter(census > 1) %>%
  group_by(treeID) %>%
  slice_min(census) %>%
  ungroup %>%
  mutate(recruit = TRUE)

recruits.adult = 
  adults %>%
  filter(census > 1) %>%
  group_by(treeID) %>%
  slice_min(census) %>%
  ungroup %>%
  mutate(recruit = TRUE)

saps =
  saplings %>% 
  left_join(
    recruits.saplings, 
    by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
  ) %>%
  replace_na(list(recruit = FALSE))%>%mutate(life='sapling')

juvs =
  juveniles %>% 
  left_join(
    recruits, 
    by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
  ) %>%
  replace_na(list(recruit = FALSE))%>%mutate(life='juvenile')

trees =
  adults %>% 
  left_join(
    recruits, 
    by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
  ) %>%
  replace_na(list(recruit = FALSE))%>%mutate(life='adult')

combined_data =
  bind_rows(trees,juvs,saps) %>%
  inner_join(
    consistent_cluster_data %>%
      filter(seed == 0), 
    by = c('census', 'sp')
  )%>%unique

reference =
  combined_data %>%
  inner_join(
    combined_data %>%
      select(census, d_cutoff, number_of_groups) %>%
      unique() %>%
      slice_max(number_of_groups, with_ties = FALSE)
  ) %>%
  count(group) %>%
  mutate(reference = ((max(rank(n))+1)-(rank(n)))) %>%
  full_join(
    combined_data 
    #    inner_join(
    #        combined_data %>%
    #         select(census, d_cutoff, number_of_groups) %>%
    #       unique() %>%
    #      slice_max(number_of_groups, with_ties = FALSE)
    # )
  ) %>%
  select(treeID, reference)



z_full = NULL
for(the_life in sort(unique(combined_data$life))){

for(thecensus in sort(unique(combined_data$census))){
  
  for(thed_cutoff in sort(unique(combined_data$d_cutoff))){
    
    z = NULL
    
    x =
      combined_data %>%
      filter(
        life==the_life,
        census == thecensus,
        d_cutoff == thed_cutoff
      ) %>%
      inner_join(
        reference,
        by = 'treeID'
      ) %>%
      count(
        census,
        d_cutoff,
        group,
        reference
      ) %>%
      mutate(group = factor(group)) %>%
      complete(group, nesting(census, d_cutoff, reference)) %>%
      replace_na(list(n = 0)) %>%
      arrange(desc(n))
    
    if(nrow(x) > 0){
      z %<>%
        bind_rows(
          tibble(
            life=the_life,
            census = thecensus,
            d_cutoff = thed_cutoff,
            group = x$group[1],
            reference = x$reference[1]
          )
        )
      
      for(i in 2:nrow(x)){
        
        g = x$group[i]
        r = x$reference[i]
        
        if(!g %in% z$group & !r %in% z$reference){
          z %<>%
            bind_rows(
              tibble(
                life=the_life,
                census = thecensus,
                d_cutoff = thed_cutoff,
                group = g,
                reference = r
              )
            )
          
          z_full %<>%
            bind_rows(z)
        }
      }
    }
  }
}
}

z_full %<>%
  unique()

combined_data_consistent =
  combined_data %>%
  mutate(group = factor(group)) %>%
  left_join(z_full) %>%
  mutate(group = reference) %>%
  select(-reference)

cluster_data =
  consistent_cluster_data %>%
  mutate(group = factor(group)) %>%
  inner_join(z_full) %>%
  mutate(group = reference) %>%
  select(-reference)



rec_df = 
  bind_rows(recruits.saplings%>%mutate(life='saplings'),
            recruits.juv%>%mutate(life='juveniles'),
            recruits.adult%>%mutate(life='adults')) %>%
  mutate(
    x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
    y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
  ) %>%
  inner_join(
    kde_full %>%
      select(
        labelling_census = census, 
        x, 
        y, 
        d_cutoff, 
        soiltype
      ) %>%
      unique(),
    by = c('x', 'y')
  ) %>%
  inner_join(
    cluster_data %>%
      rename(labelling_census = census)
  )

Pm_df = 
  kde_full %>% 
  select(x, y, d_cutoff, soiltype, census) %>% 
  unique() %>% 
  count(census, d_cutoff, soiltype) %>% 
  mutate(Pm = n / (Lx * Ly / 20 ^ 2)) %>% 
  ungroup() %>%
  rename(labelling_census = census)

res =
  rec_df %>%
  group_by(
    labelling_census,
    d_cutoff,
    group
  ) %>%
  summarize(
    recruits = n(),
    matches = sum(group == soiltype),
    .groups = 'drop'
  ) %>%
  left_join(
    Pm_df %>%
      rename(group = soiltype)
  ) %>%
  mutate(
    theta = ((1 - Pm) / Pm) / (recruits / matches - 1)
  )

res_summary = 
  res %>%
  group_by(labelling_census, d_cutoff) %>%
  mutate(weight = recruits / sum(recruits)) %>%
  summarize(
    theta_mean = weighted.mean(theta, weight),
    theta_se = sqrt(sum(weight * (theta - theta_mean) ^ 2)),
    .groups = 'drop'
  ) %>%
  group_by(d_cutoff) %>%
  summarize(
    theta_mean = mean(theta_mean),
    theta_se = sqrt(sum(theta_se ^ 2)) / n(),
    .groups = 'drop'
  )


lh_cutoff=bci_dat%>%
  group_by(sp)%>%
  summarize(juv=quantile(dbh,0.56),
            adult=quantile(dbh,0.88))%>%
  ungroup()

sap.rec=bci_dat%>%
  mutate(
    x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
    y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
  )%>%
  select(census,treeID,sp,x,y,dbh)%>%
  inner_join(lh_cutoff)%>%
  filter(dbh<juv)%>%
  select(-c(juv,adult))%>%
  filter(census>1)%>%
  group_by(treeID)%>%
  slice_min(census)%>%
  ungroup()


juv.rec=bci_dat%>%
  mutate(
    x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
    y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
  )%>%
  select(census,treeID,sp,x,y,dbh)%>%
  inner_join(lh_cutoff)%>%
  filter(dbh>=juv, dbh<adult)%>%
  select(-c(juv,adult))%>%
  filter(census>1)%>%
  group_by(treeID)%>%
  slice_min(census)%>%
  ungroup()

adult.rec=bci_dat%>%
  mutate(
    x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
    y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
  )%>%
  select(census,treeID,sp,x,y,dbh)%>%
  inner_join(lh_cutoff)%>%
  filter(dbh>=adult)%>%
  select(-c(juv,adult))%>%
  filter(census>1)%>%
  group_by(treeID)%>%
  slice_min(census)%>%
  ungroup()


pm.df= kde_full%>%
  filter(census==1)%>%
  select(x,y,census,d_cutoff,soiltype)%>%
  count(d_cutoff,soiltype)%>%
  mutate(pm=n/(Lx*Ly/400))

sap.dat=sap.rec%>%
  inner_join(
    kde_full%>%
      filter(census==1)%>%
      select(-census)
  )%>%
  inner_join(cluster_data)

sap.res=sap.dat%>%
  group_by(census,d_cutoff,group)%>%
  summarize(
    recruits=n(),
    matches=sum(group==soiltype),
    .groups = 'drop'
  )%>%
  mutate(group=as.factor(group))%>%
  left_join(
    pm.df%>%
      rename(group=soiltype)
  )%>%
  mutate(
    theta = ((1 - pm) / pm) / ((recruits / matches) - 1)
  )%>%
  mutate(life='sapling')

sap.res$theta[is.infinite(sap.res$theta)]=NA


sap.summary =
  sap.res %>%
  drop_na()%>%
  group_by(census, d_cutoff) %>%
  mutate(weight = recruits / sum(recruits)) %>%
  summarize(
    theta_mean = weighted.mean(theta, weight),
    theta_se = sqrt(sum(weight * (theta - theta_mean) ^ 2)),
    .groups = 'drop'
  ) %>%
  group_by(d_cutoff) %>%
  summarize(
    theta_mean = mean(theta_mean),
    theta_se = sqrt(sum(theta_se ^ 2)) / n(),
    .groups = 'drop'
  )



juv.dat=juv.rec%>%
  inner_join(
    kde_full%>%
      filter(census==1)%>%
      select(-census)
  )%>%
  inner_join(cluster_data)



juv.res=juv.dat%>%
  group_by(census,d_cutoff,group,sp)%>%
  summarize(
    recruits=n(),
    matches=sum(group==soiltype),
    .groups = 'drop'
  )%>%
  left_join(
    pm.df%>%
      rename(group=soiltype)
  )%>%
  mutate(
    theta = ((1 - pm) / pm) / (recruits / matches - 1)
  )%>%
  mutate(life='juvenile')

juv.res$theta[is.infinite(juv.res$theta)]=NA


juv.summary =
  juv.res %>%
  drop_na()%>%
  group_by(census, d_cutoff) %>%
  mutate(weight = recruits / sum(recruits)) %>%
  summarize(
    theta_mean = weighted.mean(theta, weight),
    theta_se = sqrt(sum(weight * (theta - theta_mean) ^ 2)),
    .groups = 'drop'
  ) %>%
  group_by(d_cutoff) %>%
  summarize(
    theta_mean = mean(theta_mean),
    theta_se = sqrt(sum(theta_se ^ 2)) / n(),
    .groups = 'drop'
  )


adult.dat=adult.rec%>%
  inner_join(
    kde_full%>%
      filter(census==1)%>%
      select(-census)
  )%>%
  inner_join(cluster_data)



adult.res=adult.dat%>%
  group_by(census,d_cutoff,group,sp)%>%
  summarize(
    recruits=n(),
    matches=sum(group==soiltype),
    .groups = 'drop'
  )%>%
  left_join(
    pm.df%>%
      rename(group=soiltype)
  )%>%
  mutate(
    theta = ((1 - pm) / pm) / (recruits / matches - 1)
  )%>%
  mutate(life='adult')

adult.res$theta[is.infinite(adult.res$theta)]=NA


adult.summary =
  adult.res %>%
  drop_na()%>%
  group_by(census, d_cutoff) %>%
  mutate(weight = recruits / sum(recruits)) %>%
  summarize(
    theta_mean = weighted.mean(theta, weight),
    theta_se = sqrt(sum(weight * (theta - theta_mean) ^ 2)),
    .groups = 'drop'
  ) %>%
  group_by(d_cutoff) %>%
  summarize(
    theta_mean = mean(theta_mean),
    theta_se = sqrt(sum(theta_se ^ 2)) / n(),
    .groups = 'drop'
  )


rec.summary=rbind(sap.summary%>%mutate(life='sapling'),
                  juv.summary%>%mutate(life="juvenile"),
                  adult.summary%>%mutate(life='adult'))%>%
  mutate(life=as.factor(life))

rec.summary%>%
  ggplot(aes(x=life,y=theta_mean,col=life))+
  geom_point()+
  geom_errorbar(aes(ymin=theta_mean-theta_se,
                    ymax=theta_mean+theta_se)
  )+
  ylab('Theta')+
  xlab('Life Stage')+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=12))+
  facet_wrap(~d_cutoff)

rec.dat=bind_rows(sap.res,
                  juv.res,
                  adult.res)%>%
  mutate(life=as.factor(life),
         census=as.factor(census),
         d_cutoff=as.factor(d_cutoff))



rec.dat%>%
  ggplot(aes(life,theta,col=life))+
  geom_violin(draw_quantiles = .5)+
  facet_wrap(~census)+
  labs(col)

plot_histograms =
  res %>%
  ggplot(aes(theta)) +
  geom_histogram() +
  facet_wrap(~ d_cutoff)


#4. Perform PCA and kmeans clustering to calculate covariation and spatial distribution
# of different soil nutrients.
##################################################################################################################
#load kde and clustering data if not already loaded

cluster_data=readRDS("bci_clustering_analysis_consistent.rds")

kde_full=readRDS('kde_full.rds')


data=kde_full%>%
      filter(census==7,d_cutoff==20)%>%
      group_by(x,y)%>%
      slice_max(density)%>%
      select(x,y,group)%>%
      inner_join(nutrients%>%
                   mutate(
                     across(
                       Al:water, 
                       function(x) scale(x)[, 1]
                     )
                   ))

set.seed(0)

pca = 
  nutrients %>% 
  select(Al:water) %>% 
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
        select(Al:water)
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
      select(Al:water),
    centers = 4,
    nstart = 100
  )

df_scores = 
  df_scores %>%
  mutate(kmeans = k4$cluster)


bci_soilclus=df_scores%>%
  mutate(kmeans=as.factor(kmeans))%>%
  ggplot(aes(x,y,fill=kmeans))+
  geom_tile()+
  labs(fill=gsub('\\s','\n',"Soil-Nutrient Cluster"))+
  #scale_fill_brewer(palette='YlGnBu')+
  scale_fill_manual(values=cbpalette[1:4])+
  theme(aspect.ratio = 0.5,
        legend.key.size = unit(0.4, 'cm'),
        legend.title=element_text(size=8))


bci_cor.kmeans=df_scores%>%
  pivot_longer(Al:water,names_to = 'nutrients',values_to = 'conc')%>%
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
        filter(census==7,d_cutoff==20)%>%
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


cluster_data=readRDS("bci_clustering_analysis_consistent.rds")

kde_full=readRDS('kde_full.rds')

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
      filter(census==7,
             d_cutoff==20)%>%
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


###########################################################################################

#Morpho/growth trait analysis
#Perform PCAs and kmeans analysis

#KDE functions:
KernelDensityEstimation =
  function(gx, gy, Lx = L, Ly = 500, quadrat_length = 20, ...){
    
    evalpoints =
      expand_grid(
        x = quadrat_length / 2 + seq(0, Lx - quadrat_length, by = quadrat_length),
        y = quadrat_length / 2 + seq(0, Ly - quadrat_length, by = quadrat_length)
      ) %>%
      arrange(y, x)
    
    dens =
      bivariate.density(
        ppp(gx, gy, xrange = c(0, Lx), yrange = c(0, Ly)),
        h0 = quadrat_length,
        xy =
          list(
            x = evalpoints$x,
            y = evalpoints$y
          )
      )
    
    evalpoints %>%
      mutate(
        density = as.numeric(t(dens$z$v))
      ) %>%
      return
  }

KDE =
  function(Census, Algorithm, Seed, D_cutoff, Group, .data){
    df =
      .data %>%
      filter(
        census == Census,
        algorithm == Algorithm,
        seed == Seed,
        d_cutoff == D_cutoff,
        group == Group
      ) %>%
      select(gx, gy) %>%
      unique()
    
    result =
      KernelDensityEstimation(
        gx = df$gx,
        gy = df$gy
      ) %>%
      mutate(
        census = Census,
        algorithm = Algorithm,
        seed = Seed,
        d_cutoff = D_cutoff,
        group = Group
      ) %>%
      return
  }



bci_dat=bci%>%
        inner_join(
        bci%>%
          group_by(sp)%>%
          summarize(baldeck=quantile(dbh,0.56))%>%
          ungroup())%>%
        filter(dbh>baldeck)

cluster_data=readRDS("bci_clustering_analysis_consistent.rds")

splist=unique(bci$sp)

size_traits =
  c(
    "DBH_AVG",
    "HEIGHT_AVG",
    "DIAM_AVG"
  )

leaf_traits =
  c(
    "LEAFAREA_AVD",
    "LEAFTHCK_AVD",
    "LMADISC_AVD",
    "LMALEAF_AVD",
    "LMALAM_AVD",
    "LDDISC_AVD",
    "LDMC_AVD",
    "LEAFAREA_AVI",
    "LEAFTHCK_AVI",
    "LMADISC_AVI",
    "LMALEAF_AVI",
    "LMALAM_AVI" ,
    "LDDISC_AVI" ,
    "LDMC_AVI"  ,
    "AVG_LAMTUF" ,
    "AVG_VEINTUF"
  )

seed_traits =
  c(
    "FRUIT_FRSH",
    "FRUIT_DRY" ,
    "DSPR_FRESH",
    "DSPR_DRY",
    "SEED_FRESH",
    "SEED_DRY"
  )

wood_traits =
  c(
    "SG60C_AVG",
    "SG100C_AVG"
  )

vital_traits =
  c(
    "RGR_10",
    "RGR_50",
    "RGR_100",
    "MORT_10",
    "MORT_100"
  )

traitlist =
  list(
    vital = vital_traits,
    leaf = leaf_traits,
    seed = seed_traits,
    wood = wood_traits,
    size = size_traits
  )


foo =
  trait_data_raw %>%
  mutate(sp = tolower(`SP$`)) %>%
  pivot_longer(-c(1:6, sp), names_to = 'trait') %>%
  filter(!is.na(value)) %>%
  filter(!str_detect(trait, '_N')) %>%
  filter(!str_detect(trait, 'N_')) %>%
  filter(!str_detect(trait, '_SE')) %>%
  filter(!str_detect(trait, 'SEM_')) %>%
  filter(value > 0) %>%
  mutate(logged_value = log(value))

normality =
  foo %>%
  group_by(trait) %>%
  summarize(
    normal = shapiro.test(value)$p.value > .05,
    lognormal = shapiro.test(logged_value)$p.value > .05,
    .groups = 'drop'
  )

trait_data =
  foo %>%
  left_join(normality, by = 'trait') %>%
  mutate(standardized = ifelse(normal, value, logged_value)) %>%
  group_by(trait) %>%
  mutate(standardized = scale(standardized)[, 1]) %>%
  ungroup %>%
  select(sp, trait, standardized)

trait_data_w=trait_data%>%
  filter(sp%in%splist)%>%
  pivot_wider(names_from=trait,values_from=standardized)


pca_data2=trait_data_w%>%
  select(-sp)%>%
  pca(nPCs=4,method='ppca',scale='none',center=FALSE)

trait_data_w=trait_data_w%>%
  mutate(PC1=pca_data2@scores[,1],
         PC2=pca_data2@scores[,2])



trait_data_w%<>%
  inner_join(
    pca_data%>%
      select(sp,trait_type,pc1)%>%
      pivot_wider(names_from = trait_type,values_from=pc1))%>%
  select(sp,vital,leaf,seed,wood,size,PC1,PC2)



kmeans_dat=bci_dat%>%
  inner_join(trait_data_w)%>%
  inner_join(
    cluster_data%>%
      filter(census==7)%>%
      select(sp,group))


abuns=kmeans_dat%>%
  group_by(sp)%>%
  summarize(abun=n())%>%
  mutate(abun=log(abun))

trdat=trait_data_w%>%
  inner_join(abuns)  

k = 
  future_pmap_dfr(
    expand_grid(
      groups = 2:10, 
      seed = 0:50
    ),
    function(groups, seed){
      foo = 
        trdat %>%
        select(PC1:PC2)
      wgt=trdat$abun
      
      if(seed > 0){
        set.seed(seed)
        foo = apply(foo, 2, sample)
      }
      wss = 
        sum(kmeansW(
          foo, 
          centers = groups, 
          nstart = 100,
          weight = wgt
        )$withinss)
      
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
  theme(aspect.ratio = 1)


#Optimal clusters=3

k3=kmeansW(trdat%>%select(PC1,PC2),centers=3,weight = trdat$abun,nstart=100)

trdat%<>%mutate(kmeans=k3$cluster)%>%mutate(kmeans=as.factor(kmeans))

PC_cluster=trdat%>%
  ggplot(aes(PC1,PC2,col=kmeans))+geom_point(aes(size=abun))+
  xlab('Trait PC1')+
  ylab('Trait PC2')+
  scale_color_manual(values=cbpalette[c(1,2,3)])+
  theme(aspect.ratio=1,
        axis.title=element_text(size=12,face='bold'),
        legend.title=element_text(size=12,face='bold'))


trdat%>%
  pivot_longer(cols = vital:size, names_to='Trait',values_to='value')%>%
  ggplot(aes(kmeans,value,fill=Trait))+
  geom_bar(stat='identity', position=position_dodge())+
  scale_fill_manual(values=cbpalette)+
  theme(
    axis.title=element_text(size=12,face='bold'),
    legend.text = element_text(size=12),
    legend.title=element_text(size=12,face='bold')
  )+
  ylab('Values')


kmeans_dat%<>%
  inner_join(
    trdat%>%
      select(sp,kmeans)
  )

plot_kmeans=kmeans_dat%>%
  select(gx,gy,sp,kmeans,group,vital,leaf,seed,wood,size)%>%
  mutate(x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
         y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)])%>%
  group_by(x,y,kmeans)%>%
  summarize(n=n())%>%
  ungroup()%>%
  pivot_wider(names_from = kmeans,values_from=n)

names(plot_kmeans)=c('x','y','k1','k2','k3')



plot_kmeans%>%mutate(total=k1+k2+k3)%>%
  mutate(k1=k1/total,
         k2=k2/total,
         k3=k3/total)%>%
  pivot_longer(k1:k3,names_to='kmeans',values_to='density')%>%
  group_by(x,y)%>%
  slice_max(density,with_ties=FALSE)%>%
  ungroup()%>%
  mutate(kmeans=as.factor(kmeans))%>%
  ggplot(aes(x,y,fill=kmeans))+geom_tile()


plot_kmeans%>%mutate(k1=cut(k1,c(0,quantile(k1,0.3),quantile(k1,0.6),quantile(k1,1))),
                     k2=cut(k2,c(0,quantile(k2,0.3),quantile(k2,0.6),quantile(k2,1))),
                     k3=cut(k3,c(0,quantile(k3,0.3),quantile(k3,0.6),quantile(k3,1))))%>%
  ggplot(aes(x,y,fill=k3))+
  geom_tile()+
  theme(aspect.ratio = 0.5)



plot_grid(plot_gap,plot.soiltypes,labels=c("(A)","(B)"),nrow=2)

#Plot kmeans clusters

kde_kmeans<-kmeans_dat %>%
  select(gx,gy,kmeans)%>%
  group_by(kmeans) %>% 
  summarize(
    density = 
      KernelDensityEstimation(gx = gx, gy = gy, Lx = 1000),
    .groups = 'drop'
  )

kde_kmeans=bind_cols(
  kde_kmeans$kmeans,
  kde_kmeans$density
)
names(kde_kmeans)=c('kmeans','x','y','density')

kde_kmeans%>%
  group_by(x,y)%>%
  slice_max(density)%>%
  ungroup()%>%
  ggplot(aes(x,y,fill=kmeans))+
  geom_tile()+
  theme(aspect.ratio=0.5)

kde_kmeans%>%
  ggplot(aes(x,y,fill=density))+
  geom_tile()+
  facet_grid(rows=vars(kmeans))+
  scale_fill_gradientn(colours = terrain.colors(10))

trdat%>%
  select(c(sp,vital:size,kmeans))%>%
  pivot_longer(vital:size,names_to = 'trait',values_to = 'Concentration')%>%
  ggplot(aes(kmeans,Concentration,fill=trait))+
  geom_col(position='dodge')

trdat%>%
  select(c(sp,vital:size,kmeans))%>%
  pivot_longer(vital:size,names_to = 'trait',values_to = 'Concentration')%>%
  ggplot(aes(kmeans,Concentration))+geom_boxplot()+facet_wrap(~trait,scales='free')


c_table=kmeans_dat%>%
  select(group,kmeans)%>%
  count(group,kmeans)%>%
  pivot_wider(names_from = kmeans,values_from = n)

c_table=as.matrix(c_table[,2:4])
c_table[4,2]=0
rcompanion::cramerV(c_table)



#Leaf trait analysis
#Perform PCAs and kmeans analysis

############################################################################################
trait_data_st$SPECIES.=tolower(trait_data_st$SPECIES.)

foo=trait_data_st%>%
  filter(LIGHT.=='SUN')%>%
  select(-LIGHT.)%>%
  rename(sp=SPECIES.)%>%
  pivot_longer(AL:N,names_to = 'nutrient',values_to = 'conc')%>%
  group_by(sp,nutrient)%>%
  summarize(conc=mean(conc))%>%
  ungroup()

foo=foo%>%
  inner_join(
    foo%>%
      group_by(nutrient)%>%
      summarize(q05=quantile(conc,0.05,na.rm=TRUE),
                q95=quantile(conc,0.95,na.rm=TRUE))%>%
      ungroup()
  )%>%
  filter(conc<q95)%>%
  select(sp,nutrient,conc)%>%
  pivot_wider(names_from = nutrient,values_from = conc)%>%
  mutate(CN=C/N,
         NP=N/P)



bci_dat=bci%>%
  inner_join(
  bci%>%
    group_by(sp)%>%
    summarize(baldeck=quantile(dbh,0.56))%>%
    ungroup()
)%>%
  filter(dbh>baldeck)

cluster_data=readRDS("bci_clustering_analysis_consisten.rds")

kde_full=readRDS('kde_full.rds')

trdat=bci%>%select(gx,gy,sp)%>%
  inner_join(foo)%>%
  inner_join(cluster_data%>%
               filter(census==7)%>%
               unique()%>%
               select(sp,group))


trdat%>%
  pivot_longer(AL:NP,names_to='nutrient',values_to='value')%>%
  ggplot(aes(group,value))+geom_boxplot()+facet_wrap(~nutrient,scales='free')

trdat_long=trdat%>%
  select(c(AL:NP,group))%>%
  unique()%>%         
  pivot_longer(AL:NP,names_to = 'trait',values_to='conc')

trdat_long%>%
  drop_na()%>%
  group_by(trait)%>%
  mutate(conc=scale(conc)[,1])%>%
  ggplot(aes(group,conc))+
  geom_violin()+
  stat_summary(fun=mean, geom="point", size=2, color="red")+
  facet_wrap(vars(trait),scale='free')


trdat_long%>%
  group_by(group,trait)%>%
  summarize(mean=mean(conc,na.rm=T),
            sd=sd(conc,na.rm=T))%>%
  ungroup()%>%
  mutate(trait=as.factor(trait))%>%
  ggplot(aes())

cmpr=list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))

ggboxplot(trdat_long, x = "group", y = "conc",
          color = "group")+
  facet_wrap(~trait,scales='free')+
  stat_compare_means(comparisons = cmpr, tip.length=0,                                                                               label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")),
                     vjusts=-0.1) +
  theme_minimal() + 
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40")) + 
  # Don't include a legend
  guides(colour=FALSE)



splist=unique(bci$sp)


trait_data_w=foo%>%
  filter(sp%in%splist)%>%
  drop_na()


pca_data=trait_data_w%>%
  select(-sp)%>%
  pcaMethods::pca(nPCs=4,method='ppca',scale='none',center=FALSE)

trait_data_w=trait_data_w%>%
  mutate(PC1=pca_data@scores[,1],
         PC2=pca_data@scores[,2])


kmeans_dat=bci%>%
  inner_join(trait_data_w)%>%
  inner_join(
    cluster_data%>%
      filter(census==7)%>%
      select(sp,group)%>%
      unique())



abuns=kmeans_dat%>%
  group_by(sp)%>%
  summarize(abun=n())%>%
  mutate(abun=log(abun))

trdat=trait_data_w%>%
  inner_join(abuns)

k = 
  future_pmap_dfr(
    expand_grid(
      groups = 2:10, 
      seed = 0:100
    ),
    function(groups, seed){
      foo = 
        trdat %>%
        select(AL:ZN)
      wgt=trdat$abun
      
      if(seed > 0){
        set.seed(seed)
        foo = apply(foo, 2, sample)
      }
      wss = 
        sum(kmeansW(
          foo, 
          centers = groups, 
          nstart = 100,
          weight = wgt
        )$withinss)
      
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
  theme(aspect.ratio = 1)


#Optimal clusters=3

kdash=kmeansW(trdat%>%select(AL:N),centers=3,weight = trdat$abun,nstart=100)

trdat%<>%mutate(kmeans=kdash$cluster)%>%mutate(kmeans=as.factor(kmeans))

PC_cluster=trdat%>%
  ggplot(aes(PC1,PC2,col=kmeans))+geom_point(aes(size=abun))+
  xlab('Trait PC1')+
  ylab('Trait PC2')


kmeans_dat%<>%
  inner_join(
    trdat%>%
      select(sp,kmeans)
  )

plot_kmeans=kmeans_dat%>%
  select(c(gx,gy,sp,kmeans,AL:ZN))%>%
  mutate(x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
         y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)])%>%
  group_by(x,y,kmeans)%>%
  summarize(n=n())%>%
  ungroup()%>%
  pivot_wider(names_from = kmeans,values_from=n)

names(plot_kmeans)=c('x','y','k1','k2','k3')


plot_kmeans%>%mutate(total=k1+k2+k3)%>%
  mutate(k1=k1/total,
         k2=k2/total,
         k3=k3/total)%>%
  pivot_longer(k1:k3,names_to='kmeans',values_to='density')%>%
  group_by(x,y)%>%
  slice_max(density,with_ties=FALSE)%>%
  ungroup()%>%
  mutate(kmeans=as.factor(kmeans))%>%
  ggplot(aes(x,y,fill=kmeans))+geom_tile()




kde_kmeans<-kmeans_dat %>%
  select(gx,gy,kmeans)%>%
  group_by(kmeans) %>% 
  summarize(
    density = 
      KernelDensityEstimation(gx = gx, gy = gy, Lx = 1000),
    .groups = 'drop'
  )

kde_kmeans=bind_cols(
  kde_kmeans$kmeans,
  kde_kmeans$density
)
names(kde_kmeans)=c('kmeans','x','y','density')

kde_kmeans%>%
  group_by(x,y)%>%
  slice_max(density)%>%
  ungroup()%>%
  ggplot(aes(x,y,fill=kmeans))+
  geom_tile()+
  theme(aspect.ratio=0.5)

kde_kmeans%>%
  ggplot(aes(x,y,fill=density))+
  geom_tile()+
  facet_grid(rows=vars(kmeans))+
  scale_fill_gradientn(colours = terrain.colors(10))

trdat%>%
  select(c(sp,vital:size,kmeans))%>%
  pivot_longer(vital:size,names_to = 'trait',values_to = 'Concentration')%>%
  ggplot(aes(kmeans,Concentration,fill=trait))+
  geom_col(position='dodge')

trdat%>%
  select(c(sp,AL:ZN,kmeans))%>%
  pivot_longer(AL:ZN,names_to = 'trait',values_to = 'Concentration')%>%
  ggplot(aes(kmeans,Concentration))+geom_boxplot()+facet_wrap(~trait,scales='free')


c_table=kmeans_dat%>%
  select(group,kmeans)%>%
  count(group,kmeans)%>%
  pivot_wider(names_from = kmeans,values_from = n)

c_table=as.matrix(c_table[,2:4])
c_table[4,2]=0
rcompanion::cramerV(c_table)



