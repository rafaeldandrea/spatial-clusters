library(tidyverse)
library(magrittr)
library(furrr)
library(parallel)
library(readxl)
library(parallelDist) ## for function parDist()
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
## a graph then find communities in the graph
library(tidyverse)
library(furrr)  ## for parallel computing
library(parallel)  
library(magrittr)
library(RANN)  ## for neighbor-finding function nn2()


theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)


#Load BCI census data

#label the folder "bci_dryad and UNZIP it"
dir<-getwd()
setwd(paste0(dir,"/","bci_dryad"))

raw.files<-list.files()

#Load data from DRYAD datasets: Census 7
mydat<- lapply(raw.files, function(x) {
  load(file = x)
  get(ls()[ls()!= "filename"])
})

setwd(dir)

all <- do.call("rbind", mydat)%>%tibble()
bci<-all%>%mutate(census=rep(1:8,sapply(mydat, nrow)))
bci<-bci%>%select(sp,gx,gy,dbh,census)%>%drop_na()


#Load BCI soil nutrient data
nutrients=readRDS('bci_nutrient_data.rds')


#Load BCI plant trait data
trait_data_raw = readRDS('bci_trait_data.rds')
trait_data_st=read.csv('bci_st_tr.csv')%>%as_tibble()

#Load soil water level data
water =
  read.table(
    url(
      'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
    ),
    header = TRUE
  ) %>%
  as_tibble() %>%
  mutate(
    x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)],
    y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
  ) %>%
  replace_na(list(x = 10, y = 10)) %>%
  group_by(x, y) %>%
  summarize(water = mean(swp), .groups = 'drop')



#List of analyses performed

#1. Spatial clustering analysis:
# Step 1:For each census, count the number of heterospecific trees in the neighborhood of each tree,
# where the neighborhood is a circle of radius = 10m, 20m and 30m. 
# Step 2: Build an adjacency matrix where each entry (Aij)  in a matrix is the number of neighbors
# of species i which belong to species j (significance is determined by the relative no. of conspecific neighbors)
# Step 3: Calculate maximum modularity of the matrix and determine the strength of clustering
# as well as no. of clusters. 
do.clustering.analysis = 1

#2. Kernel density estimation: Kernel smoothing of geographic areas of each of the inferred clusters.
do.kde.analysis = 0

#3. Perform the recruitment analysis to test the survival of the saplings, juveniles 
#and adult in and out of their inferred spatial niches (from the kde analysis).
# ## Definition of theta := P(recruit | match) / P(recruit | !match), where "match" means 
# finding an individual tree within its own spatial niche. 

do.recruitment.analysis = 0

#4. Perform PCA and kmeans clustering to calculate covariation and spatial distribution
# of different soil nutrients.
do.nutrient.analysis = 0



do.C5.analysis = 0


do.trait.analysis = 0



L = ifelse(fdp == 'bci', 1000, 500)

#
if(do.clustering.analysis){
  
  if(do.data){
    
    cores = max(4,  detectCores() - 2)
    plan(multisession, workers = cores)
    
    save_date = gsub('-', '', Sys.Date())
    
    source('clustering_functions_rann.R')
    
    if(fdp == 'bci'){
      Lx = 1000
      filename = paste0('bci_clustering_analysis_500.rds')
      
      bci = 
        readRDS(
          url(
            'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bcifullcensus.rds?raw=true'
          )
        )
      
      
      #Robin Foster cutoff
      rf_cut=function(x){
        x1<-x[order(x,decreasing=TRUE)][1:6]
        return(mean(x1)/2)
      }
      
      rf_cutoff=
        bci%>%
        group_by(sp)%>%
        summarize(
          rf=rf_cut(dbh),
          .groups='drop'
        )
      
      #Baldeck cutoffs
      baldeck_cutoff = 
        bci %>%
        group_by(sp) %>%
        summarize(
          baldeck = quantile(dbh, .56), 
          .groups ='drop'
        )
    }
    
    #thecensus= dat %>%
    # pull(census) %>%
    #unique()
 
    
    #parameters = 
    # expand_grid(
    #  thecensus,
    # algorithm = 'louvain',
    #d_cutoff = c(10, 20, 30),
    #Lx = Lx,
    #Ly = 500,
    #weighted = TRUE,
    #seed = 0:100
    #) %>%
    #filter(thecensus == 8 | d_cutoff == 20) %>%
    #filter(thecensus == 8 | seed == 0)
    
    cutoff_methods<-'baldeck'
    
    thecensus<-bci%>%pull(census)%>%unique()
    
    parameters=
      expand.grid(
        cutoff=cutoff_methods,
        thecensus=thecensus,
        algorithm='louvain',
        method='ann',
        d_cutoff = c(10,20,30),
        Lx = Lx,
        Ly = 500,
        weighted=TRUE,
        seed=0
      )
    
    
    fdp_analyzed = 
      parameters %>%
      future_pmap_dfr(
        .f = function(
          cutoff,
          thecensus,
          algorithm,
          method,
          d_cutoff,
          Lx,
          Ly,
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
              Ly = Ly,
              method=method
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
    
    dir.create(save_date, showWarnings = FALSE)
    
    
    saveRDS(fdp_analyzed, file = filename)
    
    cluster_data = readRDS(filename)
    filename_consistent = paste0(save_date, 'bci_clustering_analysis_consistent.rds')
    
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
    
    saveRDS(consistent_cluster_data, filename_consistent)
    
  }
}


if(do.kde.analysis){
  #Load bci data
  dir<-getwd()
  setwd(paste0(dir,"/","bci_dryad"))
  
  raw.files<-list.files()
  
  #Load data from DRYAD datasets: Census 7
  mydat<- lapply(raw.files, function(x) {
    load(file = x)
    get(ls()[ls()!= "filename"])
  })
  
  setwd(dir)
  
  all <- do.call("rbind", mydat)%>%tibble()
  bci<-all%>%mutate(census=rep(1:8,sapply(mydat, nrow)))
  bci<-bci%>%select(sp,gx,gy,dbh,census)%>%drop_na()
  
  
  #Baldeck cutoffs
  baldeck_cutoff = 
    bci %>%drop_na()%>%
    group_by(sp) %>%
    summarize(
      baldeck = quantile(dbh, .56), 
      .groups ='drop'
    )
  

    bci %<>%
    inner_join(baldeck_cutoff) %>%
    filter(dbh > baldeck) %>%
    mutate(fdp = 'bci')%>%select(-baldeck)
    

  census_data = bci
  
  consistent_cluster_data = readRDS('20211205bci_clustering_analysis_consistent.rds')
  
  data = 
    census_data %>% 
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
      d_cutoff = unique(data$d_cutoff)
    ) %>%
    future_pmap_dfr(
      .f = KDE,
      .data = data,
      .options = furrr_options(seed = NULL)
    )
  
  kde_full<-kde_full%>%
    inner_join(
      kde_full %>%
        group_by(census, d_cutoff, x, y) %>%
        slice_max(density, n = 1, with_ties = FALSE) %>%
        rename(soiltype = group) %>%
        ungroup() %>%
        select(-density)%>%
        mutate(fdp=fdp)
    )
  
  kde_full$soiltype<-factor(kde_full$soiltype,levels=c(1,2,3,4))
  
}

saveRDS(kde_full,"bci_kde_full.RDS")

## Definition of theta := P(recruit | match) / P(recruit | !match)

if(do.recruitment.analysis){
  
  library(sparr) # for function bivariate.density() in KDE()
  filter = dplyr::filter
  
  ## Determine whether working on SeaWulf (SBU hpc) or personal computer
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 4
  plan(multisession, workers = cores)
  
 
  dir<-getwd()
  setwd(paste0(dir,"/","bci_dryad"))
  
  raw.files<-list.files()
  
  #Load data from DRYAD datasets: Census 7
  mydat<- lapply(raw.files, function(x) {
    load(file = x)
    get(ls()[ls()!= "filename"])
  })
  
  setwd(dir)
  
  all <- do.call("rbind", mydat)%>%tibble()
  bci<-all%>%mutate(census=rep(1:8,sapply(mydat, nrow)))%>%drop_na()
  
  #Baldeck cutoffs
  baldeck_cutoff = 
    bci %>%
    group_by(sp) %>%
    summarize(
      baldeck = quantile(dbh, .56), 
      .groups ='drop'
    )
  
  
  adults= bci %>%
    inner_join(baldeck_cutoff) %>%
    filter(dbh > baldeck) %>%
    mutate(fdp = 'bci')%>%select(-baldeck)
  
  
  
  recruits = 
    adults %>%
    filter(census > 1) %>%
    group_by(treeID) %>%
    slice_min(census) %>%
    ungroup %>%
    mutate(recruit = TRUE)
  
  trees =
    adults %>% 
    left_join(
      recruits, 
      by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
    ) %>%
    replace_na(list(recruit = FALSE))
  
  combined_data =
    trees %>%
    inner_join(
      consistent_cluster_data %>%
        filter(seed == 0), 
      by = c('census', 'sp')
    )
  
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
  for(thecensus in sort(unique(combined_data$census))){
    
    for(thed_cutoff in sort(unique(combined_data$d_cutoff))){
      
      z = NULL
      
      x =
        combined_data %>%
        filter(
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
  
  parms = 
    combined_data_consistent %>%
    select(
      Census = census, 
      Algorithm = algorithm,
      Seed = seed,
      D_cutoff = d_cutoff,
      Group = group
    ) %>%
    unique %>%
    filter(
      Algorithm == 'louvain',
      Seed == 0
    )
  
  
  KernelDensityEstimation = 
    function(gx, gy, Lx = 1000, Ly = 500, quadrat_length = 20, ...){
      
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
  
  
  kde_full = 
    parms %>%
    future_pmap_dfr(
      .f = KDE,
      .data = combined_data_consistent,
      .options = furrr_options(seed = TRUE)
    )
  
  soiltype = 
    kde_full %>%
    group_by(algorithm, seed, d_cutoff, census, x, y) %>%
    slice_max(density, with_ties = FALSE) %>%
    ungroup %>%
    rename(soiltype = group) %>%
    select(-density)
  
  kde_full %<>%
    left_join(
      soiltype,
      by = c('x', 'y', 'census', 'algorithm', 'seed', 'd_cutoff')
    ) %>% 
    mutate(fdp = fdp)
  
  # saveRDS(kde_full, file = kde_filename)
  
  plot_soiltypes = 
    kde_full %>% 
    mutate(soiltype = factor(soiltype)) %>% 
    ggplot(aes(x, y, fill = soiltype)) + 
    geom_tile() + 
    facet_grid(census ~ d_cutoff) + 
    theme(aspect.ratio = ifelse(fdp == 'bci', .5, 1))
  
  plot_soiltypes %>%
    show()
  
  
  rec_df = 
    recruits %>%
    mutate(
      x = seq(10, L - 10, 20)[cut(gx, seq(20, L, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
    ) %>%
    inner_join(
      kde_full %>%
        select(
          labelling_census = census, 
          x, 
          y, 
          algorithm, 
          seed, 
          d_cutoff, 
          soiltype, 
          fdp
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
    mutate(Pm = n / (L * 500 / 20 ^ 2)) %>% 
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
  
  plot_histograms = 
    res %>% 
    ggplot(aes(theta)) + 
    geom_histogram() + 
    facet_wrap(~ d_cutoff)
  
  plot_bars = 
    res_summary %>%
    ggplot(aes(d_cutoff, theta_mean)) +
    geom_hline(yintercept = 1, color = 'grey') +
    geom_errorbar(
      aes(
        x = d_cutoff, 
        ymin = theta_mean, 
        ymax = theta_mean + 2 * theta_se
      )
    ) +
    geom_col(fill = 'plum4') +
    labs(x = 'Distance cutoff', y = 'Theta estimate') +
    theme(aspect.ratio = 1) +
    scale_y_continuous(breaks = 0:6)
  
  plot_histograms %>%
    show()
  
  plot_bars %>%
    show
  
}

if(do.nutrient.analysis){
  
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 4
  plan(multisession, workers = cores)
  
  if(do.data){
    
    library(caret)
    library(C50)
    library(readxl)
    
    if(fdp=='bci'){
      nutrients %<>%
        pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
        group_by(nutrient) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
    }
    
    if(fdp == 'lap') {
      nutrients %<>%
        mutate(gx = gx + 10, gy = gy + 10) %>%  #move from left lower to center
        pivot_longer(-c(gx, gy), names_to = 'nutrient') %>%
        group_by(nutrient) %>%
        rename(x = gx) %>%
        rename(y = gy) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
    }
    
    
    nutrients_wide = 
      nutrients %>%
      select(-value) %>%
      pivot_wider(names_from = nutrient, values_from = standardized)
    
    
    kde = 
      kde_full %>%
      select(-c(group, density)) %>%
      unique()
    
    dtf =
      nutrients_wide %>%
      full_join(kde, by = c('x', 'y'))
    
    parms = 
      kde %>%
      select(census, algorithm, d_cutoff, fdp) %>%
      unique()
    
    indices = 
      if(seawulf){
        seq(nrow(parms)) 
      }else{
        which(parms$census == 7 & parms$d_cutoff == 20) 
      } 
    
    soiltype_v_nutrients = 
      indices %>%
      future_map_dfr(
        .options = furrr_options(seed = TRUE),
        .f = function(index){ 
          data = 
            dtf %>% 
            inner_join(parms[index, ]) %>%
            {
              if (fdp == 'bci') 
                select(.,Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, `N(min)`, pH, soiltype)
              else .
            } %>% 
            {
              if (fdp == 'lap')select(.,Al, Ca, Cu, Fe, K, Mg, Mn, P, Zn, pH, soiltype)
              else .
            } %>%
            mutate(soiltype = factor(soiltype, levels = unique(soiltype)))
          
          C5_model = 
            train(
              soiltype ~ ., 
              data = data, 
              method = 'C5.0',
              trControl = trainControl(method = 'repeatedcv', repeats = 10),
              metric = 'Kappa'
            )
          
          C5_model$results %>%
            as_tibble %>%
            bind_cols(parms[index, ]) %>%
            return()
        }
      )
    
    if(seawulf){
      dir.create(save_directory, showWarnings = FALSE)
      saveRDS(
        soiltype_v_nutrients, 
        file = paste0(save_directory, fdp, '_C5_soiltype_vs_nutrients.rds')
      )
    }
    
    
  }
  
  
  if(do.plots){
    
    soiltype_v_nutrients_summary = 
      soiltype_v_nutrients %>% 
      group_by(d_cutoff) %>% 
      summarize(
        mean = mean(Kappa), 
        se = sd(Kappa) / sqrt(n()), 
        .groups = 'drop'
      )
    
    plot_bars = 
      soiltype_v_nutrients_summary %>% 
      ggplot(aes(d_cutoff, mean)) + 
      geom_errorbar(
        aes(x = d_cutoff, ymin = mean - 2 * se, ymax = mean + 2 * se), 
        width = .5
      ) +
      geom_col(fill = 'plum4') + 
      labs(x = 'Distance cutoff', y = 'Cohen\'s Kappa') +
      coord_cartesian(ylim = c(0, 1)) +
      theme(aspect.ratio = 1)+
      ggtitle('Soil nutrients VS inferred soil types')
    
    
    if (fdp=='bci'){
      nutrients = 
        read_excel(
          nutrients_filename, 
          sheet = 2
        ) %>%
        pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
        group_by(nutrient) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
      
      kde_full = readRDS(kde_filename)
    }
    
    if (fdp == 'lap') {
      nutrients = 
        read_excel(
          nutrients_filename,
          sheet = 4
        ) %>%
        #move from left lower to center
        mutate(gx = gx + 10, gy = gy + 10) %>%
        pivot_longer(-c(gx, gy), names_to = 'nutrient') %>%
        group_by(nutrient) %>%
        rename(x = gx) %>%
        rename(y = gy) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
      kde_full = readRDS(kde_filename)
    }
    
    cor_analysis =
      nutrients %>%
      select(x, y, nutrient, standardized) %>%
      inner_join(
        kde_full %>%
          select(-soiltype) %>%
          rename(soiltype = group), 
        by = c('x', 'y')
      ) %>%
      group_by(
        census, 
        algorithm,
        seed,
        d_cutoff,
        fdp,
        nutrient,
        soiltype
      ) %>%
      summarize(
        cor = cor(density, standardized, use = 'complete.obs'), 
        p.value = cor.test(density, standardized)$p.value, 
        p.adj = p.adjust(p.value, method = 'hochberg'),
        significant = (p.adj < .05),
        cor_sig = ifelse(significant == TRUE, cor, NA),
        .groups = 'drop'
      ) 
    
    ranked_soiltypes = 
      cor_analysis %>%
      group_by(
        algorithm,
        seed,
        fdp,
        d_cutoff,
        census,
        soiltype
      ) %>%
      summarize(
        order = -mean(cor),
        .groups = 'drop'
      ) %>%
      group_by(
        algorithm,
        seed,
        fdp,
        d_cutoff,
        census
      ) %>%
      mutate(ranked_soiltype = factor(rank(order))) %>%
      ungroup()
    
    plot_nutrient_soiltype_correlations =
      cor_analysis %>%
      left_join(
        ranked_soiltypes,
        by = c("census", "algorithm", "seed", "d_cutoff", "fdp", "soiltype")
      ) %>% 
      filter(d_cutoff == 20) %>%
      ggplot(aes(nutrient, cor_sig, group = ranked_soiltype, fill = ranked_soiltype)) + 
      geom_col(position = 'dodge') + 
      facet_grid(census ~ ranked_soiltype, labeller = label_both) +
      labs(fill = 'group', y = 'Pearson correlation coefficient') +
      ggtitle('Distance cutoff = 20 m')
    
    
    plot_kde = 
      kde_full %>% 
      left_join(
        ranked_soiltypes,
        by = c("census", "algorithm", "seed", "d_cutoff", "fdp", "soiltype")
      ) %>%
      filter(d_cutoff %in% c(6, 8, 10, 16, 20, 30)) %>% 
      ggplot(aes(x, y, fill = ranked_soiltype)) + 
      geom_tile() + 
      facet_grid(d_cutoff ~ census, labeller = label_both) + 
      theme(aspect.ratio = .5)
    
    plot_bars %>%
      show
    
    plot_nutrient_soiltype_correlations %>%
      show
    
    plot_kde %>%
      show
    
  }
  
}  


if(do.pca.analysis){
  plants =
    readRDS('bci_kde_nonull.rds'
    ) %>%
    filter(
      census == 7,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(plants) = c('x', 'y', 'g1', 'g2', 'g3', 'g4')
  
  water =
    read.table(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
      ),
      header = TRUE
    ) %>%
    as_tibble() %>%
    mutate(
      x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
    ) %>%
    replace_na(list(x = 10, y = 10)) %>%
    group_by(x, y) %>%
    summarize(water = mean(swp), .groups = 'drop')
  
  nutrients =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bci_nutrient_data.rds?raw=true'
      )
    ) %>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )
  
  data = 
    plants %>%
    inner_join(nutrients) %>%
    inner_join(water)
  
  data_simplified =
    data %>%
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density') %>%
    group_by(x, y) %>%
    slice_max(density) %>%
    ungroup() %>%
    select(-density)
  
  plot_BCI = 
    df_scores %>% 
    ggplot(aes(x, y, fill = group))+
    geom_tile()+
    theme(aspect.ratio = .5) +
    scale_fill_manual(
      values = 
        c(
          colors$red, 
          colors$green, 
          colors$blue, 
          colors$yellow)
    )
  
  set.seed(0)
  pca = 
    nutrients %>% 
    #inner_join(water) %>%
    select(Al:pH) %>% 
    pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)
  
  df_scores = 
    pca@scores %>% 
    as_tibble() %>% 
    bind_cols(data_simplified)  
  
  if(do.kmeans.analysis){
    k = 
      future_pmap_dfr(
        expand_grid(
          groups = 2:10, 
          seed = 0:100
        ),
        function(groups, seed){
          foo = 
            nutrients %>%
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
    
    k4 = 
      kmeans(
        nutrients %>% 
          select(Al:pH),
        centers = 6,
        nstart = 100
      )
    
    df_scores = 
      df_scores %>%
      mutate(kmeans = k4$cluster)
    
    plot_mapping = 
      df_scores %>% 
      count(group, kmeans) %>%
      ggplot(aes(kmeans, n)) +
      geom_col() +
      facet_wrap(~group) +
      labs(
        x = 'soil cluster',
        y = 'quadrat count'
      ) +
      theme(aspect.ratio = 1)
    
    
    
    
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
    
    soiltypes=df_scores%>%
      group_by(kmeans)%>%
      summarize(across(Al:pH,mean))%>%
      ungroup()%>%
      pivot_longer(-kmeans,names_to="nutrients",values_to="mean")
    
    soiltypes$kmeans=as.factor(soiltypes$kmeans)
    
    plot.soiltypes=soiltypes%>%
      ggplot(aes(nutrients,mean,fill=kmeans))+
      geom_bar(stat='identity')+
      facet_grid(~kmeans)
    
    plot.kmeans.geo=df_scores%>%
      mutate(kmeans=as.factor(kmeans))%>%
      ggplot(aes(x,y,fill=kmeans))+
      geom_tile()
    
    plot_grid(plot_gap,plot.soiltypes,labels=c("(A)","(B)"),nrow=2)
    
    library(egg)
    library(gridExtra)
    library(grid)
    
    f1=ggplotGrob(plot_gap)
    f2=ggplotGrob(plot.kmeans.geo)
    f3=ggplotGrob(plot.soiltypes)
    
    fg1 <- gtable_frame(f1)
    fg2 <- gtable_frame(f2)
    
    fg12 <-
      gtable_frame(gtable_cbind(fg1, fg2),
                   width = unit(1, "null"),
                   height = unit(3, "null"))
    
    fg3 <-
      gtable_frame(
        f3,
        width = unit(3, "null"),
        height = unit(1, "null")
      )
    
    combined <- gtable_rbind(fg12, fg3)
    
    grid.newpage()
    grid.draw(combined)
  }
  
  df_loadings = 
    pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(pca@loadings))) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  plot_PCA = 
    ggplot() + 
    geom_point(aes(PC1, PC2, color = group), data = df_scores) + 
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings)), 
        y = rep(0, nrow(df_loadings)), 
        xend = 10*PC1, 
        yend = 10*PC2
      ), 
      data = df_loadings, 
      arrow = arrow(length = unit(.25, "cm"))
    ) +
    geom_text(
      aes(10 * PC1, 10 * PC2, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC1,
      nudge_y = .5 * df_loadings$PC2,
      color = 'darkgreen'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85)) +
    scale_color_manual(
      values = 
        c(
          colors$red, 
          colors$green, 
          colors$blue, 
          colors$yellow)
    )
  
  plot_group = function(gr){
    plot = 
      df_scores %>%
      ggplot() + 
      geom_point(
        aes(PC1, PC2), 
        data = df_scores %>% filter(group != gr),
        color = 'grey',
        alpha = .3
      ) + 
      geom_point(
        aes(PC1, PC2), 
        data = df_scores %>% filter(group == gr),
        color = 'blue'
      ) +
      theme(aspect.ratio = 1) +
      geom_hline(yintercept = 0, color = 'gray') +
      geom_vline(xintercept = 0, color = 'gray') +
      labs(
        x = 'PC1',
        y = 'PC2'
      ) +
      coord_cartesian(
        xlim = range(c(10 * df_loadings$PC1, df_scores$PC1)), 
        ylim = range(c(10 * df_loadings$PC2, df_scores$PC2))
      )
    
    return(plot)
  }
  
  p1 = plot_group('g1')
  p2 = plot_group('g2')
  p3 = plot_group('g3')
  p4 = plot_group('g4')
  
  plot_grid(plot, plot_grid(p1, p2, p3, p4), axis = 'tb')  
  
  s2d = 
    plot(
      df_scores[,1:2], 
      pch = 20, 
      col = (1:4)[match(df_scores$group, c('g1', 'g2', 'g3', 'g4'))]
    )
  s3d = 
    scatterplot3d(
      df_scores[,1:3], 
      pch = 20, 
      color = (1:4)[match(df_scores$group, c('g1', 'g2', 'g3', 'g4'))]
    )
  
  plot_loadings_angle = 
    df_loadings %>%
    mutate(angle = 360 / (2 * pi) * atan(PC2/ PC1)) %>%
    mutate(feature = factor(feature, levels = feature[order(angle)])) %>%
    ggplot(aes(feature, angle, label = feature)) +
    geom_text() +
    theme(legend.position = 'none')
  
  scatterplot3d(
    df_scores$`N(min)`, 
    df_scores$P, 
    df_scores$Fe, 
    pch = 20,
    c(2, 3, 4, 6)[match(df_scores$group, c('g1', 'x', 'x','g4'))]
  )
  
  data %>% 
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density') %>% 
    pivot_longer(Al:water, names_to = 'feature') %>% 
    group_by(feature, group) %>% 
    summarize(
      cor = cor(density, value), 
      .groups = 'drop'
    ) %>% 
    ggplot(aes(feature, cor)) + 
    facet_wrap(~group) + 
    geom_col()
  
  # feature groups
  # group 1: B K Ca Zn Mg N(min) Cu
  # group 2: Mn Fe pH N
  # group 3: Al water P
  
  p124 = 
    df_scores %>% 
    filter(group != 'g3') %>% 
    ggplot(aes(pH, P, color = group)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    scale_color_manual(values = c(colors$red, colors$green, colors$yellow)) +
    theme(aspect.ratio = 1)
  
  pAl = 
    df_scores %>% 
    ggplot(aes(pH, Al, color = group)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    scale_color_manual(values = c(colors$red, colors$green, colors$blue, colors$yellow)) +
    theme(aspect.ratio = 1)
  
  p14 = 
    df_scores %>% 
    filter(group %in% c('g1', 'g4')) %>% 
    ggplot(aes(N, Fe, color = group)) + 
    geom_point() +
    ggtitle('G1 vs G4') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    scale_color_manual(values = c(colors$red, colors$yellow))
  
  p13 = 
    df_scores %>% 
    filter(group %in% c('g1', 'g3')) %>% 
    ggplot(aes(PC1, Al, color = group)) + 
    geom_point() +
    ggtitle('G1 vs G3') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')
  
  p12 = 
    df_scores %>% 
    filter(group %in% c('g1', 'g2')) %>% 
    ggplot(aes(N, P, color = group)) + 
    geom_point() +
    ggtitle('G1 vs G2') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')
  
  p24 = 
    df_scores %>% 
    filter(group %in% c('g2', 'g4')) %>% 
    ggplot(aes(N, P, color = group)) + 
    geom_point() +
    ggtitle('G2 vs G4') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')
  
  plot_grid(p12, p13, p14, p24)
  
}

## test C5.0 on Gaussian random field
if(do.C5.analysis){
  library(RandomFields)
  library(gstat)
  library(caret)
  
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 12
  plan(multisession, workers = cores)
  
  analysis =
    function(
    number_of_nulls,
    kde_full,
    nutrients,
    clusters_run,
    save.result = clusters_run
    ){
      
      dtf_wide =
        kde_full %>%
        group_by(x, y) %>%
        slice_max(density, with_ties = FALSE) %>%
        ungroup() %>%
        inner_join(nutrients%>%select(x,y,Al:water))
      
      if(fdp == 'bci'){
        soil_dtf =
          dtf_wide %>%
          select(Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, `N(min)`, pH, water, soiltype)
      }
      
      if(fdp == 'lap'){
        soil_dtf =
          dtf_wide %>%
          select(Al, Ca, Cu, Fe, K, Mg, Mn, P, Nmin, pH, soiltype)
      }
      
      RandomField =
        function(
    Lx,
    Ly,
    rangepar,
    sillpar,
    nuggetpar,
    seed = seed,
    plot = FALSE
        ) {
          stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
          
          RFoptions(seed = seed)
          
          stress =
            RFsimulate(
              RMgauss(
                scale = rangepar + 1e-16,
                var = 2 * sillpar * (1 - nuggetpar)
              ) + RMtrend(mean = 0) + RMnugget(var = 2 * sillpar * nuggetpar),
              x = 1:Lx,
              y = 1:Ly
            )@data$variable1
          
          if(plot){
            plot(
              raster::raster(matrix(stress, Lx, Ly)),
              las = 1,
              xlab = 'x-coordinate',
              ylab = 'y-coordinate'
            )
          }
          
          return(stress)
          
        }
      
      Gaussian_vgm_optim =
        function(parms, sample.vgm){
          range = parms[1]
          sill = parms[2]
          nugget = parms[3]
          dist = c(0, sample.vgm$dist)
          observed = c(0, sample.vgm$gamma)
          predicted = (sill - nugget) * (1 - exp(-dist ^ 2 / range ^ 2)) + nugget * (dist > 0)
          sum_squared_errors = sum((observed - predicted) ^ 2)
          return(sum_squared_errors)
        }
      
      sample.vgm =
        variogram(
          soiltype ~ 1,
          data = dtf_wide,
          locations = ~ x + y,
          width = 1
        )
      
      fitted.vgm =
        optim(
          par = c(range = 20, sill = .6, nugget = 20),
          fn = Gaussian_vgm_optim,
          sample.vgm = sample.vgm,
          method = "L-BFGS-B",
          lower = c(1, 1, 0)
        )$par
      
      sill = as.numeric(fitted.vgm[2])		## semivariance at the landscape level --> regional heterogeneity
      nugget = as.numeric(fitted.vgm[3])	## semivariance at distance = 1	--> local uniformity
      range = as.numeric(fitted.vgm[1])		## distance at which the semivariance reaches 63% of the sill
      range95 = sqrt(3) * range   ## distance at which the semivariance reaches 95% of the sill
      
      randomfield_test =
        function(
    seed,
    Lx,
    Ly,
    date,
    .data,
    .range95,
    .sill,
    .func
        ){
          set.seed(seed)
          
          test_data =
            expand_grid(
              y = seq(10, Ly - 10, by = 20),
              x = seq(10, Lx - 10, by = 20)
            ) %>%
            mutate(
              field =
                RandomField(
                  Lx = Lx / 20,
                  Ly = Ly / 20,
                  rangepar = range95,
                  sillpar = sill,
                  nuggetpar = .015,
                  seed = seed
                ),
              soiltype =
                cut(
                  field,
                  breaks = quantile(
                    field,
                    p = c(
                      0,
                      dtf_wide %>%
                        count(soiltype) %>%
                        mutate(p = cumsum(n) / sum(n)) %>%
                        pull(p)
                    )
                  ), ## these quantiles match the relative prevalence
                  ## of each soiltype in the real data
                  labels = FALSE,
                  include.lowest = TRUE
                ) %>%
                as.factor
            )
          
          C5_model =
            train(
              soiltype ~ .,
              data = .data %>% mutate(soiltype = test_data$soiltype),
              method = 'C5.0',
              trControl = trainControl(method = 'repeatedcv', repeats = 10),
              metric = 'Kappa'
            )
          
          res =
            bind_cols(
              seed = seed,
              C5_model$results
            ) %>%
            as_tibble
          
          return(res)
          
        }
      
      results =
        1:number_of_nulls %>%
        future_map_dfr(
          .f = randomfield_test,
          Lx = ifelse(fdp == 'bci', 1000, 500),
          Ly = 500,
          date = date,
          .data = soil_dtf,
          .sill = sill,
          .range95 = range95,
          .func = RandomField,
          .options = furrr_options(seed = NULL)
        )
      
      if(save.result){
        
        save_date = gsub('-', '', Sys.Date())
        save_directory = paste0('~/SpatialNiche/Data/', save_date)
        dir.create(save_directory, showWarnings = FALSE)
        saveRDS(
          results,
          file = paste0(save_directory, '/', fdp, '_C5_randomfields.rds')
        )
      }
      
      return(results)
      
    }
  
  results =
    analysis(
      number_of_nulls = 50,
      #number_of_nulls = ifelse(seawulf, 100, 1),
      kde_full = kde_full,
      nutrients = nutrients,
      clusters_run = seawulf
    )
  
}

if(do.st.trait.analysis){
  
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
  
  
  
  bci=readRDS('bci_raw.rds')%>%
    drop_na(dbh)
  
  bci%<>%inner_join(
    bci%>%
      group_by(sp)%>%
      summarize(baldeck=quantile(dbh,0.56))%>%
      ungroup()
  )%>%
    filter(dbh>baldeck)
  
  cluster_data=readRDS("bci_clustering_analysis_consistent_nonulls.rds")
  
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
    ggplot(aes(x,y,fill=density))+
    geom_tile()+
    facet_grid(rows=vars(kmeans))+
    scale_fill_gradientn(colours = terrain.colors(10))
  
  trdat%>%
    select(c(sp,AL:ZN,kmeans))%>%
    pivot_longer(AL:ZN,names_to = 'trait',values_to = 'Concentration')%>%
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
  
  
}

if(do.st.trait.C5.analysis){
  library(caret)
  library(C50)
  library(readxl)
  
  trait_data_st=read.csv('bci_st_tr.csv')%>%as_tibble()
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
  
  
  
  bci=readRDS('bci_raw.rds')%>%
    drop_na(dbh)
  
  bci%<>%inner_join(
    bci%>%
      group_by(sp)%>%
      summarize(baldeck=quantile(dbh,0.56))%>%
      ungroup()
  )%>%
    filter(dbh>baldeck)
  
  cluster_data=readRDS("bci_clustering_analysis_consistent_nonulls.rds")
  
  kde_full=readRDS('kde_full.rds')%>%
    filter(census==7)%>%
    select(x,y,group,density)%>%
    group_by(x,y)%>%
    slice_max(density,with_ties=FALSE)%>%
    ungroup()
  
  
  foo1=bci%>%
    filter(!is.nan(gx))%>%
    mutate(
      x = seq(10, 990, 20)[cut(gx, seq(0, 1000, length.out=51), labels = FALSE,include.lowest=T)],
      y = seq(10, 490, 20)[cut(gy, seq(0,  500, length.out=26), labels = FALSE,include.lowest=T)]
    )%>%
    select(x,y,sp)%>%
    inner_join(foo%>%
                 pivot_longer(
                   AL:ZN,names_to='trait',values_to='conc')
    )%>%
    select(-sp)%>%
    group_by(x,y,trait)%>%
    summarize(conc=mean(conc,na.rm=TRUE))%>%
    ungroup()
  
  c5dat=foo1%>%
    inner_join(kde_full%>%
                 select(x,y,group))%>%
    pivot_wider(names_from = trait,values_from = conc)%>%
    select(-c(x,y))%>%
    mutate(group=as.factor(group))
  
  set.seed(1001) 
  samp=sample(nrow(c5dat),nrow(c5dat)-125)
  
  c5_tr=c5dat[samp,]
  c5_tst=c5dat[-samp,]
  
  cv <- trainControl(method = "repeatedcv", number =10, repeats =10)
  
  
  C5_res=C5.0(c5_tr%>%select(-group),
              c5_tr$group,
              rule=TRUE,
              control=C5.0Control(bands=10,winnow=TRUE),
              trials=1)
  
  
  
  
  C5_model=train(c5_tr%>%select(-soiltype),c5_tr$soiltype,method='C5.0',metric='Kappa',tunelength=10,trControl=cv)
  
  C5.1=train(c5dat%>%select(-group),c5dat$group,method='C5.0',metric='Kappa',tunelength=10,trControl=cv)
}

if(do.trait.analysis){
  
  
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
      normal = shapiro.test(value)$p.value > .05
    )
  
  trait_data =
    foo %>%
    left_join(normality, by = 'trait') %>%
    mutate(standardized = ifelse(normal, value, logged_value)) %>%
    group_by(trait) %>%
    mutate(standardized = scale(standardized)[, 1]) %>%
    ungroup %>%
    select(sp, trait, standardized)
  
  
  pca_data =
    seq(traitlist) %>%
    map_dfr(
      .f = function(i){
        subdat =
          trait_data %>%
          filter(trait %in% traitlist[[i]]) %>%
          select(sp, trait, standardized) %>%
          pivot_wider(names_from = trait, values_from = standardized)
        
        subpca =
          subdat %>%
          pca(method = 'ppca', scale = 'none', center = FALSE)
        
        tibble(
          sp = subdat$sp,
          trait_type = names(traitlist[i]),
          pc1 = subpca@scores[, 1],
          pc2 = subpca@scores[, 2]
        ) %>%
          return
      }
    ) %>%
    left_join(
      trait_data %>%
        select(sp) %>%
        unique,
      by = 'sp'
    )
  
  
  
  trait_tibble =
    tibble(
      trait_type = rep(names(traitlist), lengths(traitlist)),
      trait = unlist(traitlist)
    )
  
  pc1_sign =
    trait_data %>%
    inner_join(trait_tibble) %>%
    inner_join(pca_data) %>%
    group_by(trait, trait_type) %>%
    summarize(
      cor_sign = sign(cor(pc1, standardized)),
      .groups = 'drop'
    ) %>%
    group_by(trait_type) %>%
    summarize(
      cor_sign = round(mean(cor_sign)),
      .groups = 'drop'
    )
  
  pca_data %<>%
    inner_join(pc1_sign) %>%
    mutate(pc1 = pc1 * cor_sign)
  
  nutrients %<>%
    pivot_longer(-c(x, y), names_to = 'nutrient') %>%
    group_by(nutrient) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup
  
  all_data =
    census_data %>%
    filter(dbh >= 100) %>%
    inner_join(
      cluster_data %>%
        filter(seed == 0)
    ) %>%
    mutate(
      x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
    ) %>%
    inner_join(
      kde_full %>%
        select(x, y, census, d_cutoff, soiltype, fdp) %>%
        unique()
    ) %>%
    filter(census == 7, d_cutoff == 20) %>%
    inner_join(nutrients) %>%
    select(census, d_cutoff, sp, gx, gy, x, y, group, soiltype, nutrient, standardized) %>%
    inner_join(pca_data)
  
  correlations =
    all_data %>%
    group_by(trait_type, nutrient) %>%
    summarize(
      cor = cor(standardized, pc1),
      pval = cor.test(standardized, pc1)$p.value,
      .groups = 'drop'
    )
  
  cor_analysis =
    nutrients %>%
    select(x, y, nutrient, standardized) %>%
    inner_join(
      kde_full %>%
        select(-soiltype) %>%
        rename(soiltype = group),
      by = c('x', 'y')
    ) %>%
    group_by(
      census,
      d_cutoff,
      fdp,
      nutrient,
      soiltype
    ) %>%
    summarize(
      cor = cor(density, standardized, use = 'complete.obs'),
      p.value = cor.test(density, standardized)$p.value,
      p.adj = p.adjust(p.value, method = 'hochberg'),
      significant = (p.adj < .05),
      cor_sig = ifelse(significant == TRUE, cor, NA),
      .groups = 'drop'
    )
  
  mean_correlation =
    cor_analysis %>%
    filter(census == 7, d_cutoff == 20, !nutrient %in% c('Al', 'pH')) %>%
    group_by(soiltype) %>%
    summarize(mean_cor = mean(cor), .groups = 'drop') %>%
    mutate(
      group_ID =
        factor(
          soiltype,
          levels = soiltype[order(mean_cor, decreasing = TRUE)]
        )
    )
  
  all_data %<>%
    inner_join(
      mean_correlation %>%
        select(group = soiltype, group_ID)
    )
  
  cor_analysis %<>%
    inner_join(
      mean_correlation %>%
        select(soiltype, group_ID)
    )
  
  species_group_v_trait =
    all_data %>%
    select(sp, trait_type, pc1, group_ID) %>%
    unique
  
  species_group_v_trait_wide =
    species_group_v_trait %>%
    pivot_wider(names_from = trait_type, values_from = pc1) %>%
    arrange(sp)
  
  abundances =
    all_data %>%
    select(sp, gx, gy) %>%
    unique() %>%
    count(sp)
  
  group_levels =
    species_group_v_trait %>%
    pull(group_ID) %>%
    levels()
  
  plot_correlations =
    correlations %>%
    filter(pval <= .05) %>%
    ggplot(aes(trait_type, nutrient, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue')
  
  plot_correlations %>%
    show()
  
  plot_group_densities =
    all_data %>%
    filter(trait_type == 'vital', nutrient == 'Fe') %>%
    ggplot(aes(gx, gy)) +
    geom_density_2d_filled() +
    theme(aspect.ratio = .5) +
    facet_wrap(~ group_ID)
  
  plot_trait_violins_trees =
    all_data %>%
    filter(nutrient == 'Fe') %>%
    ggplot(aes(group_ID, pc1, fill = group_ID)) +
    geom_violin(draw_quantiles = .5) +
    facet_wrap(~ trait_type)
  
  plot_trait_violins_species =
    all_data %>%
    filter(nutrient == 'Fe') %>%
    select(sp, group_ID, pc1, trait_type) %>%
    unique() %>%
    ggplot(aes(group_ID, pc1, fill = group_ID)) +
    geom_violin(draw_quantiles = .5) +
    facet_wrap(~ trait_type)
  
  plot_group_nutrient_correlations =
    cor_analysis %>%
    filter(census == 7, d_cutoff == 20) %>%
    ggplot(aes(nutrient, cor, fill = group_ID)) +
    geom_col() +
    facet_wrap(~ group_ID)
  
  plot_trait_space =
    species_group_v_trait_wide %>%
    ggplot(aes(vital, wood, color = group_ID)) +
    geom_point(size = 4) +
    theme(aspect.ratio = 1)
  
  
  table =
    species_group_v_trait %>%
    filter(trait_type %in% c('vital', 'wood', 'leaf')) %>%
    group_by(trait_type) %>%
    summarize(
      group_ID_x = rep(group_levels[-1], 3),
      group_ID_y = rep(group_levels[-length(group_levels)], each = 3),
      pval = as.numeric(pairwise.wilcox.test(pc1, group_ID)$p.value),
      .groups = 'drop'
    )
  
  if(do.kmeans.analysis){
    
    bci=readRDS('bci_raw.rds')%>%
      drop_na(dbh)
    
    bci%<>%inner_join(
      bci%>%
        group_by(sp)%>%
        summarize(baldeck=quantile(dbh,0.56))%>%
        ungroup()
    )%>%
      filter(dbh>baldeck)
    
    cluster_data=readRDS("bci_clustering_analysis_consistent_nonulls.rds")
    
    splist=unique(bci$sp)
    
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
    
    
    
    kmeans_dat=bci%>%
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
          seed = 0:100
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
      ylab('Trait PC2')
    
    
    
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
    
    plot_kmeans%>%
      
      
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
    
  }
  
}

#Complementary analyses

#Testing the degree distributions within each group
d_cutoff=20
Lx=1000
Ly=500
method='ann'

data<-bci%>%filter(census==7)

result<-data%>%adjacency_matrix(Lx=Lx,Ly=Ly,d_cutoff=d_cutoff,method=method)
adj<-result$adjacency
abd<-result$abundances
cdc<-combined_data_consistent%>%filter(census==7,d_cutoff=20)%>%select(sp,group)%>%
  group_by(sp)%>%summarize(group=max(group))%>%ungroup()
adj<-adj%>%inner_join(abd)%>%inner_join(cdc)



par(mfrow=c(2,2))
for(i in 1:max(adj$group)){
  temp<-adj%>%filter(group==i)%>%group_by(sp)%>%summarize(binary=sum(binary))%>%ungroup()
  numsp<-nrow(temp)
  temp<-temp$binary
  temp=temp[order(temp,decreasing=TRUE)]
  barplot(temp)
  abline(h=numsp)
  
}


for(j in 1:max(adj$group)){
  temp1<-adj%>%filter(group==j)%>%group_by(sp)%>%summarize(n=max(n))%>%ungroup()
  temp1<-temp1$n
  temp1<-temp1[order(temp1,decreasing=TRUE)]
  barplot(temp1)
  
}


#Family-level association with spatial clusters
load("bci.spptable.rdata")
bci.sp<-bci.spptable%>%as_tibble()
bci.sp<-bci.sp%>%select(sp,Genus,Family)

adj<-adj%>%inner_join(bci.sp)


adj$group<-as.factor(adj$group)
adj$Genus<-as.factor(adj$Genus)
adj$Family<-as.factor(adj$Family)




bci_cluster=readRDS('bci_clustering_analysis.rds')
yas_cluster=readRDS('yas_clustering_analysis.rds')

sps=intersect(unique(bci_cluster$sp),unique(yas_cluster$sp))

bci_dat=bci_cluster%>%
  filter(sp%in%sps)%>%
  select(sp,group)



bci_mat=expand.grid(sp1=bci_dat$sp,sp2=bci_dat$sp)%>%
        as_tibble()%>%
        bind_cols(expand_grid(x=as.numeric(bci_dat$group),y=as.numeric(bci_dat$group)))%>%
        mutate(bci=as.numeric(((x-y)==0)))%>%
        select(sp1,sp2,bci)
        
yas_dat=yas_cluster%>%
  filter(sp%in%sps)%>%
  select(sp,group)

bci_mat=bci_mat%>%
        bind_cols(expand_grid(x=as.numeric(yas_dat$group),y=as.numeric(yas_dat$group)))%>%
        mutate(yas=as.numeric(((x-y)==0)))%>%
        mutate(bci=as.factor(bci),yas=as.factor(yas))


mat=confusionMatrix(bci_mat$bci,bci_mat$yas)        


#To see the relative growth of trees within vs. outside of their clusters
dir<-getwd()
setwd(paste0(dir,"/","bci_dryad"))
raw.files<-list.files()
#Load data from DRYAD datasets: Census 7
mydat<- lapply(raw.files, function(x) {
  load(x)
  get(ls()[grep('tree',ls())])
})
setwd(dir)
all <- do.call("rbind", mydat)%>%tibble()

cluster_data=readRDS("bci_clustering_analysis_consistent_nonulls.rds")


#nutrients

water =
  read.table(
    url(
      'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
    ),
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


nutrients=readRDS("bci_nutrient_data.rds")%>%as_tibble()%>%
  mutate(
    across(
      Al:pH, 
      function(x) scale(x)[, 1]
    )
  )


nutrients%<>%inner_join(water)

k4 = 
  kmeans(
    nutrients %>% 
      select(Al:water),
    centers = 4,
    nstart = 100
  )

nutrients%<>%mutate(kmeans=k4$cluster)


bci_raw=readRDS('bci_raw.rds')
bci=bci_raw%>%filter(census%in%c(7,8))%>%select(treeID,sp,gx,gy,dbh,agb,census)%>%drop_na()%>%
        mutate(
          x = seq(10, 990, 20)[cut(gx, breaks = seq(0, 1000, 20), labels = FALSE)],
          y = seq(10, 490, 20)[cut(gy, breaks = seq(0, 500, 20), labels = FALSE)]
          )

baldeck_cutoff =
  bci %>%
  group_by(sp) %>%
  summarize(
    baldeck = quantile(dbh, .56),
    .groups ='drop'
  )

  bci %<>%
  inner_join(baldeck_cutoff) %>%
  filter(dbh > baldeck)%>%
    select(-baldeck)

  
  bci_gr=bci%>%
          select(treeID,census,dbh)%>%
         pivot_wider(names_from=census,values_from=dbh)%>%
         drop_na()
  
  
bci_gr=bci_gr%>%
        mutate(gr=log(bci_gr$`8`/bci_gr$`7`))

dat=bci_gr%>%
    select(treeID,gr)%>%
    left_join(bci%>%
                filter(census==7))

dat%<>%
  group_by(sp)%>%
  mutate(gr1=ifelse(gr<=0,0.5*min(abs(gr)),gr))%>%
  mutate(gr1=gr-min(abs(gr)))%>%ungroup()

dat=dat%>%
  left_join(nutrients)%>%
  left_join(cluster_data%>%
            filter(!duplicated(.))%>%  
            filter(census==7)%>%  
            select(sp,group))%>%drop_na()%>%
  mutate(kmeans=as.factor(kmeans))

correlations=dat%>%
              filter(gr1>-0.5 & gr1<0.5)%>%
              pivot_longer(cols=Al:water,names_to='nutrients',values_to="conc")%>%
              group_by(group,nutrients)%>%
              summarize(
                cor = cor(gr1, conc, use = 'complete.obs'),
                p.value = cor.test(gr1, conc)$p.value,
                p.adj = p.adjust(p.value, method = 'hochberg'),
                significant = (p.adj < .05),
                cor_sig = ifelse(significant == TRUE, cor, NA),
                .groups = 'drop')%>%
                mutate(cor=ifelse(significant,cor,NA))

correlations%>%
  ggplot(aes(nutrients,cor,fill=group))+
  geom_bar(stat='identity')+
  facet_wrap(~group)+
  ylab('Correlation coefficient')


dat%>%
  ggplot(aes(x,y,fill=gr1))+
      geom_tile()+
      facet_wrap(~group)+
      scale_fill_gradient(
        low='red',
        high='blue',
        na.value='grey50', 
        space = "Lab",
        guide = "colourbar",
        aesthetics = "fill"
      )

dat%>%
  filter(gr1>-0.5 & gr1<0.5)%>%
  ggplot(aes(kmeans,gr1,col=kmeans))+
  geom_violin()+
  facet_wrap(~group)
  

dat1=dat%>%
  select(gr1,group,kmeans)%>%
  group_by(group,kmeans)%>%
  summarize(means=mean(gr1),sd=sd(gr1))%>%
  ungroup()

dat1%>%ggplot(aes(kmeans,means))+
      geom_bar(stat='identity')+
      geom_errorbar(aes(ymin=means-sd,ymax=means+sd))+
      facet_wrap(~group)


#Plot soil variable distributions
bci.range=readRDS("bci_nutrient_data.rds")%>%
        as_tibble()%>%
        mutate(fdp='BCI')%>%
        select(-c(x,y))%>%
        pivot_longer(Al:pH,names_to='Variable',values_to='value')

lap.range=readRDS("lap_nutrient_data.rds")%>%
  as_tibble()%>%
  mutate(fdp='La Planada')%>%
  select(-c(gx,gy))%>%
  pivot_longer(Al:pH,names_to='Variable',values_to='value')

yas.range=read.csv(
  url("https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/resoilnutrientdatarequest/yas_20x20_soil.csv?raw=true"
  )
)%>%as_tibble()%>%
  as_tibble()%>%
  mutate(fdp='Yasuni')%>%
  select(-c(x,y))%>%
  pivot_longer(Al:pH,names_to='Variable',values_to='value')

nut.range=rbind(bci.range,yas.range,lap.range)%>%as_tibble()

nut.range%>%
  ggplot(aes(x=value,col=fdp))+
  geom_density()+
  facet_wrap(~Variable,scale='free')


#Plot heatmaps of soil concentrations
water =
  read.table(
    url(
      'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
    ),
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

nut.bci=readRDS("bci_nutrient_data.rds")%>%as_tibble()%>%
  mutate(
    across(
      Al:pH, 
      function(x) scale(x)[, 1]
    )
  )
nut.bci%<>%inner_join(water)

nut.lap=readRDS("lap_nutrient_data.rds")%>%as_tibble()%>%
  mutate(
    across(
      Al:pH, 
      function(x) scale(x)[, 1]
    )
  )

names(nut.lap)[1:2]<-c("x","y")



nut.yas=read.csv(
  url("https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/resoilnutrientdatarequest/yas_20x20_soil.csv?raw=true"
  )
)%>%as_tibble()%>%
  mutate(
    across(
      Al:pH, 
      function(x) scale(x)[, 1]
    )
  )

nut.yas%<>%
  mutate(across(
    c(x,y),
    function(x) 20*x-10
  ))

#BCI
nut.bci%>%pivot_longer(cols=Al:water,names_to='nutrient',values_to='value')%>%
  ggplot(aes(x,y,fill=value))+
  geom_tile()+
  scale_fill_gradient(low="white",high="red")+
  facet_wrap(~nutrient)

#La Planada
nut.lap%>%pivot_longer(cols=Al:pH,names_to='nutrient',values_to='value')%>%
  ggplot(aes(x,y,fill=value))+
  geom_tile()+
  scale_fill_gradient(low="white",high="red")+
  facet_wrap(~nutrient)

#Yasuni
nut.yas%>%pivot_longer(cols=Al:pH,names_to='nutrient',values_to='value')%>%
  ggplot(aes(x,y,fill=value))+
  geom_tile()+
  scale_fill_gradient(low="white",high="red")+
  facet_wrap(~nutrient)
