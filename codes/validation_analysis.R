
library(tidyverse)
library(magrittr)
library(RandomFields)
library(furrr)
library(parallel)


do.data = 1
do.plots = 0

Lx=1000
Ly=500


source('clustering_functions_rann.R')
  
  generate_landscape =
    function(
    Lx,
    Ly,
    quadrat_length,
    rangepar,
    sillpar,
    nuggetpar,
    num_soiltypes,
    seed
    ) {
      stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
      RFoptions(seed = seed)
      
      if(num_soiltypes == 1){
        soiltype = 1
      }else{
        stress = 
          RFsimulate(
            RMgauss(
              scale = rangepar + 1e-16,
              var = 2 * sillpar * (1 - nuggetpar)
            ) + 
              RMtrend(mean = 0) + 
              RMnugget(var = 2 * sillpar * nuggetpar),
            x = seq(Lx / quadrat_length),
            y = seq(Ly / quadrat_length)
          )@data$variable1
        
        soiltype = 
          stress %>%
          cut(
            num_soiltypes, 
            labels = FALSE
          )
      }
      
      landscape = 
        expand_grid(
          y = seq(Ly / quadrat_length),
          x = seq(Lx / quadrat_length)
        ) %>%
        mutate(
          soiltype = as.factor(soiltype)
        )
      
      return(landscape)
    }
  
  
  parameters=expand.grid(seed=1:10,nsoiltypes=c(2,4,5,10,15))
  
  #Create a folder to save the landscape files
  setwd( "C:/Users/Mihir/Documents/landscapes")
  
  parameters%>%
    future_pmap_dfr(
      .f=function(seed,nsoiltypes){
        
        land=generate_landscape(
          Lx,
          Ly,
          quadrat_length,
          rangepar,
          sillpar,
          nuggetpar,
          nsoiltypes,
          seed
        )%>%
          mutate(seed=seed,
                 nsoiltypes=nsoiltypes)
        
        write.csv(land,paste0("landscape_",seed,"_",nsoiltypes,".csv"))
      }
    )
  
  
  
  recruitment = 
    function(
    soiltype, 
    species_group,
    species_abundance,
    theta
    ){
      prob_mat = 
        outer(
          soiltype, 
          species_group, 
          function(s, g) ifelse(s == g, theta, 1)
        )
      
      species = 
        apply(
          prob_mat, 
          1, 
          function(probs){
            sample(
              seq_along(probs), 
              size = 1, 
              prob = probs * species_abundance
            )
          } 
        )
      
      return(species)
    }
  
  build_recruits = 
    function(
    coords, 
    landscape, 
    species_group,
    species_abundance,
    Lx,
    Ly,
    quadrat_length,
    theta
    ){
      df =
        coords %>%
        mutate(
          x = 
            cut(
              gx, 
              breaks = seq(0, Lx, quadrat_length), 
              labels = FALSE
            ),
          y = 
            cut(
              gy, 
              breaks = seq(0, Ly, quadrat_length), 
              labels = FALSE
            )
        ) %>%
        left_join(
          landscape, 
          by = c('x', 'y')
        ) %>%
        mutate(
          species = 
            soiltype %>%
            recruitment(
              species_group = species_group,
              species_abundance = species_abundance,
              theta = theta
            )
        ) %>%
        return()
    }
  
  ## average intercensus mortality and births are each 
  ## typically ~= 10% of the number of living trees
  bci_all = readRDS('bci_raw.rds')
  
  censuses = 
    expand_grid(
      census1 = 1:7,
      census2 = 1:7
    ) %>%
    filter(census2 == census1 + 1)
  
  ndeaths = 
    censuses %>%
    pmap_dbl(
      .f = function(census1, census2){
        a = 
          bci_all %>% 
          filter(
            census == census1, 
            status == 'A',
            dbh >= 100
          )
        
        d = 
          bci_all %>% 
          filter(
            census == census2, 
            status == 'D'
          )
        
        length(intersect(a$treeID, d$treeID)) %>%
          return
      }
    )
  
  nbirths = 
    censuses %>%
    pmap_dbl(
      .f = function(census1, census2){
        a1 = 
          bci_all %>% 
          filter(
            census == census1, 
            status == 'A',
            dbh < 100
          )
        
        a2 = 
          bci_all %>% 
          filter(
            census == census2, 
            status == 'A',
            dbh >= 100
          )
        
        length(intersect(a1$treeID, a2$treeID)) %>%
          return
      }
    )
  
  mean_births = mean(c(nbirths, ndeaths))
  sd_births = sd(nbirths)
  
  mean_deaths = mean(c(nbirths, ndeaths))
  sd_deaths = sd(ndeaths)
  
  
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
  
  
  abuns = 
    bci %>%
    count(census, sp)
  
  number_of_species = 
    abuns %>% 
    filter(n >= 40) %>%
    pull(sp) %>% 
    unique() %>% 
    length()
  
  average_community_size = 
    bci %>%
    count(census) %>%
    pull(n) %>%
    mean()
  
  parameters =
    expand_grid(
      average_community_size = average_community_size,
      nspecies = number_of_species,
      nsoiltypes = c(2, 3, 4, 5, 6, 10, 15),
      ncensuses = 100,
      d_cutoff = c(10, 20), 
      d_step = 1e-5,
      Lx = 1000,
      Ly = 500,
      quadrat_length = 20,
      rangepar = 20, 
      sillpar = 1, 
      nuggetpar = .001, 
      seed = 0, 
      theta = c(1, 2, 5, 10, 1e5),
      clustering_algorithm = c('louvain', 'walktrap'),
      autolinked_species_only = TRUE,
      weighted = c(FALSE, TRUE),
      self_loops = FALSE
    )
  
  if(!seawulf){
    parameters %<>%
      filter(
        nsoiltypes %in% c(1, 10, 15)
      )
  }

    parameters %<>% 
      filter(
        clustering_algorithm == 'louvain', 
        d_cutoff == 20, 
        rangepar == 20, 
        quadrat_length == 20, 
        weighted == TRUE,
        self_loops == FALSE
      )

  
  simulation = 
    function(
    average_community_size,
    nspecies,
    nsoiltypes,
    ncensuses,
    d_cutoff, 
    d_step,
    Lx,
    Ly,
    quadrat_length,
    rangepar, 
    sillpar, 
    nuggetpar, 
    seed, 
    theta,
    clustering_algorithm,
    autolinked_species_only,
    weighted,
    self_loops
    ){
      census = 0
      
      species_list = seq(nspecies)
      
      species_group = 
        rep(seq(nsoiltypes), each = nspecies / nsoiltypes)
      
      delta = nspecies - length(species_group)
      if(delta > 0){
        species_group = c(rep(1, delta), species_group)
      }
      
      landscape = 
        generate_landscape(
          Lx = Lx, 
          Ly = Ly, 
          quadrat_length = quadrat_length, 
          rangepar = rangepar, 
          sillpar = sillpar, 
          nuggetpar = nuggetpar, 
          seed = seed, 
          num_soiltypes = nsoiltypes
        )
      
      community = 
        tibble(
          gx = runif(average_community_size, min = 0, max = Lx),
          gy = runif(average_community_size, min = 0, max = Ly)
        ) %>%
        build_recruits(
          landscape = landscape,
          species_group = species_group,
          species_abundance = 1,
          Lx = Lx,
          Ly = Ly,
          quadrat_length = quadrat_length,
          theta = theta
        )
      
      data = 
        tibble(
          census = census
        ) %>%
        bind_cols(community)
      
      census_deaths = 
        rnorm(
          ncensuses, 
          mean = mean_deaths,
          sd = sd_deaths
        )
      
      census_births = 
        rnorm(
          ncensuses, 
          mean = mean_births,
          sd = sd_births
        )
      
      while(census < ncensuses){
        census = census + 1
        ndeaths = census_deaths[census]
        nbirths = census_births[census]
        
        abundances = 
          community %>%
          count(species) %>%
          right_join(
            tibble(
              species = seq(nspecies)
            ),
            by = 'species'
          ) %>%
          replace_na(list(n = 0))
        
        births = 
          tibble(
            gx = runif(nbirths, min = 0, max = Lx),
            gy = runif(nbirths, min = 0, max = Ly)
          ) %>%
          build_recruits(
            landscape = landscape,
            species_group = species_group,
            species_abundance = abundances$n,
            Lx = Lx,
            Ly = Ly,
            quadrat_length = quadrat_length,
            theta = theta
          )
        
        community %<>%
          sample_n(size = nrow(community) - ndeaths) %>%
          bind_rows(births)
        
        if(census %% 10 == 0){
          data %<>%
            bind_rows(
              tibble(
                census = census
              ) %>%
                bind_cols(community)  
            )   
        }
        
      }
      
      return(data)
    }
  
  analyze = 
    function(
    data,
    average_community_size,
    nspecies,
    nsoiltypes,
    ncensuses,
    d_cutoff, 
    d_step,
    Lx,
    Ly,
    quadrat_length,
    rangepar, 
    sillpar, 
    nuggetpar, 
    seed, 
    theta,
    clustering_algorithm,
    autolinked_species_only,
    weighted,
    self_loops
    ){
      result = 
        data %>%
        pull(census) %>%
        unique() %>%
        future_map_dfr(
          .options = furrr_options(seed = TRUE),
          .f = function(cen){
            dat = 
              data %>%
              filter(census == cen) %>%
              rename(sp = species)
            
            A = 
              adjacency_matrix(
                dat, 
                autolinked_species_only = autolinked_species_only, 
                d_cutoff = d_cutoff, 
                d_step = d_step, 
                Lx = Lx, 
                Ly = Ly
              )
            
            clusters = 
              cluster_analysis(
                A = A, 
                algorithm = clustering_algorithm, 
                weighted = weighted, 
                self_loops = self_loops
              )
            
            clusters$result %>%
              mutate(census = cen) %>%
              return
          }
        ) %>%
        bind_cols(
          tibble(
            average_community_size,
            nspecies,
            nsoiltypes,
            ncensuses,
            d_cutoff, 
            d_step,
            Lx,
            Ly,
            quadrat_length,
            rangepar, 
            sillpar, 
            nuggetpar, 
            seed, 
            theta,
            autolinked_species_only
          )
        )
      
      return(result)
      
    }
  
  wrapper = 
    function(
    average_community_size,
    nspecies,
    nsoiltypes,
    ncensuses,
    d_cutoff, 
    d_step,
    Lx,
    Ly,
    quadrat_length,
    rangepar, 
    sillpar, 
    nuggetpar, 
    seed, 
    theta,
    clustering_algorithm,
    autolinked_species_only,
    weighted,
    self_loops
    ){
      data = 
        simulation(
          average_community_size,
          nspecies,
          nsoiltypes,
          ncensuses,
          d_cutoff, 
          d_step,
          Lx,
          Ly,
          quadrat_length,
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed, 
          theta,
          clustering_algorithm,
          autolinked_species_only,
          weighted,
          self_loops
        )
      
      result = 
        analyze(
          data,
          average_community_size,
          nspecies,
          nsoiltypes,
          ncensuses,
          d_cutoff, 
          d_step,
          Lx,
          Ly,
          quadrat_length,
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed, 
          theta,
          clustering_algorithm,
          autolinked_species_only,
          weighted,
          self_loops
        )
      
      return(result)
    }
  
  wrapper_simulation = 
    function(
    average_community_size,
    nspecies,
    nsoiltypes,
    ncensuses,
    d_cutoff, 
    d_step,
    Lx,
    Ly,
    quadrat_length,
    rangepar, 
    sillpar, 
    nuggetpar, 
    seed, 
    theta,
    clustering_algorithm,
    autolinked_species_only,
    weighted,
    self_loops
    ){
      data = 
        simulation(
          average_community_size,
          nspecies,
          nsoiltypes,
          ncensuses,
          d_cutoff, 
          d_step,
          Lx,
          Ly,
          quadrat_length,
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed, 
          theta,
          clustering_algorithm,
          autolinked_species_only,
          weighted,
          self_loops
        )
      
      return(
        data %>% 
          bind_cols(
            tibble(
              average_community_size,
              nspecies,
              nsoiltypes,
              ncensuses,
              Lx,
              Ly,
              quadrat_length,
              rangepar, 
              sillpar, 
              nuggetpar, 
              seed, 
              theta
            )
          )
      )
    }
  

  
  


#Create plots
  
  theme_set(theme_bw())
  theme_update(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  result = 
    readRDS('~/SpatialNiche/Data/20210720/simulation_analysis.rds')
  
  plot_ngroups = 
    result %>%
    mutate(number_of_groups = factor(number_of_groups)) %>%
    group_by(nsoiltypes, number_of_groups) %>%
    count() %>%
    ungroup() %>%
    group_by(nsoiltypes) %>%
    mutate(p = n / sum(n) * 100) %>%
    ggplot(aes(number_of_groups, p, fill = number_of_groups)) +
    geom_col(color = 'black') +
    theme(legend.position = 'none') +
    facet_wrap(~nsoiltypes, labeller = label_both) +
    labs(x = 'number of groups', y = 'probability (%)')
  
  plot_modularity =
    result %>% 
    mutate(ngroups = factor(number_of_groups)) %>%
    ggplot(aes(ngroups, modularity, fill = ngroups)) + 
    geom_boxplot() +
    facet_wrap(~nsoiltypes, labeller = label_both) +
    theme(legend.position = 'none') +
    labs(x = 'number of groups')
  
  gridExtra::grid.arrange(plot_ngroups, plot_modularity, nrow = 1 )
  
}
