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
library(sparr)
library(dplyr)


filter = dplyr::filter

source('clustering_functions_rann.R')

if(do.clustering.analysis) {

cores = max(4, detectCores())
plan(multisession, workers = cores)

save_date = gsub('-', '', Sys.Date())

filename = paste0('lap_clustering_analysis.rds')

lap=readRDS("lap_census_data.rds")%>%as_tibble()

Lx=500
Ly=500

baldeck_cutoff = 
  lap %>%
  group_by(sp) %>%
  summarize(
    baldeck = quantile(dbh, .56), 
    .groups ='drop'
  )

dat = 
  lap %>%
  inner_join(baldeck_cutoff) %>%
  filter(dbh > baldeck) %>%
  mutate(fdp = 'lap')

thecensus=unique(lap$census)

parameters = 
  expand_grid(
    census = thecensus,
    algorithm = 'louvain',
    d_cutoff = c(10,20,30),
    Lx = Lx,
    Ly = Ly,
    weighted = TRUE,
    seed = 0
  )

fdp_analyzed = 
  parameters %>%
  future_pmap_dfr(
    .f = function(
      census,
      algorithm,
      d_cutoff,
      Lx,
      Ly,
      weighted,
      seed
    )
      {
      
      data=lap%>%filter(census==census)
      
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
          method='ann'
        ) %>%
        cluster_analysis(
          algorithm = algorithm,
          weighted = weighted
        ) %>%
        pluck('result') %>%
        mutate(
          census=census,
          d_cutoff = d_cutoff,
          seed = seed
        )%>%
        return()
    },
    .options = furrr_options(seed = TRUE)
  ) %>%
  rename(sp = name)

dir.create(save_date, showWarnings = FALSE)


saveRDS(fdp_analyzed, file = filename)

}


if(do.kde.analysis){
  cluster_data=readRDS("lap_clustering_analysis.rds")
  
  lap=readRDS("lap_census_data.rds")
  
  #Baldeck cutoffs
  baldeck_cutoff = 
    lap %>%drop_na()%>%
    group_by(sp) %>%
    summarize(
      baldeck = quantile(dbh, .56), 
      .groups ='drop'
    )
  
  
  dat <-lap%>%
    inner_join(baldeck_cutoff) %>%
    filter(dbh > baldeck) %>%
    mutate(fdp = 'lap')%>%select(-baldeck)
  
  
  census_data = dat
  
  
  data = 
    census_data %>% 
    select(sp,census, gx, gy) %>%
    full_join(
      cluster_data %>%
        filter(seed == 0) %>%
        select(sp,census, d_cutoff, group)
    ) %>%drop_na()%>%
    select( d_cutoff,census, gx, gy, group)
  
  KernelDensityEstimation = 
    function(gx, gy, Lx = 500, Ly = 500, quadrat_length = 20, ...){
      
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
    function(census, d_cutoff, .data){
      
      foo =
        .data %>% 
        inner_join(
          tibble(
            census = census,
            d_cutoff = d_cutoff
          )
        ) 
      
      if(nrow(foo) == 0) return()
      
      foo<-foo %>%
        group_by(group) %>% 
        summarize(
          density = 
            KernelDensityEstimation(gx = gx, gy = gy, Lx = 500),
          .groups = 'drop'
        )
      
      bind_cols(
        census = census, 
        d_cutoff = d_cutoff, 
        group = foo$group, 
        foo$density
      ) %>%
        return()
      
    }
  
  
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
        mutate(fdp='lap')
    )
  
  kde_full$soiltype<-factor(kde_full$soiltype,levels=c(1,2,3,4))
  saveRDS(kde_full,"lap_kde_full.RDS")
  
}


if(do.kde.AlP.analysis){
  
  ## The story here is that our plant groups se gregate in the 3-d space 
  ## whose axes are correlated nutrients + Al + P. Group 3 is associated
  ## with high nutrient concentration, while the other 3 groups are associated
  ## with low nutrient concentration. Groups 2 and 3 distinguish themselves
  ## on the P vs Al plane. Unclear how Group 1 differs from 2 and 3.
  
  plants =
    readRDS("lap_kde_full.rds"
    ) %>%
    filter(
      census == 2,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(plants) = c('x', 'y', 'g1', 'g2', 'g3')
  
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
  
  nutrients =readRDS("lap_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="lap")
  
  names(nutrients)[1:2]<-c("x","y")
  
  nutrients%<>%mutate(x=x+10,y=y+10)
  
  data = 
    plants %>%
    inner_join(nutrients) 
  
  data_simplified =
    data %>%
    pivot_longer(g1:g3, names_to = 'group', values_to = 'density') %>%
    group_by(x, y) %>%
    slice_max(density) %>%
    ungroup() %>%
    select(-density)
  
  data_simplified %>%
    ggplot(aes(Al, P)) +
    geom_point() +
    facet_wrap(~group)
  
  BD = 
    function(plantgroup, h0 = .3, data = data_simplified){
      if(plantgroup %in% unique(data_simplified$group)){
        foo = with(
          data %>% 
            filter(group == plantgroup), 
          ppp(Al, P, xrange = range(Al), yrange = range(P))
        )  
      } else{
        foo = with(
          data,
          ppp(Al, P, xrange = range(Al), yrange = range(P))
        )
      }
      
      bar = as.numeric(bivariate.density(foo, h0 = h0)$z$v)
      
      return(bar)
    }
  
  Al_min = min(data$Al)
  Al_max = max(data$Al)
  P_min = min(data$P)
  P_max = max(data$P)
  
  res = 
    expand_grid(
      Al = seq(Al_min, Al_max, length = 128),
      P = seq(P_min, P_max, length = 128)
    ) %>%
    arrange(Al, P) %>%
    bind_cols(
      tibble(
        g0 = BD('all'),
        g1 = BD('g1'),
        g2 = BD('g2'),
        g3 = BD('g3')
    ))
  
  plot_densities = 
    res %>%
    filter(g0 > 0.01) %>%
    select(-g0) %>%
    pivot_longer(g1:g3, names_to = 'group') %>%
    group_by(group) %>%
    mutate(value = scale(value), .groups = 'drop') %>%
    ggplot() +
    geom_raster(aes(Al, P, fill = value)) +
    geom_point(aes(Al, P), data = data_simplified) +
    facet_wrap(~group) +
    scale_fill_gradientn(colors = terrain.colors(2000)) +
    theme(aspect.ratio = 1) +
    labs(fill = 'scaled\ndensity')
  
  set.seed(0)
  pca = 
    nutrients %>%
    filter(x<500)%>%
    select(Al:pH) %>% 
    pcaMethods::pca(nPcs = 2, scale = 'uv', center = TRUE)
  
  
  df_scores = 
    pca@scores %>% 
    as_tibble() %>% 
    bind_cols(data_simplified) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  df_loadings = 
    pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(pca@loadings))) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  plot = 
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
      arrow = arrow(length = unit(.5, "cm"))
    ) +
    geom_text(
      aes(10 * PC1, 10 * PC2, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC1,
      nudge_y = .5 * df_loadings$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85))
  
  BD_PC = 
    function(plantgroup, h0 = .3, data = df_scores){
      if(plantgroup %in% unique(data_simplified$group)){
        foo = with(
          data %>% 
            filter(group == plantgroup), 
          ppp(PC1, PC2, xrange = range(PC1), yrange = range(PC2))
        )  
      } else{
        foo = with(
          data,
          ppp(PC1, PC2, xrange = range(PC1), yrange = range(PC2))
        )
      }
      
      bar = 
        as.numeric(
          bivariate.density(
            foo, 
            h0 = h0, 
            xy = 
              list(
                x = seq(min(df_scores$PC1), max(df_scores$PC1), length = 128),
                y = seq(min(df_scores$PC2), max(df_scores$PC2), length = 128)
              ))$z$v
        )
      
      return(bar)
    }
  
  res = 
    expand_grid(
      PC1 = seq(min(df_scores$PC1), max(df_scores$PC1), length = 128),
      PC2 = seq(min(df_scores$PC2), max(df_scores$PC2), length = 128)
    ) %>%
    arrange(PC1, PC2) %>%
    bind_cols(
      tibble(
        g0 = BD_PC('all'),
        g1 = BD_PC('g1'),
        g2 = BD_PC('g2'),
        g3 = BD_PC('g3'),
        g4 = BD_PC('g4')
      )
    ) %>%
    replace_na(list(g1 = 0, g2 = 0, g3 = 0, g4 = 0)) %>%
    mutate(across(g1:g4, function(x) (x - min(x)) / (max(x) - min(x))))
  
  
  plot_densities = 
    res %>%
    filter(g0 > 1e-6) %>%
    select(-g0) %>%
    pivot_longer(g1:g4, names_to = 'group') %>%
    ggplot() +
    geom_tile(aes(PC1, PC2, fill = value)) +
    geom_point(aes(PC1, PC2), data = df_scores, alpha = .2) +
    facet_wrap(~group) +
    scale_fill_gradientn(colors = terrain.colors(2000)) +
    theme(aspect.ratio = 1) +
    labs(fill = 'scaled\ndensity') +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings)), 
        y = rep(0, nrow(df_loadings)), 
        xend = 10*PC1, 
        yend = 10*PC2
      ), 
      data = df_loadings, 
      arrow = arrow(length = unit(.5, "cm"))
    ) +
    geom_text(
      aes(10 * PC1, 10 * PC2, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC1,
      nudge_y = .5 * df_loadings$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) 
}

if(do.pca.analysis){
  
  
  plants=readRDS(
    url("https://github.com/rafaeldandrea/Spatial-niche/tree/main/Data/YasuniData/20220115yas_kde_full.RDS?raw=true")
  )%>%
    filter(d_cutoff==20)%>%
    select(x,y,group,density)%>%
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
  
  
  nutrients=read.csv(
    url("https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/resoilnutrientdatarequest/yas_20x20_soil.csv?raw=true"
    )
  )%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )
  
  nutrients%<>%
    mutate(across(
      c(x,y),
      function(x) 20*x-10
    ))
  
  
  data = 
    plants %>%
    inner_join(nutrients) 
  
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
    inner_join(water) %>%
    select(Al:water) %>% 
    pcaMethods::pca(nPcs = 4, scale = 'uv', center = TRUE)
  
  df_scores = 
    pca@scores %>% 
    as_tibble() %>% 
    bind_cols(data_simplified) %>%
    mutate(PC1 = -PC1, PC2 = -PC2) 
  
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
    
    plot_mapping = 
      df_scores %>% 
      count(group, kmeans) %>%
      ggplot(aes(kmeans, n)) +
      geom_col() +
      facet_wrap(~group) +
      labs(
        x = 'soil group',
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
    
    grid.arrange(plot_gap, plot_mapping, ncol = 2)
    
  }
  
  df_loadings = 
    pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(pca@loadings))) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  plot_PCA = 
    ggplot() + 
    #geom_point(aes(PC1, PC2), data = df_scores) + 
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings)), 
        y = rep(0, nrow(df_loadings)), 
        xend = 10*PC2, 
        yend = 10*PC3
      ), 
      data = df_loadings, 
      arrow = arrow(length = unit(.25, "cm"))
    ) +
    geom_text(
      aes(10 * PC2, 10 * PC3, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC2,
      nudge_y = .5 * df_loadings$PC3,
      color = 'darkgreen'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC2',
      y = 'PC3'
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




if(do.nutrient.analysis){
  
  link = 'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/'
  suffix = '?raw=true'
  
  nut.bci=readRDS("bci_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="bci")
  
  nut.lap=readRDS("lap_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="lap")
  
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
    ))%>%
    mutate(fdp="yas")
  
  cols<-names(nut.lap)
  
  nut.bci%<>%select(cols)
  
  nut.yas%<>%select(cols)
  
  
  #########################
  #Plot PCAs for all plots
  
  pc_bci=nut.bci%>%select(-c(x,y))%>% 
    pcaMethods::pca(nPcs = 2, scale = 'uv', center = TRUE)
  
  df_scores_bci = 
    pc_bci@scores %>% 
    as_tibble() %>% 
    mutate(PC1 = -PC1, PC2 = -PC2) 
  
  df_loadings_bci=as.data.frame(pc_bci@loadings)
  
  df_loadings_bci%<>%rownames_to_column(var='feature')
  
  PC_bci<-ggplot() + 
    geom_point(aes(PC1, PC2), data = df_loadings_bci) +
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings_bci)), 
        y = rep(0, nrow(df_loadings_bci)), 
        xend = PC1, 
        yend = PC2
      ), 
      data = df_loadings_bci, 
      arrow = arrow(length = unit(.5, "cm"))
    ) +
    geom_text(
      aes(PC1, PC2, label = feature), 
      data = df_loadings_bci, 
      nudge_x = .1 * df_loadings_bci$PC1,
      nudge_y = .1 * df_loadings_bci$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85))
  
  
  ####
  #Yasuni
  
  pc_yas=nut.yas%>%select(-c(x,y))%>% 
    pcaMethods::pca(nPcs = 2, scale = 'uv', center = TRUE)
  
  df_scores_yas = 
    pc_yas@scores %>% 
    as_tibble()%>%
    mutate(PC1 = -PC1, PC2 = -PC2) 
  
  df_loadings_yas=as.data.frame(pc_yas@loadings)
  
  df_loadings_yas%<>%rownames_to_column(var='feature')
  
  PC_yas<-ggplot() + 
    geom_point(aes(PC1, PC2), data = df_loadings_yas) + 
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings_yas)), 
        y = rep(0, nrow(df_loadings_yas)), 
        xend = PC1, 
        yend = PC2
      ), 
      data = df_loadings_yas, 
      arrow = arrow(length = unit(.5, "cm"))
    ) +
    geom_text(
      aes(PC1, PC2, label = feature), 
      data = df_loadings_yas, 
      nudge_x = .1 * df_loadings_yas$PC1,
      nudge_y = .1 * df_loadings_yas$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85))
  
  
  
  ########
  #La Planada
  
  
  pc_lap=nut.lap%>%select(-c(x,y))%>% 
    pcaMethods::pca(nPcs = 2, scale = 'uv', center = TRUE)
  
  df_scores_lap = 
    pc_lap@scores %>% 
    as_tibble()%>%
    mutate(PC1 = -PC1, PC2 = -PC2) 
  
  df_loadings_lap=as.data.frame(pc_lap@loadings)
  
  df_loadings_lap%<>%rownames_to_column(var='feature')
  
  PC_lap<-ggplot() + 
    geom_point(aes(PC1, PC2), data = df_loadings_lap) + 
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings_lap)), 
        y = rep(0, nrow(df_loadings_lap)), 
        xend = PC1, 
        yend = PC2
      ), 
      data = df_loadings_lap, 
      arrow = arrow(length = unit(.2, "cm"))
    ) +
    geom_text(
      aes(PC1, PC2, label = feature), 
      data = df_loadings_lap, 
      nudge_x = .1 * df_loadings_lap$PC1,
      nudge_y = .1 * df_loadings_lap$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85))
  
  
}



if(do.compararison.pca){
  
  #Load kde data for all plots 
  #############################
  
  bci_plants =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/20211205/20211205_kde_full.rds?raw=true'
      )
    ) %>%
    filter(
      census == 7,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(bci_plants) = c('x', 'y', 'g1', 'g2', 'g3', 'g4')

  
  
  lap_plants<-readRDS("lap_kde_full.RDS")%>%
    filter(
      census == 2,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(lap_plants) = c('x', 'y', 'g1', 'g2', 'g3')
  
  
  yas_plants<-readRDS("yas_kde_full.rds")%>%
    filter(
      census == 7,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)

  names(yas_plants) = c('x', 'y', 'g1', 'g2', 'g3','g4')
  

  
  #Load nutrient datasets
  ##############################################################
  
  nut.bci=readRDS("bci_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="bci")
  
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
  
  
  
  nut.lap=readRDS("lap_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="lap")
  
  names(nut.lap)[1:2]<-c("x","y")
  
  nut.lap%<>%mutate(x=x+10,y=y+10)
  
  
  
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
    ))%>%
    filter(x<500)%>%
    mutate(fdp="yas")
  

  
  #BD function
  ################################
  BD = 
    function(plantgroup, h0 = .3, data = data_simplified){
      if(plantgroup %in% unique(data_simplified$group)){
        foo = with(
          data %>% 
            filter(group == plantgroup), 
          ppp(Al, P, xrange = range(Al), yrange = range(P))
        )  
      } else{
        foo = with(
          data,
          ppp(Al, P, xrange = range(Al), yrange = range(P))
        )
      }
      
      bar = as.numeric(bivariate.density(foo, h0 = h0)$z$v)
      
      return(bar)
    }
  
  
  BD_PC = 
    function(plantgroup, h0 = .3, data = df_scores){
      if(plantgroup %in% unique(data_simplified$group)){
        foo = with(
          data %>% 
            filter(group == plantgroup), 
          ppp(PC1, PC2, xrange = range(PC1), yrange = range(PC2))
        )  
      } else{
        foo = with(
          data,
          ppp(PC1, PC2, xrange = range(PC1), yrange = range(PC2))
        )
      }
      
      bar = 
        as.numeric(
          bivariate.density(
            foo, 
            h0 = h0, 
            xy = 
              list(
                x = seq(min(df_scores$PC1), max(df_scores$PC1), length = 128),
                y = seq(min(df_scores$PC2), max(df_scores$PC2), length = 128)
              ))$z$v
        )
      
      return(bar)
    }
  
  #Plot BCI data
  #############################
  set.seed(0)
  
  nutrients=nut.bci%>%
    inner_join(water)
  
  bci_pca = 
    nutrients%>%
    select(Al:water) %>% 
    pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)
  
  df_scores=bci_pca@scores%>%as_tibble()
  
  nutrients=nutrients%>%bind_cols(df_scores)
  
  data = 
    bci_plants %>%
    inner_join(nutrients)%>%
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density')
    

  
  data_simplified =
    data %>%
    group_by(x,y)%>%
    slice_max(density)%>%
    ungroup()%>%
    select(-fdp)
    
  #Correlations
  data1<-data%>%
    select(-fdp)%>%
    pivot_longer(cols=Al:water, names_to = "nutrients", values_to="conc")
  
  cor.dat=data1%>%
    group_by(group,nutrients)%>%
    summarize(cor=cor(density,conc),
              p.value = cor.test(density, conc)$p.value,
              p.adj = p.adjust(p.value, method = 'hochberg'),
              significant = (p.adj < .05),
              cor_sig = ifelse(significant == TRUE, cor, NA),
              .groups = 'drop')
  
 bci_cor=cor.dat%>%
   filter(significant)%>%
   ggplot(aes(nutrients,cor,fill=group))+
   geom_bar(stat='identity')+
   facet_wrap(~group)+ylab("Pearson's correlation coefficient")+ggtitle("BCI")
  
bci_cor_pca=
  data%>%
  select(x,y,PC1,PC2,PC3,group,density)%>%
  pivot_longer(cols=PC1:PC3, names_to = "nutrients", values_to="conc")%>%
  group_by(group,nutrients)%>%
  summarize(cor=cor(density,conc),
            p.value = cor.test(density, conc)$p.value,
            p.adj = p.adjust(p.value, method = 'hochberg'),
            significant = (p.adj < .05),
            cor_sig = ifelse(significant == TRUE, cor, NA),
            .groups = 'drop')%>%
  filter(significant)%>%
  ggplot(aes(nutrients,cor,fill=group))+
  geom_bar(stat='identity')+
  facet_wrap(~group)+ylab("Pearson's correlation coefficient")+ggtitle("BCI")

  

  #Plot Yasuni

 set.seed(0)
 
 yas_pca = 
   nut.yas%>%
   select(Al:pH) %>% 
   pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)
 
 df_scores=yas_pca@scores%>%as_tibble()
 
 nutrients=nut.yas%>%bind_cols(df_scores)
 
 data = 
   yas_plants %>%
   inner_join(nutrients)%>%
   pivot_longer(g1:g4, names_to = 'group', values_to = 'density')
 
 
 #Correlations
 data1<-data%>%
   select(-fdp)%>%
   pivot_longer(cols=Al:pH, names_to = "nutrients", values_to="conc")
 
 cor.dat=data1%>%
   group_by(group,nutrients)%>%
   summarize(cor=cor(density,conc),
             p.value = cor.test(density, conc)$p.value,
             p.adj = p.adjust(p.value, method = 'hochberg'),
             significant = (p.adj < .05),
             cor_sig = ifelse(significant == TRUE, cor, NA),
             .groups = 'drop')
 
 yas_cor=cor.dat%>%
   filter(significant)%>%
   ggplot(aes(nutrients,cor,fill=group))+
   geom_bar(stat='identity')+
   facet_wrap(~group)+ylab("Pearson's correlation coefficient")+ggtitle("Yasuni")
 
 yas_cor_pca=
   data%>%
   select(x,y,PC1,PC2,PC3,group,density)%>%
   pivot_longer(cols=PC1:PC3, names_to = "nutrients", values_to="conc")%>%
   group_by(group,nutrients)%>%
   summarize(cor=cor(density,conc),
             p.value = cor.test(density, conc)$p.value,
             p.adj = p.adjust(p.value, method = 'hochberg'),
             significant = (p.adj < .05),
             cor_sig = ifelse(significant == TRUE, cor, NA),
             .groups = 'drop')%>%
   filter(significant)%>%
   ggplot(aes(nutrients,cor,fill=group))+
   geom_bar(stat='identity')+
   facet_wrap(~group)+ylab("Pearson's correlation coefficient")+ggtitle("Yasuni")
 
  bci_df_loadings=
    bci_pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(bci_pca@loadings)))
  
  #Plot La Planada
  ####################################
 set.seed(0)
 
 lap_pca = 
   nut.lap%>%
   select(Al:pH) %>% 
   pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)
 
 df_scores=lap_pca@scores%>%as_tibble()
 
 nutrients=nut.lap%>%bind_cols(df_scores)
 
 data = 
   lap_plants %>%
   inner_join(nutrients)%>%
   pivot_longer(g1:g3, names_to = 'group', values_to = 'density')
 
 
 #Correlations
 data1<-data%>%
   select(-fdp)%>%
   pivot_longer(cols=Al:pH, names_to = "nutrients", values_to="conc")
 
 cor.dat=data1%>%
   group_by(group,nutrients)%>%
   summarize(cor=cor(density,conc),
             p.value = cor.test(density, conc)$p.value,
             p.adj = p.adjust(p.value, method = 'hochberg'),
             significant = (p.adj < .05),
             cor_sig = ifelse(significant == TRUE, cor, NA),
             .groups = 'drop')
 
 lap_cor=cor.dat%>%
   filter(significant)%>%
   ggplot(aes(nutrients,cor,fill=group))+
   geom_bar(stat='identity')+
   facet_wrap(~group)+ylab("Pearson's correlation coefficient")+ggtitle("La Planada")
 
 lap_cor_pca=
   data%>%
   select(x,y,PC1,PC2,PC3,group,density)%>%
   pivot_longer(cols=PC1:PC3, names_to = "nutrients", values_to="conc")%>%
   group_by(group,nutrients)%>%
   summarize(cor=cor(density,conc),
             p.value = cor.test(density, conc)$p.value,
             p.adj = p.adjust(p.value, method = 'hochberg'),
             significant = (p.adj < .05),
             cor_sig = ifelse(significant == TRUE, cor, NA),
             .groups = 'drop')%>%
   filter(significant)%>%
   ggplot(aes(nutrients,cor,fill=group))+
   geom_bar(stat='identity')+
   facet_wrap(~group)+ylab("Pearson's correlation coefficient")+ggtitle("La Planada")
 
 

 
 
 
 
}
  

if(do.kmeans.analysis){
  
  #Load cluster and nutrient data
  ################################
  bci_plants =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/20211205/20211205_kde_full.rds?raw=true'
      )
    ) %>%
    filter(
      census == 7,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(bci_plants) = c('x', 'y', 'g1', 'g2', 'g3', 'g4')
  
  
  bci_500.plants=readRDS("bci500_kde_full.RDS")%>%
    filter(
      census == 7,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(bci_500.plants) = c('x', 'y', 'g1', 'g2', 'g3')
  
  
  
  lap_plants<-readRDS("lap_kde_full.RDS")%>%
    filter(
      census == 2,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(lap_plants) = c('x', 'y', 'g1', 'g2', 'g3')
  
  
  yas_plants<-readRDS("yas_kde_full.rds")%>%
    filter(
      census == 7,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(yas_plants) = c('x', 'y', 'g1', 'g2', 'g3','g4')
  
  
  
  nut.bci=readRDS("bci_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="bci")
  
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
  
  
  
  
  nut.bci.500<-nut.bci%>%filter(x<500)
  
  nut.lap=readRDS("lap_nutrient_data.rds")%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )%>%
    mutate(fdp="lap")
  
  names(nut.lap)[1:2]<-c("x","y")
  
  nut.lap%<>%mutate(x=x+10,y=y+10)
  
  
  
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
    ))%>%
    filter(x<500)%>%
    mutate(fdp="yas")
  
  cols<-names(nut.lap)
  
  nut.bci%<>%select(cols)
  
  nut.yas%<>%select(cols)
  
  ######################################################
  
  plants=yas_plants
  nutrients=nut.yas
  
  set.seed(0)
  pca = 
    nutrients%>%
    select(Al:pH) %>% 
    pcaMethods::pca(nPcs = 3, scale = 'uv', center = TRUE)
  
  data = 
    plants %>%
    inner_join(nutrients) 
  
  data_simplified =
    data %>%
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density') %>%
    group_by(x, y) %>%
    slice_max(density) %>%
    ungroup() %>%
    select(-density)
  
  df_scores = 
    pca@scores %>% 
    as_tibble() %>% 
    bind_cols(data_simplified) %>%
    mutate(PC1 = -PC1, PC2 = -PC2, PC3=-PC3) 

  df_loadings = 
    pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(pca@loadings))) %>%
    mutate(PC1 = -PC1, PC2 = -PC2, PC3=-PC3)
  

  
  k = 
    future_pmap_dfr(
      expand_grid(
        groups = 2:10, 
        seed = 0:100
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
  
  k4 = 
    kmeans(
      df_scores %>% 
        select(Al:pH),
      centers = 5,
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
      x = 'soil group',
      y = 'quadrat count'
    ) +
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
  
  grid.arrange(plot.soiltypes, plot_mapping, ncol = 2)
  
}

