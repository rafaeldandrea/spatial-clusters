library(parallelDist) ## for function parDist()
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
## a graph then find communities in the graph
library(tidyverse)
library(furrr)  ## for parallel computing
library(parallel)  
library(magrittr)
library(RANN)  ## for neighbor-finding function nn2()
library(sparr)



## pdf of the distance d between two random points in a 
## square of side a, from Johan Philip 2007
euclidean2d_distance_probability = function(d, a){
  if(d > a * sqrt(2)) return(0)
  f =
    ifelse(
      d < a,
      -4 * d / a ^ 3 + pi / a ^ 2 + d ^ 2 / a ^ 4,
      -2 / a ^ 2 + 
        4 / a ^ 2 * asin(a / d) + 
        4 / a ^ 3 * sqrt(d ^ 2 - a ^ 2) - 
        pi / a ^ 2 - 
        d ^ 2 / a ^ 4
    )
  return(2 * d * f)
} 

## pdf of the distance d between two random points in a 
## rectangle of sides a and b > a, from Johan Philip 2007
euclidean2d_distance_probability_ab = function(d, a, b){
  if(a == b) return(euclidean2d_distance_probability(d, a))
  s = d ^ 2
  a2 = a ^ 2
  b2 = b ^ 2
  if(d < 0 | s > a2 + b2) return(0)
  if(a > b){
    c = a
    a = b
    b = c
  }
  if(s <= a2){
    f = - 2 * d / a2 / b - 2 * d / a / b2 + pi / a / b + s / a2 / b2
  }
  if(s > a2 & s <= b2){
    f = - 2 * d / a2 / b - 1 / b2 + 2 / a / b * asin(a / d) + 2 / a2 / b * sqrt(s - a2)
  }
  if(s > b2){
    f = 
      - 1 / b2 + 2 / a / b * asin(a / d) + 2 / a2 / b * sqrt(s - a2) - 
      1 / a2 + 2 / a / b * asin(b / d) + 2 / a / b2 * sqrt(s - b2) -
      pi / a / b - s / a2 / b2
  }
  return(2 * d * f)
}

## randomize adjacency matrix while keeping total number of edges and diagonal intact
randomize_adjacency_matrix = function(nvertices, nedges = 1, seed = 0) {
  set.seed(seed)
  
  edges.max = nvertices * (nvertices - 1) / 2
  
  stopifnot(length(nvertices) == 1 && 1 <= nvertices)
  stopifnot(0 <= nedges & nedges <= edges.max)
  
  index.edges = lapply(list(1:(nvertices - 1)), function(k) rep(k * (k + 1) / 2, nvertices - k)) 
  index.edges = index.edges[[1]] + 1:edges.max
  graph.adjacency = matrix(0, ncol = nvertices, nrow = nvertices)
  graph.adjacency[sample(index.edges, nedges)] = 1
  
  return(graph.adjacency + t(graph.adjacency) + diag(nvertices))
}


adjacency_matrix_parallelDist =
  function(
    dat,
    autolinked_species_only,
    d_cutoff,
    d_step,
    Lx,
    Ly
  ){
    
    ## int p(r) * dr from 0 to d*
    cumulative_null_prob_threshold = pi * d_cutoff^2 / (Lx *Ly)
    
    ## abundance cutoff based on E[n] >= 1, where n is the number of pairs
    ## within a circle of radius d*
    abundance_threshold = round(sqrt(1 / cumulative_null_prob_threshold))
    
    abuns =
      dat %>%
      count(sp) %>%
      arrange(desc(n))
    
    dat_filtered =
      dat %>%
      inner_join(abuns, by = 'sp') %>%
      filter(n >= abundance_threshold)
    
    ## selects species that are more correlated with *themselves* than
    ## a hypothetical random pair of species with equal abundance
    if(autolinked_species_only){
      selected_species =
        unique(dat_filtered$sp) %>%
        map_dfr(
          .f = function(char){
            
            df =
              dat_filtered %>%
              filter(sp == char)
            
            d = parallelDist::parDist(cbind(df$gx, df$gy))
            
            nn12_positive = d[d > 0]
            
            connected =
              mean(nn12_positive <= d_cutoff) -
              cumulative_null_prob_threshold
            
            return(tibble(sp = char, keep = connected > 0))
          }
        ) %>%
        filter(keep == TRUE) %>%
        pull(sp)
      
    }else{
      selected_species = unique(dat_filtered$sp)
    }
    
    dat_selected =
      dat_filtered %>%
      filter(sp %in% selected_species)
    
    distances =
      dat_selected %>%
      select(gx, gy) %>%
      as.matrix() %>%
      parallelDist::parDist() %>%
      as.numeric()
    
    n = nrow(dat_selected)
    k = which(distances <= d_cutoff) - 1
    j = n - 2 - floor(sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5)
    i = k + j + 1 - n * (n - 1) / 2 + (n - j) * ((n - j) - 1) / 2
    inds = j * n + i + 1
    
    foo =
      tibble(
        ind = inds,
        distance = distances[k + 1]
      )
    
    bar =
      arrayInd(foo$ind, .dim = rep(nrow(dat_selected), 2)) %>%
      as_tibble %>%
      rename(treeID1 = V1, treeID2 = V2) %>%
      mutate(
        sp1 = dat_selected$sp[treeID1],
        sp2 = dat_selected$sp[treeID2],
        distance = foo$distance
      ) %>%
      count(sp1, sp2)
    
    adjacency =
      bar %>%
      bind_rows(
        bar %>%
          rename(sp1 = sp2, sp2 = sp1)
      ) %>%
      group_by(sp1, sp2) %>%
      summarize(pairs = sum(n), .groups = 'drop') %>%
      left_join(
        abuns %>%
          rename(sp1 = sp, n1 = n),
        by = 'sp1'
      ) %>%
      left_join(
        abuns %>%
          rename(sp2 = sp, n2 = n),
        by = 'sp2'
      ) %>%
      mutate(
        weighted = pairs / n1 / n2,
        binary = 1 * (weighted > cumulative_null_prob_threshold)
      ) %>%
      select(-c(pairs, n1, n2)) %>%
      complete(
        sp1,
        sp2,
        fill = list(weighted = 0, binary = 0)
      )
    
    return(
      list(
        dat = dat,
        parms =
          tibble(
            distance_cutoff = d_cutoff,
            distance_step = d_step,
            Lx = Lx,
            Ly = Ly,
            abundance_cutoff = abundance_threshold
          ),
        abundances = abuns,
        adjacency = adjacency
      )
    )
  }


adjacency_matrix = 
  function(
    dat, 
    d_cutoff, 
    Lx, 
    Ly,
    method = 'ann',
    autolinked_species_only=FALSE,
    num_neighbors = 300
  ){
    ## int p(r) * dr from 0 to d*
    cumulative_null_prob_threshold = pi * d_cutoff^2 / (Lx *Ly)
    
    ## abundance cutoff based on E[n] >= 1, where n is the number of pairs
    ## within a circle of radius d*
    abundance_threshold = round(sqrt(1 / cumulative_null_prob_threshold))
  
    
    abuns = 
      dat %>%
      count(sp) %>%
      arrange(desc(n))
    
    dat_filtered = 
      dat %>%
      inner_join(abuns, by = 'sp') %>%
      filter(n >= abundance_threshold)
    
    ## selects species that are more correlated with *themselves* than
    ## a hypothetical random pair of species with equal abundance
    if(autolinked_species_only){
      selected_species =
        unique(dat_filtered$sp) %>%
        future_pmap_dfr(
          .f = function(char){
            
            df = 
              dat_filtered %>%
              filter(sp == char)
            
            d = parallelDist::parDist(cbind(df$gx, df$gy))
            
            nn12_positive = d[d > 0]
            
            connected =
              mean(nn12_positive <= d_cutoff) -
              cumulative_null_prob_threshold
            
            return(tibble(sp = char, keep = connected > 0))
          },
          .options = furrr_options(seed = TRUE)
        ) %>%
        filter(keep == TRUE) %>%
        pull(sp)
      
    }else{
      selected_species = unique(dat_filtered$sp)
    }
    

    dat_selected = 
      dat_filtered %>%
      filter(sp %in% selected_species)
    
    sp_pairs = 
      expand_grid(
        sp1 = unique(dat_selected$sp),
        sp2 = unique(dat_selected$sp)
      ) %>%
      filter(sp1 < sp2)
    
    adjacency =
      sp_pairs %>%
      future_pmap_dfr(
        .options = furrr_options(seed = NULL),
        .f = function(sp1, sp2){
          df1 = 
            dat_selected %>% 
            filter(sp == sp1)
          
          n1 = nrow(df1)
          
          df2 = 
            dat_selected %>% 
            filter(sp == sp2)
          
          n2 = nrow(df2)
          
          if(method == 'ann'){
            distances =
              nn2(
                with(df1,cbind(gx,gy)),
                with(df2,cbind(gx,gy)),
                k = min(n1, num_neighbors),
                treetype = "kd",
                searchtype = "radius",
                radius = d_cutoff
              )$nn.dists %>%
              as.vector()
          }else if(method == 'pdist'){
            distances = 
              pdist::pdist(
                X = with(df1, cbind(gx, gy)),
                Y = with(df2, cbind(gx, gy))
              )@dist
          } else stop('Need a valid method to calculate neighbors')
          
          pairs = sum(distances <= d_cutoff)
          
          tibble(
            sp1 = sp1, 
            sp2 = sp2, 
            weighted = (pairs / (n1 * n2)) / cumulative_null_prob_threshold,
            binary = 1 * (weighted > 1)
          ) 
        }
      )
    
    adjacency %<>% 
      bind_rows(
        adjacency %>% 
          rename(sp1 = sp2, sp2 = sp1)
      )
    
    return(
      list(
        dat = dat,
        parms = 
          tibble(
            distance_cutoff = d_cutoff,
            Lx = Lx,
            Ly = Ly,
            abundance_cutoff = abundance_threshold
          ),
        abundances = abuns,
        adjacency = adjacency
      )
    )
  }


cluster_analysis = 
  function(
    A, 
    algorithm,
    weighted=weighted
  ){
    
    if(!weighted){
      weight_parameter = NULL
      
      if(is.matrix(A)){
        adjacency_matrix = A
      }
      
      if(is.list(A)){
        adjacency_matrix = 
          A$adjacency %>% 
          complete(
            sp1, 
            sp2, 
            fill = list(weighted = 0, binary = 0)
          ) %>%
          select(-weighted) %>% 
          pivot_wider(
            names_from = sp2, 
            values_from = binary
          ) %>% 
          select(-sp1) %>% 
          as.matrix  
      }
    } 
    
    if(weighted){
      weight_parameter = TRUE
      
      if(is.matrix(A)){
        adjacency_matrix = A
      }
      
      if(is.list(A)){
        adjacency_matrix = 
          A$adjacency %>% 
          complete(
            sp1, 
            sp2, 
            fill = list(weighted = 0, binary = 0)
          ) %>%
          select(-binary) %>% 
          pivot_wider(
            names_from = sp2, 
            values_from = weighted
          ) %>% 
          select(-sp1) %>% 
          as.matrix 
      }
    } 
    
    graph = 
      graph_from_adjacency_matrix(
        adjacency_matrix, 
        weighted = weight_parameter,
        mode = 'undirected',
        diag = FALSE
      )  
    
    fn = get(paste0('cluster_', algorithm))
    communities = try(fn(graph), silent = TRUE)
    
    if(class(communities) == 'try-error'){
      communities = NA
      
      result = 
        tibble(
          name = seq(nrow(adjacency_matrix)),
          algorithm = algorithm, 
          weighted = weighted,
          group = factor(name),
          modularity = 0,
          number_of_groups = nrow(adjacency_matrix)
        ) 
    }else{
      result = 
        membership(communities) %>%
        enframe %>%
        rename(group = value) %>%
        mutate(
          algorithm = algorithm, 
          weighted = weighted,
          group = factor(group),
          modularity = modularity(communities),
          number_of_groups = length(communities)
        )
    }
    
    
    return(
      list(
        data = A,
        adjacency_matrix = adjacency_matrix,
        parms = 
          tibble(
            algorithm, 
            weighted
          ),
        graph = graph,
        communities = communities,
        result = result
      )
    )
    
  }


find.consistency<-function(mat){

  num.gr<-max(as.numeric(mat$group))
  ref.gr<-max(as.numeric(mat$reference_group))
  
  mat1<-as_tibble(mat)
  mat1<-mat1%>%select(group,int,reference_group)

  init<-mat1%>%
    group_by(group)%>%
    slice_max(int,with_ties=FALSE)%>%ungroup()%>%
    arrange(desc(int))%>%as.data.frame
  
  newgr<-NULL
  
  while(nrow(init)>0){
    
    newgr<-rbind(newgr,c(init[1,1],init[1,3]))
    
    mat1<-mat1%>%filter(group!=init[1,1]) %>% filter(reference_group!=init[1,3])
    
    init<-mat1%>%
      group_by(group)%>%
      slice_max(int,with_ties=FALSE)%>%ungroup()%>%
      arrange(desc(int))%>%as.data.frame
    
  }
  
  foo<-setdiff(seq(num.gr),newgr[,1])
  if(length(foo)>0){
  newgr<-rbind(newgr,cbind(foo,(ref.gr+1):num.gr))
  }
  
  newgr<-as.data.frame(newgr)
  colnames(newgr)<-c('group',"newgroup")
  newgr$group<-as.factor(newgr$group)
  newgr$newgroup<-as.factor(newgr$newgroup)
  newgr<-as_tibble(newgr)
  
  mat<-mat%>%full_join(newgr)
  
return(mat)

}



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

kmeans_manual=function(df,K,weight,nstarts){

  ssdat=NULL
  clusdat=NULL
  nc=NULL
  
  #initial centroids by random
  
  centroid = as.matrix(df[sample(1:nrow(df),k),])

  for(ns in 1:nstarts){
    
    #initiate containers
    distance_matrix = matrix( nrow = nrow(df), ncol = K)
    cluster = NULL

    #compute distance between points and initial centroids
    
    for (i in c(1:K)) {
      distance_matrix[,i]=apply(df,1, function(x) sqrt(sum((x-centroid[i,])^2)))
    }
    
    #initial cluster assignment for each city
    cluster=apply(distance_matrix,1,which.min)
    
    
    #iteration baseline
    old_cluster = vector(length = length(cluster))
    new_cluster = cluster
    
    #iterations
    iter=0
    while (!all(old_cluster == new_cluster) && length(unique(new_cluster==K))) {
      if(iter>=10000) break
      
      #update old cluster assignment
      old_cluster = new_cluster
      
      #update the centroids using weighted average
      
      for (j in c(1:K)){
        
        ind=which(cluster==j)
        
        sub_k=df[ind,]
        
        centroid[j,]=apply(sub_k,2,weighted.mean,w=weight[ind])
        
      }
      
      #compute distance between data points and centroids
      
      for (i in c(1:K)) {
        distance_matrix[,i]=apply(df,1, function(x) sqrt(sum((x-centroid[i,])^2)))
      }
      
      #update new_cluster
      new_cluster = apply(distance_matrix,1,which.min)
      
      iter=iter+1
      nc=table(new_cluster)
     
    }
    
  }
  
  return(
    list(
      sse=min(ssdat),
      cluster=clusdat[,which.min(ssdat)]
    )
  )
}