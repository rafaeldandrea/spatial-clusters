## DO NOT SKIP THIS!

## Instructions for downloading required datasets

###############################################################################################################
# Download census data
###############################################################################################################

# Download the dryad repository https://doi.org/10.15146/5xcp-0d46

# Unzip the folder labelled "bci.tree.zip" and

# 1. move the eight .rdata files to the working directory
# 2. move file called "bci.spptable.rdata" to the working directory. If applicable, change the extension of the file from ".gz" to ".rdata"

# Ref: Condit, Richard et al. (2019), Complete data from the Barro Colorado
# 50-ha plot: 423617 trees, 35 years, Dryad, Dataset, https://doi.org/10.15146/5xcp-0d46


###############################################################################################################
# Download morphological/growth trait data
###############################################################################################################

# Navigate to source for trait data: https://doi.org/10.6084/m9.figshare.c.3303654.v1

# Download a text file called "Supplement_20100505.txt" from the repo. and move it to the working directory.

# Ref: Joseph Wright, S.; Kitajima, Kaoru; J. B. Kraft, Nathan; Reich,
# Peter B.; J. Wright, Ian; Bunker, Daniel E.; et al. (2016).
# Functional traits and the growth–mortality trade-off in tropical trees.
# Wiley. Collection. https://doi.org/10.6084/m9.figshare.c.3303654.v1


###############################################################################################################
# Download foliar elemental trait data
###############################################################################################################

# Navigate to source for leaf stoichiometry data: https://doi.org/10.25573/data.23463254

# Move the tab-separated text file called "BCI LEAF ELEMENTAL COMPOSITION with binomials.txt" to the working directory

# Ref: Wright, Joseph (2023). Foliar elemental composition and carbon and nitrogen isotope
# values for 339 woody species from Barro Colorado Island, Panama.
# Smithsonian Tropical Research Institute. Dataset. https://doi.org/10.25573/data.23463254.v1


###############################################################################################################
# Download soil nutrient data
###############################################################################################################

# Navigate to http://ctfs.si.edu/webatlas/datasets/bci/soilmaps/bci.block20.data.xls

# Download the file to the working directory.

# Ref:Hubbell, S.P., Condit, R., and Foster, R.B. 2005. Barro Colorado Forest Census
# Plot Data. URL http://ctfs.si.edu/webatlas/datasets/bci.


###############################################################################################################
# Download soil water data
###############################################################################################################

# Navigate to https://doi.org/10.6084/m9.figshare.c.4372898

# Download full data in a zipped folder named "Kupers_et_al", unzip it, go to the folder called "Output" and
# move the file called "BCI_SWP_map_mid_dry_season_regular.txt" in an R working directory

# Ref:Kupers, S.J., Wirth, C., Engelbrecht, B.M.J. et al. Dry season soil water potential maps of a 50
# hectare tropical forest plot on Barro Colorado Island, Panama. Sci Data 6, 63 (2019).
# https://doi.org/10.1038/s41597-019-0072-z


###############################################################################################################
# Download custom functions
###############################################################################################################

# Navigate to https://github.com/rafaeldandrea/spatial-clusters/blob/main/codes/clustering_functions_rann.R

# Download the raw file to the working directory.

###############################################################################################################


## RUN TIME: Once all raw and processed data files are present in the working directory, the total run time 
## of this script is circa 8 min on a 36-core Intel(R) Xeon(R) W-2295 CPU @ 3.00GHz, 64 GB RAM.


###############################################################################################################
# ======================================== BEGIN ==============================================================

## Required libraries
library(tidyverse)
library(openxlsx)
library(magrittr)
library(furrr)
library(readxl)
library(parallelDist) ## for function parDist()
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into a graph then find communities in the graph
library(RANN)  ## for neighbor-finding function nn2()
library(FactoClass)
library(C50)
library(caret)
library(sparr) # for function bivariate.density() in KDE()
library(pcaMethods)
library(rcompanion)
library(parallel)
library(ggpubr)

# Source important functions
source('clustering_functions_rann.R')

# Parallelize jobs
cores = detectCores() - 2
plan(multisession, workers = cores)

# Plotting settings
theme_set(theme_bw())
theme_update(panel.grid = element_blank(),
             strip.background = element_rect(fill = 'orange'))

cbpalette <-
  c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#C5E772",
    "#4F2CAA"
  )


## Read and process census data

#First, define the dimensions of the plot,since they will be used frequently. 
Lx=1000
Ly=500
raw.files = list.files()[intersect(grep(".rdata", list.files()), grep("bci.tree", list.files()))]
for (i in raw.files) {
  load(i)
}

tree.files = ls()[grep("bci.tree", ls())]

mydat = lapply(tree.files, function(x) get(x))

all = do.call("rbind", mydat) %>% tibble()
bci_raw = all %>% mutate(census = rep(1:8, sapply(mydat, nrow)))
bci = bci_raw %>% select(sp, gx, gy, dbh, census) %>% drop_na() %>% unique()

load('bci.spptable.rdata')

# Read morphological/growth trait data
trait_data_raw = 
  read.table(
    "Supplement_20100505.txt",
    skip = 25,
    sep = "\t",
    header = TRUE
  )


# Read foliar elemental trait data
trait_data_st = 
  read.table(
    "BCI LEAF ELEMENTAL COMPOSITION with binomials.txt",
    sep = "\t",
    header = TRUE
  )


splist = 
  bci %>%
  inner_join(
    bci %>%
      group_by(sp) %>%
      summarize(
        baldeck = quantile(dbh, .56),
        .groups = 'drop'
      ),
    by = 'sp'
  ) %>%
  filter(dbh > baldeck) %>%
  select(sp) %>%
  unique() %>%
  pull()

bci.sp = 
  bci.spptable %>% 
  as_tibble() %>%
  filter(sp %in% splist) %>%
  select(sp, Genus, Species) %>%
  rename(genus = Genus,
         species = Species) %>%
  mutate(latin = paste0(genus, " ", species))


# Read soil nutrient data
nutrients =
  read_xls("bci.block20.data.xls", sheet = "means (20 X 20 m)") %>%
  as_tibble()


# Read soil water data
water =
  read.table('BCI_SWP_map_mid_dry_season_regular.txt', header = TRUE) %>%
  as_tibble() %>%
  mutate(x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)],
         y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]) %>%
  replace_na(list(x = 10, y = 10)) %>%
  group_by(x, y) %>%
  summarize(
    water = mean(swp), 
    .groups = 'drop'
  ) %>%
  mutate(water = scale(water)[, 1])

#Treat water levels as "nutrient"
nutrients = 
  nutrients %>% 
  inner_join(water, by = c('x', 'y'))



##################################################################################################################
#1. Spatial clustering analysis:

# Step 1:For each census, count the number of heterospecific trees in the neighborhood of each tree,
# where the neighborhood is a circle of radius = 10m, 20m and 30m.
# Step 2: Build an adjacency matrix where each entry (Aij)  in a matrix is the number of neighbors
# of species i which belong to species j (significance is determined by the relative no. of conspecific neighbors)
# Step 3: Calculate maximum modularity of the matrix and determine the strength of clustering
# as well as no. of clusters.
#Caution- lengthy execution time! One iteration of parameter combination takes ~5-6 minutes.

Lx = 1000
Ly = 500



#Cutoffs to determine the adult individuals of each species

#Baldeck cutoffs
baldeck_cutoff =
  bci %>%
  group_by(sp) %>%
  summarize(baldeck = quantile(dbh, .56),
            .groups = 'drop')

cutoff_methods <- 'baldeck'

thecensus <- bci %>% pull(census) %>% unique()

parameters =
  expand.grid(
    cutoff = cutoff_methods,
    thecensus = thecensus,
    algorithm = 'louvain',
    d_cutoff = c(10, 20, 30),
    weighted = TRUE,
    seed = 0:100
  )

## WARNING! Long running time! Each parameter combination requires 2-3 minutes of run time. There
# are 2400 total parameter combinations! Reduce the no. of combinations above (no.of seeds)
#if you just want to test the code.

if (!file.exists("bci_clustering_analysis_consistent.csv")) {
  fdp_analyzed =
    parameters[1:10,] %>%
    future_pmap_dfr(
      .f = 
        function(
          cutoff,
          thecensus,
          algorithm,
          d_cutoff,
          weighted,
          seed
        ) {
            
          if (cutoff == 'raw') {
            dat <-
              bci %>%
              filter(dbh > 100) %>%
              mutate(fdp = 'bci')
          }
          if (cutoff == 'baldeck') {
            dat =
              bci %>%
              inner_join(baldeck_cutoff, by = 'sp') %>%
              filter(dbh > baldeck) %>%
              mutate(fdp = 'bci')
          }
          
          data =
            dat %>%
            filter(census == thecensus)
          
          if (seed > 0) {
            data %<>%
              mutate(sp = sample(sp))
          }
          
          result =
            data %>%
            adjacency_matrix(
              d_cutoff = d_cutoff,
              Lx = Lx,
              Ly = Ly) %>% 
            cluster_analysis(
              algorithm = algorithm,
              weighted = weighted) %>% 
            pluck('result') %>%
            mutate(
              census = thecensus,
              d_cutoff = d_cutoff,
              seed = seed,
              cutoff = cutoff
            ) %>%
            return()
        },
      .options = furrr_options(seed = TRUE)
    ) %>%
    rename(sp = name)
  
  
  #Make sure that the inferred clusters across censuses are consistent in terms of their species identity
  
  ref = 
    cluster_data %>%
    filter(seed == 0) %>%
    select(sp, algorithm, census, d_cutoff, number_of_groups) %>%
    unique() %>%
    slice_min(number_of_groups, n = 1, with_ties = FALSE)
  
  ref = 
    cluster_data %>%
    filter(census == 8, d_cutoff == 20) %>% unique()
  
  x =
    cluster_data %>%
    filter(seed == 0) %>%
    select(algorithm, census, d_cutoff, sp, group, number_of_groups)
  
  x0 =
    x %>%
    filter(census == 8, d_cutoff == 20)
  
  x0 = 
    x0 %>%
    inner_join(
      x0 %>%
        count(group) %>%
        mutate(newgr = (length(n) + 1) - rank(n))) %>% 
    select(-group) %>%
    mutate(group = newgr)
 
  par.cal = x %>% select(census, d_cutoff) %>%
    unique() %>% as.matrix()
  
  
  res = NULL
  
  for (i1 in 1:nrow(par.cal)) {
    x1 = x %>% filter(census == par.cal[i1, 1], d_cutoff == par.cal[i1, 2])
    res1 = NULL
    
    for (i2 in 1:length(unique((x0$group)))) {
      foo = x1 %>%
        group_by(algorithm, census, d_cutoff, group) %>%
        summarize(
          int = length(intersect(
            sp, x0 %>% filter(group == i2) %>% pull(sp)
          )) / (x0 %>% filter(group == i2) %>% pull(sp) %>% length()),
          .groups = 'drop'
        ) %>%
        mutate(reference_group = i2)
      
      res1 = res1 %>% bind_rows(foo)
    }
    res1 = find.consistency(res1)
    
    res = res %>% bind_rows(res1)
  }
  
  consistent_cluster_data = 
    cluster_data %>%
    inner_join(res) %>%
    select(-group) %>%
    rename(group = newgroup) %>% 
    select(-int) %>% 
    select(-reference_group)
  
  saveRDS(consistent_cluster_data, "bci_clustering_analysis_consistent.rds")
  
  write.csv(consistent_cluster_data, "bci_clustering_analysis_consistent.csv")
  
}

consistent_cluster_data =
  read.csv('bci_clustering_analysis_consistent.csv') %>%
  as_tibble() %>%
  select(-X) %>%
  unique()


##################################################################################################################
#2. Kernel density estimation: Kernel smoothing of geographic areas of each of the inferred clusters.

if(!file.exists('bci_kde_full.csv')){
  for (.d_cutoff in unique(consistent_cluster_data$d_cutoff)) {
    data_cutoff =
      bci %>%
      select(sp, census, gx, gy) %>%
      full_join(
        consistent_cluster_data %>%
          filter(seed == 0, d_cutoff == .d_cutoff) %>%
          select(sp, census, d_cutoff, group) %>% 
          unique(),
        by = c('sp', 'census')
      ) %>% 
      drop_na() %>%
      select(census, d_cutoff, gx, gy, group) %>% 
      unique()
    
    assign(paste0('data_', .d_cutoff), data_cutoff)
  }
  
  data  =
    data_10 %>%
    bind_rows(data_20) %>%
    bind_rows(data_30)
  
  
  kde_full =
    expand_grid(
      census = unique(data$census),
      d_cutoff = unique(data$d_cutoff)) %>%
    future_pmap_dfr(
      .f = KDE,
      .data = data,
      .options = furrr_options(seed = NULL)
    )
  
  
  saveRDS(kde_full, "bci_kde_full.rds")
  
  write.csv(kde_full, "bci_kde_full.csv")
}

bci_kde_full = 
  read.csv('bci_kde_full.csv') %>% 
  rename(soiltype = group) %>% 
  as_tibble()

#Plot spatial clusters by census and d_cutoffs

bci_clustmap_all = 
  bci_kde_full %>%
  group_by(census, d_cutoff, x, y) %>%
  slice_max(density, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(soiltype = as.factor(soiltype)) %>%
  ggplot(aes(x, y, fill = soiltype)) +
  geom_tile() +
  theme(
    aspect.ratio = 0.5,
    strip.text = element_text(size = 15, face = 'bold'),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  ) +
  labs(fill = 'Species cluster') +
  scale_fill_manual(values = cbpalette[c(1, 2, 3, 4, 5)]) +
  #theme(legend.position="top")+
  facet_grid(census ~ d_cutoff)

bci_clustmap_7 = 
  bci_kde_full %>%
  group_by(census, d_cutoff, x, y) %>%
  slice_max(density, with_ties = FALSE) %>%
  ungroup() %>%
  filter(d_cutoff == 20, census == 7) %>%
  group_by(census, d_cutoff, x, y) %>%
  slice_max(density, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(soiltype = as.factor(soiltype)) %>%
  ggplot(aes(x, y, fill = soiltype)) +
  geom_tile() +
  theme(aspect.ratio = 0.5) +
  labs(fill = 'Species cluster') +
  scale_fill_manual(values = cbpalette[c(1, 2, 3, 4, 5)]) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )




##################################################################################################################
#3. Perform the recruitment analysis to test the survival of the saplings, juveniles
#and adult in and out of their inferred spatial niches (from the kde analysis).
# ## Definition of theta := P(recruit | match) / P(recruit | !match), where "match" means
# finding an individual tree within its own spatial niche.

#Data filters-
#Species- Use only those species for which kde and clustering analysis was performed.

bci_dat = 
  bci_raw %>%
  select(census, sp, treeID, gx, gy, dbh)%>%
  drop_na()

bci_kde =
  read.csv('bci_kde_full.csv') %>%
  as_tibble() %>%
  group_by(census, d_cutoff, x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  rename(soiltype = group) %>%
  select(census, d_cutoff, x, y, soiltype) %>%
  unique()

cluster_data =
  read.csv('bci_clustering_analysis_consistent.csv') %>%
  select(census, d_cutoff, sp, group) %>%
  unique()

lh_cutoff = 
  bci_dat %>%
  group_by(sp) %>%
  summarize(
    juv = quantile(dbh, 0.56),
    adult = quantile(dbh, 0.88)
  ) %>% 
  ungroup()

recruits = 
  bci_dat %>%
  mutate(
    x = seq(10, 990, 20)[cut(gx,seq(0, 1000, 20),labels = FALSE,include.lowest = TRUE)],
    y = seq(10, 490, 20)[cut(gy,seq(0,  500, 20),labels = FALSE,include.lowest = TRUE)]
  ) %>%
  select(census, treeID, sp, x, y, dbh) %>%
  inner_join(lh_cutoff, by = 'sp') %>%
  filter(census > 1) %>%
  mutate(stage = ifelse(dbh < juv, 'sapling', ifelse(dbh < adult, 'juvenile', 'adult'))) %>%
  select(-c(juv, adult)) %>%
  group_by(treeID) %>%
  slice_min(census) %>%
  ungroup() %>%
  inner_join(
    bci_kde %>% 
      filter(census == 1) %>% 
      select(-census) %>%
      unique(),
    by = c('x', 'y'),
    relationship ='many-to-many'
  ) %>%
  inner_join(
    cluster_data %>%
      unique(),
    by = c('census', 'sp', 'd_cutoff'),
    relationship = 'many-to-many'
  )

niche_area = 
  bci_kde %>%
  filter(census == 1) %>%
  select(x, y, census, d_cutoff, soiltype) %>%
  count(d_cutoff, soiltype, name = 'number_of_quadrats') %>%
  mutate(niche_area = 400 * number_of_quadrats) %>% 
  rename(niche = 'soiltype')

matches = 
  recruits %>%
  group_by(stage, census, d_cutoff, group) %>%
  summarize(
    recruits = n(),
    matches = sum(group == soiltype),
    .groups = 'drop'
  ) %>%
  rename(niche = 'group')

theta = 
  matches %>%
  inner_join(
    niche_area,
    by = c('d_cutoff', 'niche')
  ) %>%
  mutate(
    theta = (matches / (recruits) / (niche_area / (Lx * Ly)))
  )

theta_summary = 
  theta %>%
  group_by(stage, d_cutoff) %>%
  summarize(
    n = n(),
    stheta = sd(theta) / sqrt(n),
    theta = mean(theta),
    .groups = 'drop'
  )

plot_theta = 
  theta %>%
  rename(`distance cutoff` = d_cutoff) %>%
  ggplot() +
  geom_boxplot(aes(stage, theta)) +
  facet_wrap(~`distance cutoff`, labeller = label_both) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
  theme(aspect.ratio = 1)

plot_theta_summary = 
  theta_summary %>%
  rename(`distance cutoff` = d_cutoff) %>%ggplot() +
  geom_point(aes(stage, theta)) +
  geom_errorbar(aes(x = stage, ymin = theta - stheta, ymax = theta + stheta), width = .2) +
  facet_wrap(~`distance cutoff`, labeller = label_both) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
  theme(aspect.ratio = 1)


##################################################################################################################
#4. Perform PCA and kmeans clustering to calculate covariation and spatial distribution
# of different soil nutrients.
# Analyse 7th census

#load kde data if not already loaded
kde_full = 
  read.csv('bci_kde_full.csv') %>% 
  as_tibble()

data = 
  kde_full %>%
  filter(census == 7, d_cutoff == 20) %>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  select(x, y, group) %>%
  inner_join(
    nutrients %>% mutate(across(Al:water, function(x) scale(x)[, 1])),
    by = c('x', 'y')
  )

set.seed(0)

pcadat =
  nutrients %>%
  select(Al:water) %>%
  pcaMethods::pca(nPcs = 4, scale = 'uv', center = TRUE)

#Plot pca loadings to visualize the relationships between nutrients
pcaplot = pcadat@loadings %>% as_tibble()
nutnames = names(data %>% select(Al:water))

ggplot() +
  geom_hline(
    yintercept = 0,
    lty = 2,
    color = "grey",
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2,
    color = "grey",
    alpha = 0.9
  ) +
  geom_segment(
    data = pcaplot,
    aes(
      x = 0,
      xend = PC1,
      y = 0,
      yend = PC2
    ),
    arrow = arrow(length = unit(0.025, "npc"), type = "open")
  ) +
  geom_text(
    data = pcaplot,
    aes(
      x = PC1 * 1.15,
      y =  PC2 * 1.15,
      label = nutnames
    ),
    check_overlap = F,
    size = 3
  ) +
  xlab("PC 1") +
  ylab("PC 2")

ggplot() +
  geom_hline(
    yintercept = 0,
    lty = 2,
    color = "grey",
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2,
    color = "grey",
    alpha = 0.9
  ) +
  geom_segment(
    data = pcaplot,
    aes(
      x = 0,
      xend = PC2,
      y = 0,
      yend = PC3
    ),
    arrow = arrow(length = unit(0.025, "npc"), type = "open")
  ) +
  geom_text(
    data = pcaplot,
    aes(
      x = PC2 * 1.15,
      y =  PC3 * 1.15,
      label = nutnames
    ),
    check_overlap = F,
    size = 3
  ) +
  xlab("PC 2") +
  ylab("PC 3")

ggplot() +
  geom_hline(
    yintercept = 0,
    lty = 2,
    color = "grey",
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2,
    color = "grey",
    alpha = 0.9
  ) +
  geom_segment(
    data = pcaplot,
    aes(
      x = 0,
      xend = PC1,
      y = 0,
      yend = PC3
    ),
    arrow = arrow(length = unit(0.025, "npc"), type = "open")
  ) +
  geom_text(
    data = pcaplot,
    aes(
      x = PC1 * 1.15,
      y =  PC3 * 1.15,
      label = nutnames
    ),
    check_overlap = F,
    size = 3
  ) +
  xlab("PC 1") +
  ylab("PC 3")

df_scores =
  pcadat@scores %>%
  as_tibble() %>%
  bind_cols(data)


# WARNING: Long running time! Reduce number of seeds to check code execution
k =
  expand_grid(
    groups = 2:10,
    seed = 0:100
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
    .data = df_scores %>% select(Al:water),
    .options = furrr_options(seed = NULL)
  )

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
  
#4 clusters

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

df_scores =
  df_scores%>%
  inner_join(
    df_scores%>%
      group_by(kmeans)%>%
      summarize(xm=mean(x),ym=mean(y))%>%
      ungroup()%>%
      arrange(xm)%>%
      mutate(newkmeans=1:length(kmeans)),
    by = 'kmeans'
  )
  

bci_soilclus = df_scores %>%
  select(-kmeans)%>%
  rename(kmeans=newkmeans)%>%
  mutate(kmeans = as.factor(kmeans)) %>%
  ggplot(aes(x, y, fill = kmeans)) +
  geom_tile() +
  labs(fill = gsub('\\s', '\n', "Soil-Nutrient Cluster")) +
  scale_fill_manual(values = cbpalette[1:4]) +
  theme(
    aspect.ratio = 0.5,
    legend.key.size = unit(0.4, 'cm'),
    legend.title = element_text(size = 8)
  )

df_scores%>%
  select(-kmeans)%>%
  rename(kmeans=newkmeans)%>%
  mutate(kmeans = as.factor(kmeans)) %>%
  ggplot(aes(PC1,PC2,col=kmeans))+
  geom_point()+
  scale_color_manual(values = cbpalette[c(1, 2, 3, 4)])


bci_cor.kmeans = 
  df_scores %>%
  select(-kmeans)%>%
  rename(kmeans=newkmeans)%>%
  mutate(kmeans = as.factor(kmeans)) %>%
  pivot_longer(Al:water, names_to = 'nutrients', values_to = 'conc') %>%
  group_by(kmeans, nutrients) %>%
  summarize(mean = mean(conc), .groups = 'drop') %>%
  mutate(kmeans = as.factor(kmeans)) %>%
  ggplot(aes(nutrients, mean, fill = kmeans)) +
  geom_col(position = 'dodge') +
  facet_grid(cols = vars(kmeans)) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    aspect.ratio = 0.75
  ) +
  scale_fill_manual(values = cbpalette[c(1, 2, 3, 4)])

#Plot the correlations of spatial clusters and soil nutrients
cordat = 
  kde_full %>%
  filter(census == 7, d_cutoff == 20) %>%
  inner_join(
    nutrients %>%
      pivot_longer(
        cols = Al:water,
        names_to = 'nutrient',
        values_to = 'value'
      ) %>%
      group_by(nutrient) %>%
      mutate(scaled_value = (value - min(value)) / (max(value) - min(value))) %>%
      ungroup(),
    by = c('x', 'y'),
    relationship = 'many-to-many'
  ) %>%
  group_by(group, nutrient) %>%
  summarize(
    cor = cor(density, scaled_value, use = 'complete.obs'),
    p.value = cor.test(density, scaled_value)$p.value,
    p.adj = p.adjust(p.value, method = 'hochberg'),
    significant = (p.adj < .05),
    cor_sig = ifelse(significant == TRUE, cor, NA),
    .groups = 'drop'
  )

correlation_data = 
  kde_full %>%
  filter(census == 7, d_cutoff == 20) %>%
  inner_join(nutrients, by = c('x', 'y')) %>%
  pivot_longer(
    cols = Al:water,
    names_to = 'nutrient',
    values_to = 'value'
  ) %>%
  group_by(nutrient) %>%
  mutate(scaled_value = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup() %>%
  group_by(group, nutrient) %>%
  summarize(
    cor = cor(density, scaled_value, use = 'complete.obs'),
    p.value = cor.test(density, scaled_value)$p.value,
    p.adj = p.adjust(p.value, method = 'hochberg'),
    significant = (p.adj < .05),
    cor_sig = ifelse(significant == TRUE, cor, NA),
    .groups = 'drop'
  )

#Ignore the warning message Removed 5 rows containing missing values (`geom_bar()`)
#as some of the species had missing trait values.
cor_plot=correlation_data %>%
  mutate(group = as.factor(group)) %>%
  ggplot(aes(nutrient, cor_sig, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(vars(group)) +
  theme(
    axis.text = element_text(size = 10, angle = 90),
    axis.title = element_text(size = 12, face = 'bold'),
    legend.title = element_text(size = 12, face = 'bold')
  ) +
  xlab('Nutrients') +
  ylab("Pearson's correlation coefficient") +
  labs(fill = 'Spatial cluster') +
  scale_fill_manual(values = cbpalette[c(1, 2, 3, 6)])



#############################################################################################
#5. Perform C5.0 analysis to find associations between spatial clusters and soil nutrient variables

cluster_data = read.csv("bci_clustering_analysis_consistent.csv") %>% as_tibble()

kde_full = read.csv('bci_kde_full.csv') %>% as_tibble()

kde_full = 
  kde_full %>%
  filter(census == 7, d_cutoff == 20) %>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  unique()

nutrient = nutrients %>%
  pivot_longer(-c(x, y), names_to = 'nutrient') %>%
  group_by(nutrient) %>%
  mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

dtf =
  nutrient %>%
  select(-value) %>%
  pivot_wider(names_from = nutrient, values_from = standardized) %>%
  full_join(
    kde_full %>%
    select(x, y, group),
    by = c('x', 'y')
  ) %>%
  mutate(group = as.factor(group))


#Start the analysis
#Vary the parameter 'mincases' to determine the minimum no. of samples per split
#of the dataset.
#Record the results in terms of kappa
C5_res = C5.0(
  dtf %>% select(Al:water),
  dtf$group,
  rule = TRUE,
  control = C5.0Control(
    bands = 10,
    winnow = TRUE,
    minCases = 10
  ),
  trials = 1
)

cv = 
  trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 10
  )


# WARNING: long running time; will issue warnings that setting row names  
#          on a tibble is deprecated
C5_model = 
  train(
    dtf %>% select(Al:water),
    dtf$group,
    method = 'C5.0',
    metric = 'Kappa',
    tunelength = 10,
    trControl = cv
  )

#Record the kappa results using summary(c5_model)


###########################################################################################
#6. Morphological/growth trait analysis
# Perform PCAs and k-means analysis

trait_data =
  trait_data_raw %>%
  rename(
    genus = GENUS.,
    species = SPECIES.
  ) %>%
  as_tibble() %>%
  inner_join(
    bci.sp %>% select(species, genus, sp),
    by = c('genus', 'species')
  ) %>%
  select(-c('genus', 'FAMILY.', 'species'))

bci_dat = 
  bci %>%
  inner_join(
    bci %>%
      group_by(sp) %>%
      summarize(
        baldeck = quantile(dbh, 0.56),
        .groups = 'drop'
      ),
    by = 'sp'
  ) %>%
  filter(dbh > baldeck)

cluster_data = read.csv("bci_clustering_analysis_consistent.csv")

kde_full=read.csv("bci_kde_full.csv")%>%as_tibble

splist = unique(bci$sp)

size_traits = c('HEIGHT')

leaf_traits = c('LMA')

seed_traits = c('SEEDMASS')

wood_traits = c('WSG')

vital_traits = c("RGR95SAP",
                 "MRT25SAP")

traitlist =
  list(
    vital = vital_traits,
    leaf = leaf_traits,
    seed = seed_traits,
    wood = wood_traits,
    size = size_traits
  )


foo =
  trait_data %>%
  pivot_longer(-sp, names_to = 'trait') %>%
  filter(!is.na(value)) %>%
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


pca_data =
  seq(traitlist) %>%
  map_dfr(
    .f = function(i) {
      subdat =
        trait_data %>%
        filter(trait %in% traitlist[[i]]) %>%
        select(sp, trait, standardized) %>%
        pivot_wider(names_from = trait, values_from = standardized)
      
      if (ncol(subdat) > 2) {
        subpca =
          subdat %>%
          pca(method = 'ppca',
              scale = 'none',
              center = FALSE)
        
        res = tibble(
          sp = subdat$sp,
          trait_type = names(traitlist[i]),
          pc1 = subpca@scores[, 1],
          pc2 = subpca@scores[, 2]
        )
      } else{
        res = tibble(
          sp = subdat$sp,
          trait_type = names(traitlist[i]),
          pc1 = as.vector(unlist(subdat[, 2])),
          pc2 = NA
        )
      }
      
      return(res)
    }
  )

trait_data_w = 
  trait_data %>%
  filter(sp %in% splist) %>%
  pivot_wider(names_from = trait, values_from = standardized)


pca_data2 = 
  trait_data_w %>%
  select(-sp) %>%
  pca(
    nPCs = 4,
    method = 'ppca',
    scale = 'none',
    center = FALSE
  )

trait_data_pca = 
  trait_data_w %>%
  mutate(
    PC1 = pca_data2@scores[, 1],
    PC2 = pca_data2@scores[, 2]
  ) %>%
  inner_join(
    pca_data %>%
      select(sp, trait_type, pc1) %>%
      pivot_wider(names_from = trait_type, values_from = pc1),
    by = 'sp'
  ) %>%
  select(sp, vital, leaf, seed, wood, size, PC1, PC2)

#kmeans analysis to see if the clustering in species traits is associated with
#spatial clustering in species
kmeans_dat = 
  bci_dat %>%
  filter(census == 7) %>%
  inner_join(
    trait_data_pca, 
    relationship = "many-to-many",
    by = 'sp'
  ) %>%
  inner_join(
    cluster_data %>%
      filter(census == 7, d_cutoff == 20) %>%
      select(sp, group),
    relationship = "many-to-many",
    by = 'sp'
  )


abuns = 
  kmeans_dat %>%
  group_by(sp) %>%
  summarize(abun = n()) %>%
  mutate(abun = log(abun))

trdat = 
  trait_data_pca %>%
  inner_join(abuns, by = 'sp')

k_traits =
  expand_grid(
    groups = 2:10,
    seed = 0:100
  ) %>%
  future_pmap_dfr(
    .f = function(groups, seed) {
      foo =
        trdat %>%
        select(PC1:PC2)
      
      wgt = trdat$abun
      
      if (seed > 0) {
        set.seed(seed)
        foo = apply(foo, 2, sample)
      }
      wss =
        sum(
          kmeansW(
            foo,
            centers = groups,
            nstart = 100,
            weight = wgt
          )$withinss
        )
      
      return(
        tibble(
          groups = groups,
          seed = seed,
          wss = wss,
          lwss = log(wss)
        )
      )
    },
    .options = furrr_options(seed = NULL)
  ) 

k_traits %>%
  filter(seed == 0) %>%
  inner_join(
    k_traits %>%
      filter(seed > 0) %>%
      group_by(groups) %>%
      summarize(
        null = mean(lwss),
        sd = sd(lwss) * sqrt(1 + (1 / 100)),
        .groups = 'drop'
      ),
    by = 'groups'
  ) %>%
  mutate(gap = null - lwss) %>%
  ggplot(aes(groups, gap)) +
  geom_line() +
  geom_point() +
  xlab('number of groups') +
  theme(aspect.ratio = 1)

k_traits %>%
  filter(seed == 0) %>%
  inner_join(
    k_traits %>%
      filter(seed > 0) %>%
      group_by(groups) %>%
      summarize(
        null = mean(wss),
        lnull = mean(lwss),
        sdlnull = sd(lwss) * (sqrt(1 + (1 / 100))),
        .groups = 'drop'
      ),
    by = 'groups'
  ) %>%
  mutate(gap = lnull - lwss) %>%
  mutate(
    diffgap = c(NA, diff(gap)),
    res = diffgap - sdlnull
  ) %>%
  ggplot(aes(groups, res)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab('number of groups') +
  theme(aspect.ratio = 1)


#Optimal clusters=4

k4 = 
  kmeansW(
    trdat %>% select(PC1, PC2),
    centers = 4,
    weight = trdat$abun,
    nstart = 100
  )

trdat_clustered =
  trdat %>%
  mutate(kmeans = as.factor(k4$cluster))

trdat_clustered_summary = 
  trdat_clustered %>% 
  group_by(kmeans) %>% 
  summarize(mean_vital = mean(vital, na.rm= TRUE), .groups = 'drop') %>%
  mutate(group = as.character(rank(mean_vital)))

trdat_clustered =
  trdat_clustered %>%
  inner_join(trdat_clustered_summary, by = 'kmeans')

kmeans_dat_clustered =
  kmeans_dat %>%
  inner_join(
    trdat_clustered %>%
    select(sp, kmeans, ordered_group = group),
    by = 'sp'
  )

PC_cluster = 
  trdat_clustered %>%
  ggplot(aes(PC1, PC2, col = group)) + 
  geom_point() +
  xlab('Trait PC1') +
  ylab('Trait PC2') +
  scale_color_manual(values = cbpalette[c(1, 2, 3, 4)]) +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 12, face = 'bold'),
    legend.title = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 12)
  ) +
  ggtitle('A')

#Ignore the warning message about 'removed rows containing non-finite values'. They are due to 
#missing data on some of the species.
boxplot_clusters = 
  trdat_clustered %>%
  pivot_longer(
    cols = vital:size,
    names_to = 'Trait',
    values_to = 'value'
  ) %>%
  group_by(sp, Trait) %>%
  ggplot(aes(Trait, value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = cbpalette) +
  theme(
    axis.title = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = 'bold'),
    aspect.ratio = .5
  ) +
  ylab('Values') +
  ggtitle('B')

#Ignore the warning message about 'removed rows containing non-finite values'. They are due to 
#missing data on some of the species.
trdat_clustered %>%
  pivot_longer(cols = vital:size,
               names_to = 'Trait',
               values_to = 'value') %>%
  ggplot(aes(kmeans, value, fill = Trait)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = cbpalette) +
  theme(
    axis.title = element_text(size = 12, face = 'bold'),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = 'bold')
  ) +
  ylab('Values')


#Plot the PC loadings for the vital rate traits
pca_vital =
  trait_data %>%
  filter(trait %in% traitlist[[1]]) %>%
  select(sp, trait, standardized) %>%
  pivot_wider(names_from = trait, values_from = standardized) %>%
  pca(
    nPCs = 4,
    method = 'ppca',
    scale = 'none',
    center = FALSE
  )

dat1 = 
  trait_data %>%
  filter(trait %in% traitlist[[1]]) %>%
  select(sp, trait, standardized) %>%
  pivot_wider(names_from = trait, values_from = standardized)


trdat_group = 
  trdat_clustered %>%
  inner_join(
    cluster_data %>%
    filter(census == 7, d_cutoff == 20) %>%
    select(sp, group) %>% unique(),
    by = 'sp'
  ) %>%
  select(group, vital, leaf, seed, wood, size) %>%
  pivot_longer(
    cols = vital:size,
    names_to = "trait",
    values_to = "value"
  )


cmpr = list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))

## WARNING: will return warnings about computing exact p-values with ties.
ggboxplot(
  trdat_group,
  x = "group",
  y = "value",
  color = "group") +
  facet_wrap( ~ trait, scales = 'free') +
  stat_compare_means(
    comparisons = cmpr,
    tip.length = 0,
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    ),
    vjusts = -0.1
  ) +
  # Add a border around each facet
  theme(panel.border = element_rect(fill = NA, colour = "black")) +
  xlab("Spatial cluster") +
  # Don't include a legend
  guides(colour = FALSE)

## WARNING: will return warnings about missing values
trdat_clustered %>%
  select(c(sp, vital:size, kmeans)) %>%
  pivot_longer(vital:size, names_to = 'trait', values_to = 'Concentration') %>%
  ggplot(aes(kmeans, Concentration, fill = trait)) +
  geom_col(position = 'dodge')

## WARNING: will return warnings about missing values
trdat_clustered %>%
  select(c(sp, vital:size, kmeans)) %>%
  pivot_longer(vital:size, names_to = 'trait', values_to = 'Concentration') %>%
  ggplot(aes(kmeans, Concentration)) + geom_boxplot() + facet_wrap( ~ trait, scales =
                                                                      'free')

#Plot KDE of k-means clusters

#Ignore the warning messages 'data contain duplicated points'.
kde_kmeans = 
  kmeans_dat_clustered %>%
  select(gx, gy, kmeans) %>%
  group_by(kmeans) %>%
  summarize(
    density =
      KernelDensityEstimation(gx = gx, gy = gy, Lx = 1000),
    .groups = 'drop'
  )


kde_kmeans = 
  bind_cols(
    kde_kmeans$kmeans,
    kde_kmeans$density
  )

names(kde_kmeans) = c('kmeans', 'x', 'y', 'density')

kde_kmeans %>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  ggplot(aes(x, y, fill = kmeans)) +
  geom_tile() +
  theme(aspect.ratio = 0.5)

kde_kmeans %>%
  ggplot(aes(x, y, fill = density)) +
  geom_tile() +
  facet_grid(rows = vars(kmeans)) +
  scale_fill_gradientn(colours = terrain.colors(10))

#Density plots of each spatial cluster in the areas where they were dominant
kde_full %>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  mutate(density = density / max(density)) %>%
  ggplot(aes(x, y, fill = density)) +
  geom_tile() +
  facet_grid(rows = vars(group)) +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme(aspect.ratio = 0.5)


#Density plots of each trait cluster in the areas where they were dominant

kde_kmeans %>%
  group_by(x, y) %>% #get the dominant group from each quadrant
  slice_max(density) %>%
  ungroup() %>%
  mutate(density = density / max(density)) %>% #Scale the values between 0 and 1
  ggplot(aes(x, y, fill = density)) +
  geom_tile() +
  facet_grid(rows = vars(kmeans)) +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme(aspect.ratio = 0.5)

kde_dat = 
  kde_full %>%
  filter(census==7,d_cutoff==20)%>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  inner_join(
    kde_kmeans %>%
      group_by(x, y) %>%
      slice_max(density) %>%
      ungroup() %>%
      rename(density_kmeans = density),
    by = c('x', 'y')
  ) 




#Create a discrepancy table between spatial clusters and trait clusters across quadrats
c_table = 
  kde_dat %>%
  select(group, kmeans) %>%
  count(group, kmeans) %>%
  pivot_wider(names_from = kmeans, values_from = n)

c_table = as.matrix(c_table[, -1])
rcompanion::cramerV(c_table,ci=TRUE,type="norm")



############################################################################################
#7. Leaf elemental concentration traits analysis

#Perform PCAs and k-means analysis
cluster_data = read.csv("bci_clustering_analysis_consistent.csv")

kde_full = read.csv('bci_kde_full.csv')

splist = unique(bci$sp)

bci_dat = 
  bci %>%
  inner_join(
    bci %>%
      group_by(sp) %>%
      summarize(baldeck = quantile(dbh, 0.56), .groups = 'drop'),
    by = 'sp'
  ) %>%
  filter(dbh > baldeck)

foo_st = 
  trait_data_st %>%
  rename(sp = sp6) %>%
  mutate(sp = tolower(sp)) %>%
  filter(LIGHT == 'SUN') %>%
  select(-LIGHT) %>%
  pivot_longer(AL:N, names_to = 'nutrient', values_to = 'conc') %>%
  group_by(sp, nutrient) %>%
  summarize(conc = mean(conc), .groups = 'drop')

# remove outliers
bar = 
  foo_st %>%
  inner_join(
    foo_st %>%
      group_by(nutrient) %>%
      summarize(
        q05 = quantile(conc, 0.05, na.rm = TRUE),
        q95 = quantile(conc, 0.95, na.rm = TRUE),
        .groups = 'drop'
      ),
    by = 'nutrient'
  ) %>%
  filter(conc < q95) %>%
  select(sp, nutrient, conc) %>%
  pivot_wider(names_from = nutrient, values_from = conc) %>%
  mutate(
    CN = C / N,
    NP = N / P
  )


trdat_st = 
  bar %>%
  inner_join(
    cluster_data %>%
      filter(census == 7) %>%
      unique() %>%
      select(sp, group),
    relationship = "many-to-many",
    by = 'sp'
  )


trdat_long_st = 
  trdat_st %>%
  select(c(AL:NP, group)) %>%
  unique() %>%
  pivot_longer(AL:NP, names_to = 'trait', values_to = 'conc')


cmpr = list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))


#Ignore the warning messages
ggboxplot(
  trdat_long_st,
  x = "group",
  y = "conc",
  color = "group"
) +
  facet_wrap( ~ trait, scales = 'free') +
  stat_compare_means(
    comparisons = cmpr,
    tip.length = 0,
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    ),
    vjusts = -0.1
  ) +
  theme_minimal() +
  # Add a border around each facet
  theme(panel.border = element_rect(fill = NA, colour = "grey40")) +
  # Don't include a legend
  guides(colour = FALSE)



splist = unique(bci$sp)


trait_data_w_st = 
  bar %>%
  filter(sp %in% splist) %>%
  drop_na()


pca_data_st = 
  trait_data_w_st %>%
  select(-sp) %>%
  pcaMethods::pca(
    method = 'ppca',
    scale = 'none',
    center = FALSE
  )

trait_data_pca_st = 
  trait_data_w_st %>%
  mutate(
    PC1 = pca_data_st@scores[, 1],
    PC2 = pca_data_st@scores[, 2]
  )


kmeans_dat_st = 
  bci %>%
  inner_join(trait_data_pca_st, by = 'sp') %>%
  inner_join(
    cluster_data %>%
      filter(census == 7) %>%
      select(sp, group) %>%
      unique(),
    relationship = 'many-to-many',
    by = 'sp'
  )



abuns = 
  kmeans_dat_st %>%
  group_by(sp) %>%
  summarize(abun = n()) %>%
  mutate(abun = log(abun))

trdat_st = 
  trait_data_pca_st %>%
  inner_join(abuns, by = 'sp')

k_st =
  expand_grid(
    groups = 2:10,
    seed = 0:100
  ) %>%
  future_pmap_dfr(
    .f = function(groups, seed) {
      foo =
        trdat_st %>%
        select(AL:ZN)
      wgt = trdat_st$abun
      
      if (seed > 0) {
        set.seed(seed)
        foo = apply(foo, 2, sample)
      }
      wss =
        sum(
          kmeansW(
            foo,
            centers = groups,
            nstart = 100,
            weight = wgt
          )$withinss
        )
      
      return(
        tibble(
          groups = groups,
          seed = seed,
          wss = wss,
          lwss = log(wss)
        )
      )
    },
    .options = furrr_options(seed = TRUE)
  )

plot_gap =
  k_st %>%
  filter(seed == 0) %>%
  inner_join(
    k_st %>%
      filter(seed > 0) %>%
      group_by(groups) %>%
      summarize(
        null = mean(lwss),
        .groups = 'drop'),
    by = 'groups'
  ) %>%
  mutate(gap = null - lwss) %>%
  ggplot(aes(groups, gap)) +
  geom_line() +
  geom_point() +
  xlab('number of groups') +
  theme(aspect.ratio = 1)

k_st %>%
  filter(seed == 0) %>%
  inner_join(
    k_traits %>%
      filter(seed > 0) %>%
      group_by(groups) %>%
      summarize(
        null = mean(wss),
        lnull = mean(lwss),
        sdlnull = sd(lwss) * (sqrt(1 + (1 / 50))),
        .groups = 'drop'
      ),
    by = 'groups'
  ) %>%
  mutate(gap = lnull - lwss) %>%
  mutate(
    diffgap = c(NA, diff(gap)),
    res = diffgap - sdlnull
  ) %>%
  ggplot(aes(groups, res)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab('number of groups') +
  theme(aspect.ratio = 1)


#Optimal clusters = 3

kdash = 
  kmeansW(
    trdat_st %>% select(AL:ZN),
    centers = 3,
    weight = trdat$abun,
    nstart = 100
  )

trdat_st %<>% 
  mutate(kmeans = kdash$cluster) %>% 
  mutate(kmeans = as.factor(kmeans))

PC_cluster = 
  trdat_st %>%
  ggplot(aes(PC1, PC2, col = kmeans)) + geom_point(aes(size = abun)) +
  xlab('Trait PC1') +
  ylab('Trait PC2') +
  scale_colour_manual(values = cbpalette[c(1, 2, 3)]) +
  theme(aspect.ratio = 1)

trdat_st %>%
  pivot_longer(
    cols = AL:ZN,
    names_to = 'Trait',
    values_to = 'value'
  ) %>%
  group_by(sp, Trait) %>%
  ggplot(aes(Trait, value, fill = kmeans)) +
  geom_boxplot() +
  scale_fill_manual(values = cbpalette) +
  theme(
    axis.title = element_text(size = 12, face = 'bold'),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = 'bold')
  ) +
  ylab('Values')

trdat_st %>%
  pivot_longer(
    cols = AL:ZN,
    names_to = 'Element',
    values_to = 'Concentration'
  ) %>%
  ggplot(aes(kmeans, Concentration, fill = Element)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme(
    axis.title = element_text(size = 12, face = 'bold'),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = 'bold')
  )

kmeans_dat_clustered_st =
  kmeans_dat_st %>%
  inner_join(
    trdat_st %>%
      select(sp, kmeans),
    by = 'sp'
  )


#Ignore the warning messages 'data contain duplicated points'.
kde_kmeans_st =
  kmeans_dat_clustered_st %>%
  select(gx, gy, kmeans) %>%
  group_by(kmeans) %>%
  summarize(
    density =
      KernelDensityEstimation(gx = gx, gy = gy, Lx = 1000),
    .groups = 'drop'
  )

kde_kmeans_st = 
  bind_cols(
    kde_kmeans_st$kmeans,
    kde_kmeans_st$density
  )

names(kde_kmeans_st) = c('kmeans', 'x', 'y', 'density')

kde_kmeans_st %>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  ggplot(aes(x, y, fill = kmeans)) +
  geom_tile() +
  theme(aspect.ratio = 0.5) +
  scale_fill_manual(values = cbpalette[c(1, 2, 3)])



#Density plots of each spatial cluster in the areas where they were dominant
kde_full %>%
  filter(census == 7, d_cutoff == 20) %>%
  select(group, x, y, density)%>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  mutate(density = density / max(density)) %>%
  ggplot(aes(x, y, fill = density)) +
  geom_tile() +
  facet_grid(rows = vars(group)) +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme(aspect.ratio = 0.5)

#Density plots of each trait cluster in the areas where they were dominant
kde_kmeans_st %>%
  group_by(x, y) %>% #get the dominant group from each quadrant
  slice_max(density) %>%
  ungroup() %>%
  mutate(density = density / max(density)) %>% #Scale the values between 0 and 1
  ggplot(aes(x, y, fill = density)) +
  geom_tile() +
  facet_grid(rows = vars(kmeans)) +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme(aspect.ratio = 0.5)

kde_dat_st = 
  kde_full %>%
  filter(census == 7, d_cutoff == 20) %>%
  select(group, x, y, density)%>%
  group_by(x, y) %>%
  slice_max(density) %>%
  ungroup() %>%
  inner_join(
    kde_kmeans_st %>%
      group_by(x, y) %>%
      slice_max(density) %>%
      ungroup() %>%
      rename(density_kmeans = density),
    by = c('x', 'y')
  ) 

c_table = 
  kde_dat_st %>%
  select(group, kmeans) %>%
  count(group, kmeans) %>%
  pivot_wider(names_from = kmeans, values_from = n)

c_table = as.matrix(c_table[, -1])
rcompanion::cramerV(c_table,ci=TRUE,type="norm")




###############################################################################################################
# ========================================= END ===============================================================

