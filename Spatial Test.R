# Spatial packages
install.packages("sf")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("spdep")
install.packages("raster")
install.packages("tidyverse")
install.packages("igraph")
install.packages("tmap")
install.packages("units")
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(spdep)
library(tidyverse)
library(igraph)
library(tmap)
library(units)





MN_sf <- read_sf("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\shp_bdry_state\\shp_bdry_state\\Boundaries_of_Minnesota.shp")
plot(st_geometry(MN_sf))
MN_map <- st_read("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\shp_bdry_state\\shp_bdry_state\\Boundaries_of_Minnesota.shp")

Zones <- read_sf("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\ne_mn_moose_zones_2007\\ne_mn_moose_zones_2007\\ne_mn_moose_zones_2007.shp")
plot(st_geometry(Zones))
Zones_map <- st_read("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\ne_mn_moose_zones_2007\\ne_mn_moose_zones_2007\\ne_mn_moose_zones_2007.shp")

Aerial_sf <- read_sf("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\ne_mn_aerial_plots_2024\\ne_mn_aerial_plots_2024\\2024NEMR.shp")
plot(st_geometry(Aerial_sf))
Aerial_map <- st_read("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\ne_mn_aerial_plots_2024\\ne_mn_aerial_plots_2024\\2024NEMR.shp")

g <- ggplot() + geom_sf(data = MN_map) + geom_sf(data = Zones_map)
g

# WHOLE CLUSTERING FUNCTION
# CLUSTERING
Cluster <- function(sf, k, seed){
  sf$area <- st_area(sf)
  
  # Neighbors using Queen contiguity
  nb <- poly2nb(sf, queen = TRUE)
  adj_list <- nb2listw(nb, style = "B", zero.policy = TRUE)
  
  # Convert to igraph
  edges <- nb2lines(nb, coords = st_coordinates(st_centroid(sf)))
  g <- graph_from_adj_list(nb, mode = "all")
  
  V(g)$name <- as.character(1:length(sf$area))
  V(g)$area <- as.numeric(sf$area)
  
  
  # Total area and target area per cluster
  total_area <- sum(V(g)$area)
  target_area <- total_area / k
  
  # Start with k seed nodes
  set.seed(seed)
  seeds <- sample(V(g), k)
  
  # Initialize cluster assignment
  cluster_assignment <- rep(NA, length(V(g)))
  cluster_area <- rep(0, k)
  queue_list <- vector("list", k)
  
  # Assign seeds
  for (i in 1:k) {
    cluster_assignment[seeds[i]] <- i
    cluster_area[i] <- V(g)$area[seeds[i]]
    queue_list[[i]] <- neighbors(g, seeds[i])
  }
  
  # Greedy expansion loop
  while(any(is.na(cluster_assignment))) {
    for (i in 1:k) {
      if (length(queue_list[[i]]) == 0) next
      
      # Find unassigned neighbors
      candidates <- queue_list[[i]]
      candidates <- candidates[is.na(cluster_assignment[candidates])]
      
      if (length(candidates) == 0) next
      
      # Pick one with lowest impact on balance
      best <- candidates[which.min(abs(cluster_area[i] + V(g)$area[candidates] - target_area))]
      
      cluster_assignment[best] <- i
      cluster_area[i] <- cluster_area[i] + V(g)$area[best]
      
      # Add this node's neighbors to queue
      new_neighbors <- neighbors(g, best)
      queue_list[[i]] <- unique(c(queue_list[[i]], new_neighbors))
    }
  }
  
  sf$cluster <- cluster_assignment
  
  
  tm_shape(sf) +
    tm_polygons("cluster", palette = "Set3")
}
Cluster(Zones, 3, 69)






# junk




# CLUSTERING
Zones$area <- st_area(Zones)

# Neighbors using Queen contiguity
nb <- poly2nb(Zones, queen = TRUE)
adj_list <- nb2listw(nb, style = "B", zero.policy = TRUE)

# Convert to igraph
edges <- nb2lines(nb, coords = st_coordinates(st_centroid(Zones)))
g <- graph_from_adj_list(nb, mode = "all")

V(g)$name <- as.character(1:length(Zones$area))
V(g)$area <- as.numeric(Zones$area)

# Chat code...
# Number of clusters you want
k <- 4

# Total area and target area per cluster
total_area <- sum(V(g)$area)
target_area <- total_area / k

# Start with k seed nodes
set.seed(69)
seeds <- sample(V(g), k)

# Initialize cluster assignment
cluster_assignment <- rep(NA, length(V(g)))
cluster_area <- rep(0, k)
queue_list <- vector("list", k)

# Assign seeds
for (i in 1:k) {
  cluster_assignment[seeds[i]] <- i
  cluster_area[i] <- V(g)$area[seeds[i]]
  queue_list[[i]] <- neighbors(g, seeds[i])
}

# Greedy expansion loop
while(any(is.na(cluster_assignment))) {
  for (i in 1:k) {
    if (length(queue_list[[i]]) == 0) next
    
    # Find unassigned neighbors
    candidates <- queue_list[[i]]
    candidates <- candidates[is.na(cluster_assignment[candidates])]
    
    if (length(candidates) == 0) next
    
    # Pick one with lowest impact on balance
    best <- candidates[which.min(abs(cluster_area[i] + V(g)$area[candidates] - target_area))]
    
    cluster_assignment[best] <- i
    cluster_area[i] <- cluster_area[i] + V(g)$area[best]
    
    # Add this node's neighbors to queue
    new_neighbors <- neighbors(g, best)
    queue_list[[i]] <- unique(c(queue_list[[i]], new_neighbors))
  }
}

Zones$cluster <- cluster_assignment


tm_shape(Zones) +
  tm_polygons("cluster", palette = "Set3") +
  tm_layout(title = "Contiguous Clusters with Balanced Area")



# WHOLE CLUSTERING FUNCTION
# CLUSTERING
Cluster <- function(sf, k, seed){
sf$area <- st_area(sf)

# Neighbors using Queen contiguity
nb <- poly2nb(sf, queen = TRUE)
adj_list <- nb2listw(nb, style = "B", zero.policy = TRUE)

# Convert to igraph
edges <- nb2lines(nb, coords = st_coordinates(st_centroid(sf)))
g <- graph_from_adj_list(nb, mode = "all")

V(g)$name <- as.character(1:length(sf$area))
V(g)$area <- as.numeric(sf$area)


# Total area and target area per cluster
total_area <- sum(V(g)$area)
target_area <- total_area / k

# Start with k seed nodes
set.seed(seed)
seeds <- sample(V(g), k)

# Initialize cluster assignment
cluster_assignment <- rep(NA, length(V(g)))
cluster_area <- rep(0, k)
queue_list <- vector("list", k)

# Assign seeds
for (i in 1:k) {
  cluster_assignment[seeds[i]] <- i
  cluster_area[i] <- V(g)$area[seeds[i]]
  queue_list[[i]] <- neighbors(g, seeds[i])
}

# Greedy expansion loop
while(any(is.na(cluster_assignment))) {
  for (i in 1:k) {
    if (length(queue_list[[i]]) == 0) next
    
    # Find unassigned neighbors
    candidates <- queue_list[[i]]
    candidates <- candidates[is.na(cluster_assignment[candidates])]
    
    if (length(candidates) == 0) next
    
    # Pick one with lowest impact on balance
    best <- candidates[which.min(abs(cluster_area[i] + V(g)$area[candidates] - target_area))]
    
    cluster_assignment[best] <- i
    cluster_area[i] <- cluster_area[i] + V(g)$area[best]
    
    # Add this node's neighbors to queue
    new_neighbors <- neighbors(g, best)
    queue_list[[i]] <- unique(c(queue_list[[i]], new_neighbors))
  }
}

sf$cluster <- cluster_assignment


tm_shape(sf) +
  tm_polygons("cluster", palette = "Set3")
}
Cluster(Zones, 5, 69)


#To make new shapefile
st_write




# Other info if want to try out

# g + coord_sf(ylim = c(5200000, 5400000), xlim = c(9600000,8900000))


P <- st_crop(MN_sf, xmin = -96, xmax = -89)
ggplot() + geom_sf(data = P)


plot(st_geometry(MN_sf))
plot(st_geometry(Zones_sf), add = TRUE)
plot(st_geometry(Aerial_sf), add = TRUE)





