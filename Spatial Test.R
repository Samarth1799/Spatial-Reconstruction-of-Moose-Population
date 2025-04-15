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

Plots <- read_sf("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\ne_mn_aerial_plots_2024\\ne_mn_aerial_plots_2024\\2024NEMR.shp")
plot(st_geometry(Plots))
Aerial_map <- st_read("C:\\Users\\Victoria\\OneDrive - University of St. Thomas\\Desktop\\DASC 460\\Data\\Re_ DASC 460 - Moose Data for SPR\\ne_mn_aerial_plots_2024\\ne_mn_aerial_plots_2024\\2024NEMR.shp")

g <- ggplot() + geom_sf(data = Aerial_map)  + geom_sf(data = Zones_map) + geom_sf(MN_map)
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
  
  # Perimeter
  shared_perimeters <- function(sf, nb) {
    perim_matrix <- matrix(0, nrow = nrow(sf), ncol = nrow(sf))
    for (i in seq_along(nb)) {
      for (j in nb[[i]]) {
        if (i < j) {
          shared <- st_intersection(sf[i,], sf[j,])
          if (nrow(shared) > 0) {
            length <- st_length(st_union(shared))
            perim_matrix[i, j] <- length
            perim_matrix[j, i] <- length
          }
        }
      }
    }
    perim_matrix
  }
  
  # Compute shared perimeter matrix once
  perim_mat <- shared_perimeters(sf, nb)
  
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
      
      # Score each candidate by balance + compactness
      scores <- sapply(candidates, function(cand) {
        # Shared perimeter with current cluster members
        cluster_members <- which(cluster_assignment == i)
        shared <- sum(perim_mat[cand, cluster_members])
        # Area balance score (smaller is better)
        area_diff <- abs(cluster_area[i] + V(g)$area[cand] - target_area)
        # Weighted score: you can tune weight
        -shared + 0.0001 * area_diff  # negative shared so we prefer large values
      })
      
      best <- candidates[which.min(scores)]
      
      cluster_assignment[best] <- i
      cluster_area[i] <- cluster_area[i] + V(g)$area[best]
      
      # Add this node's neighbors to queue
      new_neighbors <- neighbors(g, best)
      queue_list[[i]] <- unique(c(queue_list[[i]], new_neighbors))
    }
  }
  
  sf$cluster <- cluster_assignment
  
  
  tm_shape(sf) +
    tm_polygons("cluster", palette = "Dark2")
}
Cluster(Zones, 3, 69)

# Pretty picture
g <- ggplot() +
  geom_sf(data = Aerial_map, aes(fill = "Aerial Area"), color = "black") +
  geom_sf(data = Zones_map, aes(fill = "Zone Area"), color = "black", alpha = 0.7) +
  scale_fill_manual(name = "Map Layers",
                    values = c("Aerial Area" = "lightgreen", "Zone Area" = "lightblue")) +
  theme_minimal() +
  labs(title = "Overlay of Aerial and Moose Zones Maps") +
  labs(x = "Longitude") +
  labs(y = "Latitude")
g

# Putting plots into aerial zones

library(sf)
library(dplyr)


# 3. Compute the intersection between plots and aerial zones.
#    This will return features that have attributes from both layers.
intersections <- st_intersection(Plots, Zones)

# 4. Compute the area of each intersection polygon.
intersections <- intersections %>%
  mutate(intersection_area = st_area(geometry))

# 5. Compute the original area of each plot (from the plots shapefile).
#    Assumes that plots have a unique identifier "plot_id".
plot_areas <- Plots %>%
  mutate(plot_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  dplyr::select(PLOTNUMBER, plot_area)

# 6. Join the original plot areas into the intersection data by the "plot_id" field.
intersections <- intersections %>%
  left_join(plot_areas, by = "PLOTNUMBER")

# 7. Compute the proportion of the plot that falls within the aerial zone.
intersections <- intersections %>%
  mutate(proportion = as.numeric(intersection_area / plot_area))

# 8. Filter intersections where the proportion is at least 50%.
#    This gives the plot-to-zone assignments.
#    Change "aerial_zone_id" to the correct field name in your aerial_zones dataset.
 plot_assignments <- intersections %>%
  filter(proportion >= 0.5) %>%
  dplyr::select(PLOTNUMBER, MOOSE_ZONE)
plot_assignments <- st_drop_geometry(plot_assignments) %>% as.data.frame()

# 9. Join these assignments back to the original plots data.

# Test
# For each plot, determine which aerial zone captures the maximum fraction of its area
plot_assignments <- intersections %>%
  group_by(PLOTNUMBER) %>%
  slice_max(order_by = proportion, with_ties = FALSE) %>%
  ungroup() %>%
  # Use > 0.5 if you want strictly more than 50%; use >= 0.5 if you want to allow exactly 50%.
  filter(proportion > 0.5) %>%  
  dplyr::select(PLOTNUMBER, MOOSE_ZONE)  # Replace aerial_zone_id with your actual column name in the aerial zone data

# Now join these assignments back to the original plots
plots_with_zones <- Plots %>%
  left_join(st_drop_geometry(plot_assignments), by = "PLOTNUMBER")

# 10. (Optional) Check the result by plotting:
plot(st_geometry(plots_with_zones), col = plots_with_zones$MOOSE_ZONE, main = "Plots Assigned to Aerial Zones")




# Brooke's work
### Aerial Zones by Moose Zones

# Calculate intersection between aerial plots and zones
intersection <- st_intersection(Plots, Zones)

# Calculate area of intersection
intersection$int_area <- st_area(intersection)

# Assign ID and area to original Aerial_map BEFORE intersection
Plots <- Plots %>%
  mutate(PLOTNUMBER = row_number(),
         aerial_area = st_area(geometry))

#Do the intersection and keep aerial ID and geometry
intersection <- st_intersection(
  Plots %>% dplyr::select(PLOTNUMBER, aerial_area), 
  Zones
)

#Calculate intersection area
intersection <- intersection %>%
  mutate(int_area = st_area(geometry))

#Calculate overlap ratio (intersection area / original aerial area)
intersection <- intersection %>%
  mutate(overlap_ratio = as.numeric(int_area) / as.numeric(aerial_area))

#Filter: only keep aerials that are at least 50% covered by a zone
over_50 <- intersection %>%
  filter(overlap_ratio >= 0.5)

#Instead of filtering for 50%, just find the top overlap for *every* aerial_id
over_max <- intersection %>%
  group_by(PLOTNUMBER) %>%
  slice_max(overlap_ratio, n = 1) %>%
  ungroup()

# Join back to full aerial shapes
zone_lookup <- over_max %>%
  st_drop_geometry() %>%
  dplyr::select(PLOTNUMBER, MOOSE_ZONE)

aerial_full_labeled <- Plots %>%
  left_join(zone_lookup, by = "PLOTNUMBER")

sum(is.na(plots_with_zones[,28]))


###Graph Aerial Plots by Moose Zones

install.packages("Polychrome") 
library(Polychrome)


# Ensure MOOSE_ZONE is a factor (with NA preserved)
aerial_full_labeled <- aerial_full_labeled %>%
  mutate(MOOSE_ZONE = as.factor(MOOSE_ZONE))

# Get all defined zone levels (excluding NA)
zone_ids <- levels(droplevels(aerial_full_labeled$MOOSE_ZONE))

# Generate 30 distinct colors
set.seed(42)
zone_colors <- Polychrome::createPalette(length(zone_ids), seedcolors = "#000000")

# Name the colors according to zone IDs
names(zone_colors) <- zone_ids

# Add gray for NA
final_colors <- c(zone_colors, "NA" = "white")

ggplot() +
  geom_sf(data = aerial_full_labeled, aes(fill = MOOSE_ZONE), color = "black", size = 0.2) +
  scale_fill_manual(
    values = final_colors,
    na.value = "white",
    name = "Moose Zone"
  ) +
  theme_minimal() +
  labs(
    title = "Aerial Plots by Moose Zone",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8)
  )


    







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



# Redo of plots in zones
# Add an identifier and compute area for each aerial plot.
# Replace "geometry" if your geometry column is named differently.
Plots <- Plots %>%
  mutate(PLOTNUMBER = row_number(),
         aerial_area = st_area(geometry))

# Compute the intersection between aerial plots and zones.
intersection <- st_intersection(Plots %>% select(PLOTNUMBER, aerial_area), Zones)

# Calculate the area of each intersection.
intersection <- intersection %>%
  mutate(int_area = st_area(geometry))

# Calculate the overlap ratio (intersection area divided by original aerial plot area).
intersection <- intersection %>%
  mutate(overlap_ratio = as.numeric(int_area) / as.numeric(aerial_area))

# For each aerial plot, pick the zone with the highest overlap ratio.
over_max <- intersection %>%
  group_by(PLOTNUMBER) %>%
  slice_max(order_by = overlap_ratio, with_ties = FALSE) %>%  # get the maximum overlap per plot
  ungroup()

# If you want to require >50% overlap, you can filter now:
# over_max <- over_max %>% filter(overlap_ratio > 0.5)

# Create a lookup table of aerial plot assignments to moose zones.
# Here we assume that Zones_map has a field named "MOOSE_ZONE".
zone_lookup <- over_max %>%
  st_drop_geometry() %>%
  dplyr::select(PLOTNUMBER, MOOSE_ZONE)

# Join the zone assignment back to the original aerial plot features.
aerial_full_labeled <- Plots %>%
  left_join(zone_lookup, by = "PLOTNUMBER")



