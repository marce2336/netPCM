
# HEATMAPS FIGURES 2A-B AND SUPPLEMENTAL FIGURES 1A-B

#
## Note 1: To create the heatmap figures it's necessary to create the file 'Lipid_profiles4heatmaps.xlsx'. 
#          For this run the matlab script named 'formatData4heatmaps.m'.

#
## Note 2: set as working directory the path where the script 'HeatMapGeneration.R' is located:
#          '../Figures_source_data/Figure_2_&_Supplementary_Figure_1'.


###
##### 1). Create heatmap for in-house dataset (Figure 2A):
###


# load lipid dataset:
library(readxl)
lipids_mpi <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "mpi_dataset")
row.names(lipids_mpi) <- lipids_mpi$species
dataMPI_asmatrix <- data.matrix(lipids_mpi)
dataMPI_asmatrix <- dataMPI_asmatrix[,-1]


# Place asterisks on selected cells according to custom threshold
asterisks_mpi <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "scaffold_asterisk_mpi")
row.names(asterisks_mpi) <- asterisks_mpi$species
new_matrix <- data.matrix(asterisks_mpi)
new_matrix <- new_matrix[,-1]
new_matrix <- abs(new_matrix)
new_matrix[new_matrix == 0] <- ""
new_matrix[new_matrix == 1] <- "*"
new_matrix[new_matrix == 2] <- "**"


# create heatmap:
library(pheatmap)
# Rotate names of columns
library(grid)
d <- matrix(rnorm(25), 5, 5)
colnames(d) = paste("bip", 1:5, sep = "")
rownames(d) = paste("blob", 1:5, sep = "")
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}


## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Group treatments by sampling day and differentiate by color
my_sample_col <- data.frame(class = rep(c("concordant", "outlier"), c(6,12)))
row.names(my_sample_col) <- colnames(dataMPI_asmatrix)
my_colour <- list(class = c(concordant = "#00994c", outlier = "#6600CC"))


# Generate Heatmap
pheatmap(dataMPI_asmatrix, scale = "row", cluster_cols = F, cluster_rows = F, display_numbers = new_matrix, annotation_col = my_sample_col, annotation_colors = my_colour, fontsize_row = 8)

# Save as PNG (500 x 1500), or PDF (8 x 12 inches)


#-------------------------------------------------------------------------------------------------------------------

###
##### 2). Create heatmap for Lusk dataset (Figure 2B):
###


# load lipid dataset:
library(readxl)
lipids_lusk <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "lusk_dataset")
row.names(lipids_lusk) <- lipids_lusk$species
dataLusk_asmatrix <- data.matrix(lipids_lusk)
dataLusk_asmatrix <- dataLusk_asmatrix[,-1]


# Place asterisks on selected cells according to custom threshold
asterisks_lusk <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "scaffold_asterisk_lusk")
row.names(asterisks_lusk) <- asterisks_lusk$species
new_matrix <- data.matrix(asterisks_lusk)
new_matrix <- new_matrix[,-1]
new_matrix <- abs(new_matrix)
new_matrix[new_matrix == 0] <- ""
new_matrix[new_matrix == 1] <- "*"
new_matrix[new_matrix == 2] <- "**"
#View(new_matrix)


# create heatmap:
library(pheatmap)
# Rotate names of columns
library(grid)
d <- matrix(rnorm(25), 5, 5)
colnames(d) = paste("bip", 1:5, sep = "")
rownames(d) = paste("blob", 1:5, sep = "")
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}


## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Group treatments by sampling day and differentiate by color
my_sample_col <- data.frame(class = rep(c("concordant", "outlier"), c(9,7)))
row.names(my_sample_col) <- colnames(dataLusk_asmatrix)
my_colour <- list(class = c(concordant = "#00994c", outlier = "#6600CC"))


# Generate Heatmap
pheatmap(dataLusk_asmatrix, scale = "row", cluster_cols = F, cluster_rows = F, display_numbers = new_matrix, annotation_col = my_sample_col, annotation_colors = my_colour, fontsize_row = 7)

# Save as PNG (500 x 1500), or PDF (5 x 16 inches)




#-------------------------------------------------------------------------------------------------------------------



###
##### 3). Create second heatmap for MPI dataset included in Supplemental Figure 1B:
###


# load lipid dataset:
library(readxl)
lipids_mpi <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "second_heatmap_mpi")
row.names(lipids_mpi) <- lipids_mpi$species
dataMPI_asmatrix <- data.matrix(lipids_mpi)
dataMPI_asmatrix <- dataMPI_asmatrix[,-1]


# Place asterisks on selected cells according to custom threshold
asterisks_mpi <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "second_scaffold_mpi")
row.names(asterisks_mpi) <- asterisks_mpi$species
new_matrix <- data.matrix(asterisks_mpi)
new_matrix <- new_matrix[,-1]
new_matrix <- abs(new_matrix)
new_matrix[new_matrix == 0] <- ""
new_matrix[new_matrix == 1] <- ""
new_matrix[new_matrix == 2] <- "**"
new_matrix[new_matrix == 3] <- "*"
#View(new_matrix)



# create heatmap:
library(pheatmap)
# Rotate names of columns
library(grid)
d <- matrix(rnorm(25), 5, 5)
colnames(d) = paste("bip", 1:5, sep = "")
rownames(d) = paste("blob", 1:5, sep = "")
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}


## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Group treatments by sampling day and differentiate by color
my_sample_col <- data.frame(class = rep(c("concordant", "outlier", "unknown"), c(2,2,1)))
row.names(my_sample_col) <- colnames(dataMPI_asmatrix)
my_colour <- list(class = c(concordant = "#00994c", outlier = "#6600CC", unknown = "grey"))


# Generate Heatmap
pheatmap(dataMPI_asmatrix, scale = "row", cluster_cols = F, cluster_rows = F, display_numbers = new_matrix, annotation_col = my_sample_col, annotation_colors = my_colour, fontsize_row = 8)

# Save as PDF (3.95 x 12 inches)







#-------------------------------------------------------------------------------------------------------------------


###
##### 4). Create second heatmap for Lusk dataset included in Supplemental Figure 1A:
###


# load lipid dataset:
library(readxl)
lipids_lusk <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "second_heatmap_lusk")
row.names(lipids_lusk) <- lipids_lusk$species
dataLusk_asmatrix <- data.matrix(lipids_lusk)
dataLusk_asmatrix <- dataLusk_asmatrix[,-1]


# Place asterisks on selected cells according to custom threshold
asterisks_mpi <- read_excel("Lipid_profiles4heatmaps.xlsx", sheet = "second_scaffold_lusk")
row.names(asterisks_mpi) <- asterisks_mpi$species
new_matrix <- data.matrix(asterisks_mpi)
new_matrix <- new_matrix[,-1]
new_matrix <- abs(new_matrix)
new_matrix[new_matrix == 0] <- ""
new_matrix[new_matrix == 1] <- ""
new_matrix[new_matrix == 2] <- "**"
new_matrix[new_matrix == 3] <- "*"
#View(new_matrix)



# create heatmap:
library(pheatmap)
# Rotate names of columns
library(grid)
d <- matrix(rnorm(25), 5, 5)
colnames(d) = paste("bip", 1:5, sep = "")
rownames(d) = paste("blob", 1:5, sep = "")
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}


## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Group treatments by sampling day and differentiate by color
my_sample_col <- data.frame(class = rep(c("concordant", "outlier", "unknown"), c(7,6,4)))
row.names(my_sample_col) <- colnames(dataLusk_asmatrix)
my_colour <- list(class = c(concordant = "#00994c", outlier = "#6600CC", unknown = "grey"))


# Generate Heatmap
pheatmap(dataLusk_asmatrix, scale = "row", cluster_cols = F, cluster_rows = F, display_numbers = new_matrix, annotation_col = my_sample_col, annotation_colors = my_colour, fontsize_row = 8)

# Save as PDF (7.40 x 12 inches)
# Save as PDF (8 x 12 inches)



