
####################################################################
# OPTION 3: Boxplot using 'ggplot' function: With solid colored boxes
# This part of the script creates a boxplot for the log2 intensities

library(tidyverse)
library(reshape2)
library(readxl)
library(ggplot2)


# logRawIntensities / logMedianNormIntensities / logMeanNormIntensities

dataSet <- read_excel("LipidomicProfiles.xlsx", sheet = "logMedianNormIntensities")
dataSet <- as.data.frame(dataSet)
dataSet <- subset(dataSet, select = -1)
mergedSet <- melt(list(dataSet), id = c("sampleNumber", "Measurement_batch"))


# lock in factor level order

ggplot(mergedSet, aes(sampleNumber, value, fill = Measurement_batch)) + 
  geom_boxplot(show.legend = FALSE) + 
  labs(fill = "Measurement_batch") + 
  scale_color_manual(values = c("darkblue", "darkred")) +
  theme(legend.position = "bottom", 
        #axis.text.x = element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())


################################################################################
################################################################################
#
# This part of the script is to analyze (qualitatively) if the normalization 
# method eliminates the batch effect
#
################################################################################
################################################################################

# Create boxPlot for intensities:

library(tidyverse)
library(reshape2)
library(readxl)
library(ggplot2)



dataSet <- read_excel("intensities_normalized.xlsx", sheet = "medianNormalizedIntensities")
dataSet <- as.data.frame(dataSet)


# Load dplyr 
library('dplyr')
subSet <- dataSet %>% select(-c(1,2,3))

# Calculate log
logData <- log(subSet)
logData <- as.data.frame(logData)


# Add row Names:
rowNames <- dataSet %>% select(c(2,3))
rowNames <- as.data.frame(rowNames)


# Merge dataframes and merge into one single list
mergedData <- cbind(rowNames, logData)
mergedSet <- melt(list(mergedData), id = c("sampleNumber", "Measurement_batch"))



ggplot(mergedSet, aes(sampleNumber, value, fill = Measurement_batch)) + 
  geom_boxplot(show.legend = FALSE) + 
  labs(fill = "Measurement_batch") + 
  scale_color_manual(values = c("darkblue", "darkred")) +
  theme(legend.position = "bottom", 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())



