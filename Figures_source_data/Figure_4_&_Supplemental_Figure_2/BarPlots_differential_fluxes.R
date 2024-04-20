#
# BAR PLOTS FIGURE 4 AND SUPPLEMENTAL FIGURE 2

# This script builds bar plots that summarize the results of the number of reactions with
# differential fluxes for Arabidopsis T-DNA lines that were classified as "full rescue", that is 
# that the mean of the fluxes for WT and T-DNA lines were equal. 

#
## Note: set as working directory the path where the script 'BarPlots_differential_fluxes.R' is located:
# '../Figures_source_data/Figure_4_&_Supplemental_Figure_2'.

# get path were source data was saved:
pthw <- getwd( )
pthw <- gsub("Figures_source_data/Figure_4_&_Supplemental_Figure_2", "SamplingResults/Ttest", pthw)
setwd(pthw)


library(ggplot2)
library(viridis)
library(readxl)



###
##### 1). Figure 4:
###

# Create figure without grid, for proportions including all metabolite classes:

differential_fluxes <- read_excel("barplots_differential_fluxes.xlsx", sheet = "merged")

c <- ggplot(differential_fluxes) + 
  geom_bar(aes(x = metabolite_class, y = proportion, fill = range),
           position = "stack",
           stat = "identity") +
  facet_grid(~ locus, switch = "x") +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white")) +
  scale_fill_viridis(discrete = T) +
  xlab("") +
  ylab("Proportion of reactions with differential fluxes grouped by
metabolite categories") +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))


c
# size print figure: 12.21 x 9
#------------------------------------------------------



###
##### 2). Supplemental Figure 2:
###

# Create figure without grid, for proportions including just lipid-related subsystems:

differential_fluxes <- read_excel("barplots_differential_fluxes.xlsx", sheet = "full_lipids_pruned")

c <- ggplot(differential_fluxes) + 
  geom_bar(aes(x = metabolite_class, y = proportion, fill = range),
           position = "stack",
           stat = "identity") +
  facet_grid(~ locus, switch = "x") +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white")) +
  scale_fill_viridis(discrete = T) +
  xlab("") +
  ylab("Proportion of lipid-related reactions with differential fluxes grouped by
subsystems") +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))


c
# size print figure: 16.23 x 9




