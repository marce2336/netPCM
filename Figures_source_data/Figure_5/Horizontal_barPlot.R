
#
# HORIZONTAL BAR PLOTS FIGURE 5
#

#
## Note: set as working directory the path where the script 'Horizontal_barPlot.R' is located:
# '../Figures_source_data/Figure_5'.

# get path were source data was saved:
pthw <- getwd( )
pthw <- gsub("Figures_source_data/Figure_5", "SamplingResults/Ttest", pthw)
setwd(pthw)



# Create figure for FC-transcripts/mean-fluxes:

library(readxl)
library(ggplot2)
library(viridis)


degs <- read_excel("Validated_TDNAs_full_rescue.xlsx", sheet = "collapsed_validated_locus")



c <- ggplot(degs, aes(fill=dataType, y=locus, x=log2_FC)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette = "Set2") +
  xlab("log2 fold change") +
  ylab("locus") +
  theme_classic() +
  geom_vline(xintercept = 0) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank())

c



# size print figure: 9 x 12



##############################################################################

