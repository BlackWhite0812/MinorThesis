##############################################################
# Yuting Zou
# February
# Caluclate the varibility of the water level data
##############################################################

# Step1: install and load libraries
library(ggplot2)

# Step2: Load reorganized file of the variability of the water level data 
#         (ordered by different climate zones)
setwd('D:/college/Minor thesis/DataAnalysis/Compare_correlation_WLVariability/data')
Variability_org <- read.csv('Variability_org.csv', header = T)
setwd('D:/college/Minor thesis/DataAnalysis/Compare_correlation_WLVariability/Result')

# Step3: Scatter plot to show the relationship between the MODIS-derived indicies 
#         and the water level data
for(i in 5:16) {
  
  yValue <- Variability_org[,i]
  plot_name <- names(Variability_org)[i]
  ggplot(Variability_org, aes(x=SD, y=yValue, color=Climate)) + 
    geom_point(size=2.5) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(size = 12,hjust = 0.5, vjust = 2, face="bold", family = "serif",angle = 90),
          axis.text.x = element_text(size = 12,face="bold", family = "serif"),
          legend.position = "none",
          panel.background = element_blank(),
          strip.background = element_rect(colour=NA, fill=NA),
          panel.border = element_rect(fill = NA, color = "black"))
  ggsave(paste0(plot_name,".png"))

}
