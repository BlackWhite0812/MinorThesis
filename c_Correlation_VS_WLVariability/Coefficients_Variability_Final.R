#########################################################################################
# Yuting Zou
# Compare the Pearson correlation coefficients and the varibility of the water level data
#########################################################################################

# Step1: load the library 
library('dplyr')

# Step2: Set workspace
setwd('D:/college/Minor thesis/DataAnalysis/MostSuitableIndex/Data/WaterLevel/')

# Step3: Laod the water level data 
filenames_WL <- list.files(pattern="*.csv", full.names=TRUE)
List_WL_new <- lapply(filenames_WL, read.csv)
names(List_WL_new[[9]]) <- c('Date', 'WaterLevel') 
names(List_WL_new) <- filenames_WL
List_WL_new <- List_WL_new[-21]

# Step4: Reorganize the water level data
List_new = list()
for(i in List_WL_new){
  i = na.omit(i)
  y = substring(i$Date,1,4)
  i$Year <- y
  i$Date<- as.Date(i$Date, "%Y/%m/%d")
  List_new[[length(List_new)+1]] <- i
}
names(List_new) <- names(List_WL_new)

# Step5: Find the max, min and varibility of the water level data
list_Variability <- list()
for(i in List_new) {
  # Caluculate the max, min, mean and the variability of the water level data
  MaxValue <- max(i$WaterLevel)
  MinValue <- min(i$WaterLevel)
  MeanValue <- abs(mean(i$WaterLevel))
  Variability <- (MaxValue-MinValue)/MeanValue
  SD <- sd(i$WaterLevel)
  SD_Mean <- SD/MeanValue
  
  # Write the results into a dataframe 
  Variability_df <- data.frame(MaxValue,MinValue,MeanValue,Variability,SD,SD_Mean)
  list_Variability[[length(list_Variability)+1]] <- Variability_df
}

# Step6: Change the name of the list
List_name <- c( "AssadLake","AtaturkLake","BalatonLake","BeysehirLake","ChiemseeLake",
                "ConstanceLake","DoiraniLake","DukanLake","ErcekLake","IseoLake","KarakayaReservoir","KastoriaLake",
                "KoroniaLake","MosulLake","OzerosLake","PamvotidaLake","PuenteNuevoReservoir","QadisiyahLake",
                "QarunLake","RazazzaLake","SevanLake","ThartharLake","TshchikskoyeLake","TsimlyanskReservoir",
                "UrmiaLake","VanLake","VegoritidaLake","VolviLakes","YesaLake")
names(list_Variability) <- names(List_name)

# Step7: Convert the list to a dataframe
DF_Vari <- bind_rows(list_Variability)
DF_Vari$Lakename <- List_name

# Step8: Export the results
write.csv(DF_Vari,"D:/college/Minor thesis/DataAnalysis/Compare_correlation_WLVariability/Result/Variability.csv", row.names = TRUE)
