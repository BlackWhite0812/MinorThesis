#########################################################################
# Find the most suitable MODIS spectral index to detect water dynamics
# Answer research question 1:
 #
# Yuting Zou
# August 2020
#########################################################################

################################## Step1: Load library and data ####################################
# Install packages
install.packages('purr')
install.packages('gridExtra')
install.packages('stlplus')

# Load libraries
library('dplyr') #Table merge
library('stats') # Correlation
library('zoo') # Interporate time series data
library('forecast') # Auto.arima() function
library('ggplot2')
library('BBmisc') #Normalize
loadfonts(device = "win")
library('ggpubr') # ggarange() plot
library('grDevices')
library('stlplus')
library('graphics') # outliers
library('datasets')
library('outliers')
library('plotrix')
library('lattice')
library('latticeExtra')

# Set workspace
setwd('D:/college/Minor thesis/DataAnalysis/MostSuitableIndex/')

# Load Spectral Index data
AWEI29 <- read.csv('Data/SpectralIndex_GEE/AWEI.csv', header = T)
AWEI_add <- read.csv('Data/SpectralIndex_GEE/AWEI_Added.csv', header = T)
MNDWI29 <- read.csv('Data/SpectralIndex_GEE/MNDWI.csv', header = T)
MNDWI_add <- read.csv('Data/SpectralIndex_GEE/MNDWI_Added.csv', header = T)
LSWI29 <- read.csv('Data/SpectralIndex_GEE/LSWI.csv', header = T)
LSWI_add <- read.csv('Data/SpectralIndex_GEE/LSWI_Added.csv', header = T)
CWI29 <- read.csv('Data/SpectralIndex_GEE/CWI.csv', header = T)
CWI_add <- read.csv('Data/SpectralIndex_GEE/CWI_Added.csv', header = T)
TCBI29 <- read.csv('Data/SpectralIndex_GEE/TCBI.csv', header = T)
TCBI_add <- read.csv('Data/SpectralIndex_GEE/TCBI_Added.csv', header = T)
NDMI29 <- read.csv('Data/SpectralIndex_GEE/NDMI.csv', header = T)
NDMI_add <- read.csv('Data/SpectralIndex_GEE/NDMI_Added.csv', header = T)

# Load all water level data of 30 study sites
setwd('D:/college/Minor thesis/DataAnalysis/MostSuitableIndex/Data/WaterLevel/')
filenames_WL <- list.files(pattern="*.csv", full.names=TRUE)
List_WL_new <- lapply(filenames_WL, read.csv)
names(List_WL_new[[9]]) <- c('Date', 'WaterLevel') 
names(List_WL_new) <- filenames_WL
List_WL_new <- List_WL_new[-21]

####################################################################################################
############### Step 2: reorganize original spectral index time series #############################
Index.Org <- function(index29, index_add) {
  # Delete unused columns
  drop <- c("imageID",".geo","ResiaLake","NeuchatelLake","KremenchukReservoir")
  index29 <- index29[,!(names(index29) %in% drop)]
  index_add <- index_add[,!(names(index_add) %in% drop)]
  
  # Adjust names to make the table readable
  names(index29)[1] <- "Date"
  names(index29)[3] <- "AtaturkLake"
  names(index_add)[1] <- "Date"
  
  # Adjust date format 
  index29$Date <- substr(index29$Date, 0, 10)
  index_add$Date <- substr(index_add$Date, 0, 10)
  index29$Date<- as.Date(index29$Date, "%Y_%m_%d")
  index_add$Date<- as.Date(index_add$Date, "%Y_%m_%d")
  
  # Join the added wetland to all other study sites
  IndexTS <- cbind(index29, index_add)
  IndexTS <- IndexTS[, !duplicated(colnames(IndexTS))]
  
  # Add year column to spectral index
  IndexTS$Year <- substring(IndexTS$Date,1,4)
  IndexTS <- replace(IndexTS, IndexTS == -9999, NA)
  
  # Change column order of data
  order_name <- c( "AssadLake","AtaturkLake","BalatonLake","BeysehirLake","ChiemseeLake",
                   "ConstanceLake","DoiraniLake","DukanLake","ErcekLake","IseoLake","KarakayaReservoir","KastoriaLake",
                   "KoroniaLake","MosulLake","OzerosLake","PamvotidaLake","PuenteNuevoReservoir","QadisiyahLake",
                   "QarunLake","RazazzaLake","SevanLake","ThartharLake","TshchikskoyeLake","TsimlyanskReservoir",
                   "UrmiaLake","VanLake","VegoritidaLake","VolviLakes","YesaLake","Date","Year" )
  SI_data <- IndexTS[, order_name]

  return(SI_data)
}
AWEI <- Index.Org(AWEI29, AWEI_add)
MNDWI <- Index.Org(MNDWI29, MNDWI_add)
CWI <- Index.Org(CWI29, CWI_add)
LSWI <- Index.Org(LSWI29, LSWI_add)
TCBI <- Index.Org(TCBI29, TCBI_add)
NDMI <- Index.Org(NDMI29, NDMI_add)

####################################################################################################
######### Step3. check the missing value of each index time series ########
# # Step1: Count the number of observations without NA value
# SI_NA <- function(si){
#   NB_Observations <- sapply(si[2:33], function(x) sum(is.na(x)), USE.NAMES = TRUE)
#   # Step2: Select wetland which the observation has less than 200 NA value
#   SI_F1 <- data.frame(NB_Observations[NB_Observations >= 200])
#   return(SI_F1)
# }
# 
# AWEI_NA <- SI_NA(AWEI)
# LSWI_NA <- SI_NA(LSWI)
# MNDWI_NA <- SI_NA(MNDWI)
# TCBI_NA <- SI_NA(TCBI)
# CWI_NA <- SI_NA(CWI)
# NDPI_NA <- SI_NA(NDPI)
# NDMI_NA <- SI_NA(NDMI)
# ResiaLake has more than 200 NA value, 
# Remove this lake in the following steps

###############################################################################
######### Step 4: Reorganize water level data ######### 
# Function to add year for water level data of each study site
List_new = list()
for(i in List_WL_new){
  i = na.omit(i)
  y = substring(i$Date,1,4)
  i$Year <- y
  i$Date<- as.Date(i$Date, "%Y/%m/%d")
  List_new[[length(List_new)+1]] <- i
}
names(List_new) <- names(List_WL_new)

################################################################################
######### Step 5: Peason correlation for Iseo lake (daily water level) #########
################################################################################

Iseo_PC <- function(SI){
  # Step a. check outliers
  outliers_check <- outlier(SI$AssadLake)
  
  # Step b. remove the rows containing the outliers
  outliers_rm <- na_if(as.numeric(SI$IseoLake),outliers_check)
  SI$IseoLake <- outliers_rm
  
  # Step c. extract data by common year between water level and spectral index
  SI$Iseo_WL <- List_new[[10]]$WaterLevel[match(SI$Date,List_new[[10]]$Date)]
  Iseo_common <- intersect(SI$Year, List_new[[10]]$Year)
  Iseo_F1 <- SI %>% filter(Year %in% Iseo_common)
  
  # Step d. check whether the frequency of time series is same for each year
  Iseo_fre <- aggregate(cbind(count = Year) ~ Year,
                        data = Iseo_F1,
                        FUN = function(x){NROW(x)})
  
  # Step e. check the difference for year in the frequency dataframe
  diff_check_year <- diff(as.numeric(Iseo_fre$Year))
  
  # Step f. compute the lengths and values of runs of equal values in difference dataframe
  rle_year <- rle(diff_check_year)
  
  # Step g. find the continuous vector with more years
  index_year <- which(rle_year$lengths == max(rle_year$lengths[rle_year$values==1])):(max(rle_year$lengths[rle_year$values==1])+1)
  
  # Step h. create a new dataframe based on the year frequency
  Iseo_count_new <- as.numeric(Iseo_fre$count[index_year]) 
  Iseo_year_new <- as.numeric(Iseo_fre$Year[index_year])
  Iseo_fre_new <- data.frame(Iseo_year_new,Iseo_count_new)
  names(Iseo_fre_new) <- c('Year', 'count')
  
  # Step i. Check difference for the count component in the update frequency dataframe
  diff_check_count <- diff(as.numeric(Iseo_fre_new$count))
  rle_count <- rle(diff_check_count)
  index_count <- which(rle_count$lengths == max(rle_count$lengths[rle_count$values==0])):(max(rle_count$lengths[rle_count$values==0])+1)
  
  # Step j. extract the year from frequency dataframe
  Year_Fre46 <- Iseo_fre_new$Year[index_count]
  
  # Step k. filter dataframe based on the frequency of time series
  Iseo_DataFrame <- Iseo_F1 %>% filter(Year %in% Year_Fre46)
  Iseo_DataFrame$Year <- as.numeric(as.character(Iseo_DataFrame$Year))
  
  # Step l: convert dataframe into time series data
  SI_TS <- ts(Iseo_DataFrame$IseoLake, frequency=46, start=c(min(Iseo_DataFrame$Year),1))
  WL_TS <- ts(Iseo_DataFrame$Iseo_WL, frequency=46, start=c(min(Iseo_DataFrame$Year),1))
  
  # Step m. using stl to detrend the time series data into seasonal, trend components
  SI_STL <- stlplus(SI_TS, s.window = "periodic")
  SI_Seasonal = SI_STL$data[,2]
  SI_Trend = SI_STL$data[,3]
  SI_Detrend = SI_TS - SI_Trend
  WL_STL <- stlplus(WL_TS, s.window = "periodic")
  WL_Trend = WL_STL$data[,3]
  WL_Seasonal = WL_STL$data[,2]
  WL_Detrend = WL_TS - WL_Trend
  
  # Step n. calculate pearson correlation coefficients for Iseo lake
  # For detrend time series (remove trend component)
  PearsonCor_Iseo_Detrend <- cor(SI_Detrend,WL_Detrend,use = "complete.obs")
  # For seasonal component
  PearsonCor_Iseo_Seasonal <- cor(SI_Seasonal,WL_Seasonal,use = "complete.obs")
  # For trend component
  PearsonCor_Iseo_Trend <- cor(SI_Trend,WL_Trend,use = "complete.obs")
  
  # Step o. comine all pearson correlation coefficients into one dataframe
  PearsonCor <- data.frame(PearsonCor_Iseo_Detrend, PearsonCor_Iseo_Seasonal,
                           PearsonCor_Iseo_Trend)
  return(PearsonCor)
}
  
PC_Iseo_AWEI <- Iseo_PC(AWEI) # r = -0.440488
PC_Iseo_LSWI <- Iseo_PC(LSWI) # r = 0.5514626
PC_Iseo_TCBI <- Iseo_PC(TCBI) # r = 0.496838
PC_Iseo_CWI <- Iseo_PC(CWI) # r = 0.02404845
PC_Iseo_NDPI <- Iseo_PC(NDPI) # r = -0.4130759
PC_Iseo_MNDWI <- Iseo_PC(MNDWI) # r = 0.4197414
PC_Iseo_NDMI <- Iseo_PC(NDMI) 
Iseo_PC <- rbind(PC_Iseo_AWEI,PC_Iseo_LSWI,PC_Iseo_TCBI,
                 PC_Iseo_CWI,PC_Iseo_NDPI,PC_Iseo_MNDWI,
                 PC_Iseo_NDMI)
# Export the results for Iseo lake (daily) into csv file
write.csv(Iseo_PC,"D:/college/Minor thesis/DataAnalysis/MostSuitableIndex/Result/Iseo_PearsonResults(8-days)_new.csv", row.names = TRUE)

# # Plot for time series for Iseo lake
# Plot_IrregularLake <- function(TS) {
#   d1 <- as.Date(TS$Date)
#   d2 <- as.numeric(TS$IseoLake)
#   d3 <- as.numeric(TS$Iseo_WL)
#   obj1 <- xyplot(d2 ~ d1,
#                  mgp = c(3, 1, 0),
#                  col = c("red", "black"),
#                  xlab=list("Date", fontsize = 16,fontfamily="serif",font = 2),
#                  ylab = list("Spectral Index Time Series", fontsize = 16,fontfamily="serif",font = 2, col = 'red'),
#                  ylab.right = list("Water Level Time Series", fontsize = 16,fontfamily="serif",font = 2, col = 'blue',rot=270),
#                  par.settings = simpleTheme(col = 1),
#                  type = c("b","g"),
#                  layout.widths=list(ylab.axis.padding = 12.5),
#                  lty=1,
#                  lwd=2.5,
#                  pch = 8,
#                  scales=list(x=list(rot=45,tick.number=20,
#                                     cex=0.9,fontfamily="serif",font = 2),
#                              y=list(rot = 90,cex=1, fontfamily="serif"))
# 
#   )
# 
#   obj2 <- xyplot(d3 ~ d1,type = "o",col="blue",
#                  lty=1, lwd = 2.5,pch = 16,
#                  scales=list(x=list(rot=45,tick.number=20,
#                                     cex=0.9,fontfamily="serif",font = 2),
#                              y=list(rot = 90,cex=1, fontfamily="serif"))
#   )
# 
#   trellis.par.set(axis.text=list(cex=1.1, fontfamily = "serif",font = 2))
#   final_plot <- doubleYScale(obj1, obj2)
# 
#   return(final_plot)
# }
# IrregularLake_AWEI <- Plot_IrregularLake(PC_Iseo_AWEI)
# IrregularLake_LSWI <- Plot_IrregularLake(PC_Iseo_LSWI)
# IrregularLake_CWI <- Plot_IrregularLake(PC_Iseo_CWI)
# IrregularLake_NDPI <- Plot_IrregularLake(PC_Iseo_NDPI)
# IrregularLake_MNDWI <- Plot_IrregularLake(PC_Iseo_MNDWI)
# IrregularLake_TCBI <- Plot_IrregularLake(PC_Iseo_TCBI)


########################################################################################
########## Step 6: Peason correlation for other lakes (irregular water level) ######### 
########################################################################################
OtherLakes_Preprocess <- function(SI){
  # Step a. Average MODIS observation data into monthly frequency
  SI_Mon_List <- list()
  for(i in 1:29) {
    
    # Step1. check outliers
    outliers_check <- outlier(SI[[names(SI[i])]])
    
    # Step2. find in which rows the outliers are and replace the outliers with NA
    # and update the original dataset
    outliers_rm <- na_if(as.numeric(SI[[names(SI[i])]]),outliers_check)
    SI[[names(SI[i])]] <- outliers_rm
    
    # Step3. add month column and extract each lake
    SI$Month <- substring(SI$Date,1,7)
    Lake <- SI[[names(SI[i])]]
    
    # Step4. aggregate original time series from 8-day frequency into Monthly frequency
    SI_Monthly <- aggregate(Lake~Month,SI,mean)
    names(SI_Monthly)[2] <- paste(names(SI[i]),'Mon', sep="_")
    SI_Monthly$Year <- substring(SI_Monthly$Month,1,4)
    
    # Step5. add the result to a list
    SI_Mon_List[[length(SI_Mon_List)+1]] <- SI_Monthly
  }

  # step b. Average water level data into monthly frequency
  # Create a new list to save all Pearson correlation coefficients
  WL_OtherLakes <- list()
  for(i in 1:29) {
    
    #Step 1. find the common year between spectral indieces and water level data
    OtherLake_common <- intersect(SI_Mon_List[[i]]$Year, List_new[[i]]$Year)
    
    #Step 2. select spectral index data and water level data based on common year
    Lake_SI_F1 <- SI_Mon_List[[i]] %>% filter(Year %in% OtherLake_common)
    Lake_WL_F1 <- List_new[[i]] %>% filter(Year %in% OtherLake_common)
    
    #Step3. add month value to water level dataframe
    Lake_WL_F1$Month <- substring(Lake_WL_F1$Date,1,7)
    
    #Step4. average water level data to monthly water level data
    OtherLake_WL_Monthly <- aggregate(Lake_WL_F1$WaterLevel~Lake_WL_F1$Month,Lake_WL_F1,mean)
    names(OtherLake_WL_Monthly)[2] <- paste(names(Lake_SI_F1)[2],'WL', sep="_")
    names(OtherLake_WL_Monthly)[1] <- 'Month'
    OtherLake_WL_Monthly$Year <- substring(OtherLake_WL_Monthly$Month,1,4)
    
    #Step5. add the result to a list
    WL_OtherLakes[[length(WL_OtherLakes)+1]] <- OtherLake_WL_Monthly
    
  }
  
  # Step c. Create full monthly date for time series and merge it with the original one
  SI_OtherLakes_full <- list()
  for(i in 1:29) {
    
    #Step1. Create a new dataframe which contains full monthly date from 2000 to 2017
    Date_min <- min(as.Date(as.yearmon(SI_Mon_List[[i]]$Month)))
    Date_max <- max(as.Date(as.yearmon(SI_Mon_List[[i]]$Month)))
    Date_full <- seq(Date_min, Date_max, by="month")
    Date_full_DF <- data.frame(Date_full)
    
    #Step2. Add a month column for full time series dataframe
    Date_full_DF$Month <- substring(Date_full_DF$Date_full,1,7)
    names(Date_full_DF) <- c('Date','Month')
    
    #Step3. Merge full time series list and the original time series
    Full_Monthly_ts <- merge(Date_full_DF, SI_Mon_List[[i]], by = "Month",
                             all.x = TRUE)
    #Step4. Adjust the new dataframe (adjust year column)
    Full_Monthly_ts$Year <- substring(Full_Monthly_ts$Month,1,4)
    
    #Step5. Add the result to a list
    SI_OtherLakes_full[[length(SI_OtherLakes_full)+1]] <- Full_Monthly_ts
    
  }
  
  # Step d. Check continuous years and the frequency for each year
  OtherLakes_SIWL <- list()
  for(i in 1:29) {
    
    # Step1. extract data by common year between water level and spectral index
    SI_OtherLakes_full[[i]]$Lake_WL <- WL_OtherLakes[[i]][,2][match(SI_OtherLakes_full[[i]]$Month,WL_OtherLakes[[i]]$Month)]
    Year_common <- intersect(SI_OtherLakes_full[[i]]$Year, WL_OtherLakes[[i]]$Year)
    Lake_SI_F1 <- SI_OtherLakes_full[[i]] %>% filter(Year %in% Year_common)
    
    # Step2. check whether the frequency of time series is same for each year
    Lake_fre_SI <- aggregate(cbind(count = Year) ~ Year,
                             data = Lake_SI_F1,
                             FUN = function(x){NROW(x)})
    
    # Step3. check the difference for year in the frequency dataframe
    diff_check_year <- diff(as.numeric(Lake_fre_SI$Year))
    
    # Step4. compute the lengths and values of runs of equal values in difference dataframe
    rle_year <- rle(diff_check_year)
    
    # Step5. find the continuous vector with more years
    index_year <- which(rle_year$lengths == max(rle_year$lengths[rle_year$values==1])):(max(rle_year$lengths[rle_year$values==1])+1)
    
    # Step6. create a new dataframe based on the year frequency
    Lake_count_new <- as.numeric(Lake_fre_SI$count[index_year]) 
    Lake_year_new <- as.numeric(Lake_fre_SI$Year[index_year])
    Lake_fre_new <- data.frame(Lake_year_new,Lake_count_new)
    names(Lake_fre_new) <- c('Year', 'count')
    
    # Step7. Check difference for the count component in the update frequency dataframe
    diff_check_count <- diff(as.numeric(Lake_fre_new$count))
    rle_count <- rle(diff_check_count)
    index_count <- which(rle_count$lengths == max(rle_count$lengths[rle_count$values==0])):(max(rle_count$lengths[rle_count$values==0])+1)
    
    # Step8. extract the year from frequency dataframe
    Year_Fre12 <- Lake_fre_new$Year[index_count]
    
    # Step9. filter dataframe based on the frequency of time series
    Lake_DataFrame <- Lake_SI_F1 %>% filter(Year %in% Year_Fre12)
    Lake_DataFrame$Year <- as.numeric(as.character(Lake_DataFrame$Year))
    
    # Step10. Add the result to a list
    OtherLakes_SIWL[[length(OtherLakes_SIWL)+1]] <- Lake_DataFrame
  }
  return(OtherLakes_SIWL)
}
# Output for preprocessed time series data
  Pre_AWEI <- OtherLakes_Preprocess(AWEI)
  Pre_LSWI <- OtherLakes_Preprocess(LSWI)
  Pre_CWI <- OtherLakes_Preprocess(CWI)
  Pre_MNDWI <- OtherLakes_Preprocess(MNDWI)
  Pre_TCBI <- OtherLakes_Preprocess(TCBI)
  Pre_NDMI <- OtherLakes_Preprocess(NDMI)

#################################################################
# Step f. Caluculate Peason correlation coefficients for all lakes
OtherLakes_PC <- function(PreData){
  
  # create empty lists
  PC_OtherLakes <- list()
  SI_Trend_OtherLakes <- list()
  SI_Seasonal_OtherLakes <- list()
  WL_Trend_OtherLakes <- list()
  WL_Seasonal_OtherLakes <- list()
  
  for(i in 1:29) {
    # Step1: convert dataframe into time series data
    TS_SI <- ts(PreData[[i]][,3], frequency=12, start=c(min(PreData[[i]]$Year),1))
    TS_WL <- ts(PreData[[i]]$Lake_WL, frequency=12, start=c(min(PreData[[i]]$Year),1))
    
    # Step2: using stl to detrend the time series data into seasonal, trend components
    SI_STL <- stlplus(TS_SI, s.window = "periodic")
    SI_Seasonal = SI_STL$data[,2]
    SI_Trend = SI_STL$data[,3]
    SI_Detrend = TS_SI - SI_Trend
    WL_STL <- stlplus(TS_WL, s.window = "periodic")
    WL_Trend = WL_STL$data[,3]
    WL_Seasonal = WL_STL$data[,2]
    WL_Detrend = TS_WL - WL_Trend
    seasonal_check <- data.frame(SI_Seasonal,WL_Seasonal)
    
    # Step3. calculate pearson correlation coefficients for Iseo lake
    # For detrend time series (remove trend component)
    PearsonCor_Detrend <- cor(SI_Detrend,WL_Detrend,use = "complete.obs")
    # For seasonal component
    PearsonCor_Seasonal <- cor(SI_Seasonal,WL_Seasonal,use = "complete.obs")
    # For trend component
    PearsonCor_Trend <- cor(SI_Trend,WL_Trend,use = "complete.obs")
    
    # Step4. comine all pearson correlation coefficients into one dataframe
    PearsonCor <- data.frame(PearsonCor_Detrend, PearsonCor_Seasonal,
                             PearsonCor_Trend)
    
    #Step5. Add results to the new list
    PC_OtherLakes[[length(PC_OtherLakes)+1]] <- PearsonCor
    SI_Trend_OtherLakes[[length(SI_Trend_OtherLakes)+1]] <- SI_Trend
    SI_Seasonal_OtherLakes[[length(SI_Seasonal_OtherLakes)+1]] <- SI_Seasonal
    WL_Trend_OtherLakes[[length(WL_Trend_OtherLakes)+1]] <- WL_Trend
    WL_Seasonal_OtherLakes[[length(WL_Seasonal_OtherLakes)+1]] <- WL_Seasonal
    
  }
  
  List_Otherlake <-list(PC_OtherLakes, SI_Trend_OtherLakes,SI_Seasonal_OtherLakes,
              WL_Trend_OtherLakes,WL_Seasonal_OtherLakes)
  

  return(List_Otherlake)
}

Coe_AWEI <- OtherLakes_PC(Pre_AWEI)
Coe_LSWI <- OtherLakes_PC(Pre_LSWI)
Coe_CWI <- OtherLakes_PC(Pre_CWI)
Coe_MNDWI <- OtherLakes_PC(Pre_MNDWI)
Coe_TCBI <- OtherLakes_PC(Pre_TCBI)
Coe_NDMI <- OtherLakes_PC(Pre_NDMI)

# aame for each element in the list
List_name <- c( "AssadLake","AtaturkLake","BalatonLake","BeysehirLake","ChiemseeLake",
                "ConstanceLake","DoiraniLake","DukanLake","ErcekLake","IseoLake","KarakayaReservoir","KastoriaLake",
                "KoroniaLake","MosulLake","OzerosLake","PamvotidaLake","PuenteNuevoReservoir","QadisiyahLake",
                "QarunLake","RazazzaLake","SevanLake","ThartharLake","TshchikskoyeLake","TsimlyanskReservoir",
                "UrmiaLake","VanLake","VegoritidaLake","VolviLakes","YesaLake")
names(Coe_AWEI) <- List_name
names(Coe_LSWI) <- List_name
names(Coe_CWI) <- List_name
names(Coe_NDPI) <- List_name
names(Coe_MNDWI) <- List_name
names(Coe_TCBI) <- List_name
names(Coe_NDMI) <- List_name


# combine all dataframe of all lists into one dataframe
OtherLakes_PC_comb <- function(Coe_SI){
  Coe_df <- bind_rows(Coe_SI)
  Name_df <- data.frame(List_name)
  Full_df <- cbind(Name_df,Coe_df)
  return(Full_df)
}
Pc_df_AWEI <- OtherLakes_PC_comb(Coe_AWEI)
Pc_df_LSWI <- OtherLakes_PC_comb(Coe_LSWI)
Pc_df_CWI <- OtherLakes_PC_comb(Coe_CWI)
Pc_df_NDPI <- OtherLakes_PC_comb(Coe_NDPI)
Pc_df_MNDWI <- OtherLakes_PC_comb(Coe_MNDWI)
Pc_df_TCBI <- OtherLakes_PC_comb(Coe_TCBI)
Pc_df_NDMI <- OtherLakes_PC_comb(Coe_NDMI)

# combine all results for each spectral index and each lake
fULL_DF <- cbind(Pc_df_AWEI,Pc_df_LSWI,Pc_df_CWI,
                 Pc_df_NDPI,Pc_df_MNDWI,Pc_df_TCBI,Pc_df_NDMI)

# Export the results into csv file
# write.csv(fULL_DF,"D:/college/Minor thesis/DataAnalysis/MostSuitableIndex/Result/Intermediate Plot/1108/PearsonCorrelation_NEW_NDMI.csv", row.names = TRUE)

##################################################################################
########################### Step 7: Plot analysis ################################
##################################################################################
# plot the original and the trend of time series data 
# and water level for comparison
# Plot_IrregularLake <- function(TS) {
  d1 <- as.Date(Pre_MNDWI[[23]]$Date)
  d2 <- as.numeric(Pre_MNDWI[[23]]$TshchikskoyeLake)
  d3 <- as.numeric(Pre_MNDWI[[23]]$Lake_WL)
  d4 <- as.numeric(Coe_MNDWI[[2]][[23]])
  d5 <- as.numeric(Coe_MNDWI[[4]][[23]])
  d6 <- as.numeric(Coe_MNDWI[[3]][[23]])
  d7 <- as.numeric(Coe_MNDWI[[5]][[23]])
   
  obj1 <- xyplot(d2 + d4 ~ d1,
                 mgp = c(3, 1, 0),
                 col = c("brown1", "brown3"),
                 xlab=list("Date", fontsize = 16,fontfamily="serif",font = 2.8),
                 ylab = list("Spectral Index Time Series", fontsize = 16,fontfamily="serif",font = 2.8, col = 'red'),
                 ylab.right = list("Water Level Time Series", fontsize = 16,fontfamily="serif",font = 2.8, col = 'blue',rot=270),
                 par.settings = simpleTheme(col = 1),
                 type = c("b","g"),
                 layout.widths=list(ylab.axis.padding = 12.5),
                 lty=1,
                 lwd=2.5,
                 pch = 8,
                 scales=list(x=list(rot=45,tick.number=15,
                                    cex=1.1,fontfamily="serif",font = 2.8),
                             y=list(rot = 90,cex=2, fontfamily="serif"))
                 
  )
  
  obj2 <- xyplot(d3 + d5 ~ d1,
                 type = "o",
                 col=c("deepskyblue1",'deepskyblue4'),
                 lty=1, lwd = 2.5,pch = 16,
                 scales=list(x=list(rot=45,tick.number=15,
                                    cex=1.1,fontfamily="serif",font = 2.8),
                             y=list(rot = 90,cex=2, fontfamily="serif"))
  )
  
  trellis.par.set(axis.text=list(cex=1.1, fontfamily = "serif",font = 2))
  final_plot <- doubleYScale(obj1, obj2)
  
  # return(final_plot)
# }
# IrregularLake_AWEI <- Plot_IrregularLake(Coe_AWEI)
# IrregularLake_LSWI <- Plot_IrregularLake(Coe_LSWI)
# IrregularLake_CWI <- Plot_IrregularLake(Coe_CWI)
# IrregularLake_NDPI <- Plot_IrregularLake(Coe_NDPI)
# IrregularLake_MNDWI <- Plot_IrregularLake(Coe_MNDWI)
# IrregularLake_TCBI <- Plot_IrregularLake(Coe_TCBI)
# IrregularLake_NDMI <- Plot_IrregularLake(Coe_NDMI)


# plot the seasonal component of time series data 
# and water level data for comparison
d1 <- as.Date(Pre_LSWI[[9]]$Date)
d2 <- as.numeric(Pre_LSWI[[9]]$ErcekLake)
d3 <- as.numeric(Pre_LSWI[[9]]$Lake_WL)
d4 <- as.numeric(Coe_LSWI[[2]][[9]])
d5 <- as.numeric(Coe_LSWI[[4]][[9]])
d6 <- as.numeric(Coe_LSWI[[3]][[9]])
d7 <- as.numeric(Coe_LSWI[[5]][[9]])

obj1 <- xyplot(d6 ~ d1,
               mgp = c(3, 1, 0),
               col = c("brown3"),
               xlab=list("Date", fontsize = 16,fontfamily="serif",font = 2.8),
               ylab = list("Spectral Index Time Series", fontsize = 16,fontfamily="serif",font = 2.8, col = 'red'),
               ylab.right = list("Water Level Time Series", fontsize = 16,fontfamily="serif",font = 2.8, col = 'blue',rot=270),
               par.settings = simpleTheme(col = 1),
               type = c("b","g"),
               layout.widths=list(ylab.axis.padding = 12.5),
               lty=1,
               lwd=2.5,
               pch = 8,
               scales=list(x=list(rot=45,tick.number=15,
                                  cex=1.1,fontfamily="serif",font = 2.8),
                           y=list(rot = 90,cex=2, fontfamily="serif"))
               
)

obj2 <- xyplot(d7 ~ d1,
               type = "o",
               col=c('deepskyblue4'),
               lty=1, lwd = 2.5,pch = 16,
               scales=list(x=list(rot=45,tick.number=15,
                                  cex=1.1,fontfamily="serif",font = 2.8),
                           y=list(rot = 90,cex=2, fontfamily="serif"))
)

trellis.par.set(axis.text=list(cex=1.1, fontfamily = "serif",font = 2))
final_plot <- doubleYScale(obj1, obj2)
