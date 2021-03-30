#########################################################################
# Harmonic analysis
# Yuting Zou
# August 2020
#########################################################################

############### Step 1. Install and load library and data ##############
# Step a. Install packages and load libraries
install.packages('locfit')
install.packages('devtools')
install.packages('extrafont')
install.packages('BBmisc')
install.packages('xts')
install.packages('survMisc')
loadfonts(device = "win")
library('outliers')
library('dplyr') # na_if
library('stlplus')
library('zoo') # as.yearmonth()
library('ggplot2')
library('Rfast')
library('lattice')
library('stlplus')
library('latticeExtra')
library('ggpubr')
library('BBmisc')

# Step b. Set workspace
setwd('D:/college/Minor thesis/DataAnalysis/Harmonic/')

# Step c. Load Spectral Index data
MNDWI <- read.csv('Data/SpectralIndex/MNDWI.csv', header = T)

# Step d. Load water level data
path_WL <- 'D:/college/Minor thesis/DataAnalysis/Harmonic/Data/WaterLevel/'
filenames_WL <- list.files(pattern="*.csv",path = path_WL, full.names=TRUE)
List_WL_new <- lapply(filenames_WL, read.csv)
names(List_WL_new[[9]]) <- c('Date', 'WaterLevel') 
names(List_WL_new) <- filenames_WL
List_WL_new <- List_WL_new[-21]

############### Step2: Reorganize water level data ###############
# Step a. Function to add year for water level data of each study site
List_new = list()
for(i in List_WL_new){
  i = na.omit(i)
  y = substring(i$Date,1,4)
  i$Year <- y
  i$Date<- as.Date(i$Date, "%Y/%m/%d")
  List_new[[length(List_new)+1]] <- i
  
}
List_name <- c( "AssadLake","AtaturkLake","BalatonLake","BeysehirLake","ChiemseeLake",
                "ConstanceLake","DoiraniLake","DukanLake","ErcekLake","IseoLake","KarakayaReservoir","KastoriaLake",
                "KoroniaLake","MosulLake","OzerosLake","PamvotidaLake","PuenteNuevoReservoir","QadisiyahLake",
                "QarunLake","RazazzaLake","SevanLake","ThartharLake","TshchikskoyeLake","TsimlyanskReservoir",
                "UrmiaLake","VanLake","VegoritidaLake","VolviLakes","YesaLake")
names(List_new) <- List_name

##################################################################################################################################
#################### Step 3: Reorganize water level and spectral index data for Iseo lake (8-day water level) ####################
##################################################################################################################################
# Step a. Check outliers
outliers_check <- outlier(MNDWI$IseoLake)

# Step b. remove the rows containing the outliers
outliers_rm <- na_if(as.numeric(MNDWI$IseoLake),outliers_check)
MNDWI$IseoLake <- outliers_rm

# Step c. extract data by common year between water level and spectral index
MNDWI$Date <- as.Date(MNDWI$Date, "%Y-%m-%d")
MNDWI$Iseo_WL <- List_new[[10]]$WaterLevel[match(MNDWI$Date,List_new[[10]]$Date)]
Iseo_common <- intersect(MNDWI$Year, List_new[[10]]$Year)
Iseo_F1 <- MNDWI %>% filter(Year %in% Iseo_common)

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

# Step l: Normalized dataframe
Iseo_Norm <- normalize(Iseo_DataFrame) 

# Step m: convert dataframe into time series data
SI_TS <- ts(Iseo_Norm$IseoLake, frequency=46, start=c(min(Iseo_DataFrame$Year),1))
WL_TS <- ts(Iseo_Norm$Iseo_WL, frequency=46, start=c(min(Iseo_DataFrame$Year),1))


# Step n. using stl to detrend the time series data into seasonal, trend components
SI_STL <- stlplus(SI_TS, s.window = "periodic")
SI_Seasonal = SI_STL$data[,2]
WL_STL <- stlplus(WL_TS, s.window = "periodic")
WL_Seasonal = WL_STL$data[,2]

# Step o. Calculate date difference
date1 <- strptime(Iseo_DataFrame$Date, format="%Y-%m-%d")
date2 <- strptime("01.01.1970", format="%d.%m.%Y")
Dis_Date <- difftime(date1,date2,units="days")
t <- as.numeric(Dis_Date/365)

# Step p. 3th order harmonic model
co <- cos(2 * pi * t)
si <- sin(2 * pi * t)
co2 <- cos(2 * pi * t * 2)
si2 <- sin(2 * pi * t * 2)
co3 <- cos(2 * pi * t * 3)
si3 <- sin(2 * pi * t * 3)

# Step q. Convert ts data into numeric
Harm_SI <- as.numeric(SI_Seasonal)
Harm_WL <- as.numeric(WL_Seasonal)

# Step r. fit the 3-order seasonal model using linear regression
fitm3_SI<- lm(Harm_SI ~ t + co + si + co2 + si2 + co3 + si3) 
fitm3_WL<- lm(Harm_WL ~ t + co + si + co2 + si2 + co3 + si3) 

# Step s. Calculate phase and amplitude metrics for order 3
phase1_or3_SI <- atan2((fitm3_SI$coefficients[3]),fitm3_SI$coefficients[4])
phase2_or3_SI <- atan2((fitm3_SI$coefficients[5]),fitm3_SI$coefficients[6])
phase3_or3_SI <- atan2((fitm3_SI$coefficients[7]),fitm3_SI$coefficients[8])
amplitude1_or3_SI <- sqrt((fitm3_SI$coefficients[3])^2 + (fitm3_SI$coefficients[4])^2)
amplitude2_or3_SI <- sqrt((fitm3_SI$coefficients[5])^2 + (fitm3_SI$coefficients[6])^2)
amplitude3_or3_SI <- sqrt((fitm3_SI$coefficients[7])^2 + (fitm3_SI$coefficients[8])^2)
phase1_or3_WL <- atan2((fitm3_WL$coefficients[3]),fitm3_WL$coefficients[4])
phase2_or3_WL <- atan2((fitm3_WL$coefficients[5]),fitm3_WL$coefficients[6])
phase3_or3_WL <- atan2((fitm3_WL$coefficients[7]),fitm3_WL$coefficients[8])
amplitude1_or3_WL <- sqrt((fitm3_WL$coefficients[3])^2 + (fitm3_WL$coefficients[4])^2)
amplitude2_or3_WL <- sqrt((fitm3_WL$coefficients[5])^2 + (fitm3_WL$coefficients[6])^2)
amplitude3_or3_WL <- sqrt((fitm3_WL$coefficients[7])^2 + (fitm3_WL$coefficients[8])^2)

combine_or3 <- cbind(phase1_or3_SI, phase2_or3_SI, phase3_or3_SI,
                     amplitude1_or3_SI, amplitude2_or3_SI,amplitude3_or3_SI,
                     phase1_or3_WL, phase2_or3_WL, phase3_or3_WL,
                     amplitude1_or3_WL, amplitude2_or3_WL,amplitude3_or3_WL)
names(combine_or3) <- c('phase1_SI','phase2_SI','phase3_SI','amlitude1_SI','amlitude2_SI','amlitude3_SI',
                        'phase1_WL','phase2_WL','phase3_WL','amlitude1_WL','amlitude2_WL','amlitude3_WL')

####################################################################################################
#################### Step 4: Calculate harmonic metrics for other lakes (Monthly frequency) ####################
####################################################################################################
# Step a. Average MODIS observation data into monthly frequency
SI_Mon_List <- list()
for(i in 2:30) {
  
  # Step a. Check outliers
  outliers_check <- outlier(MNDWI[[names(MNDWI[i])]])
  
  # Step b. Find in which rows the outliers are and replace the outliers with NA
  # and update the original dataset
  outliers_rm <- na_if(as.numeric(MNDWI[[names(MNDWI[i])]]),outliers_check)
  MNDWI[[names(MNDWI[i])]] <- outliers_rm
  
  # Step c. add month column and extract each lake
  MNDWI$Month <- substring(MNDWI$Date,1,7)
  Lake <- MNDWI[[names(MNDWI[i])]]
  
  # Step d. aggregate original time series from 8-day frequency into Monthly frequency
  SI_Monthly <- aggregate(Lake~Month,MNDWI,mean)
  names(SI_Monthly)[2] <- paste(names(MNDWI[i]),'Mon', sep="_")
  SI_Monthly$Year <- substring(SI_Monthly$Month,1,4)
  
  # Step e. add the result to a list
  SI_Mon_List[[length(SI_Mon_List)+1]] <- SI_Monthly
}

# step 2. Average water level data into monthly frequency
# Create a new list to save all Pearson correlation coefficients
WL_OtherLakes <- list()
for(i in 1:29) {
  
  #Step a. Find the common year between spectral indieces and water level data
  OtherLake_common <- intersect(SI_Mon_List[[i]]$Year, List_new[[i]]$Year)
  
  #Step b. Select spectral index data and water level data based on common year
  Lake_SI_F1 <- SI_Mon_List[[i]] %>% filter(Year %in% OtherLake_common)
  Lake_WL_F1 <- List_new[[i]] %>% filter(Year %in% OtherLake_common)
  
  #Step c. Add month value to water level dataframe
  Lake_WL_F1$Month <- substring(Lake_WL_F1$Date,1,7)
  
  #Step d. Average water level data to monthly water level data
  OtherLake_WL_Monthly <- aggregate(Lake_WL_F1$WaterLevel~Lake_WL_F1$Month,Lake_WL_F1,mean)
  names(OtherLake_WL_Monthly)[2] <- paste(names(Lake_SI_F1)[2],'WL', sep="_")
  names(OtherLake_WL_Monthly)[1] <- 'Month'
  OtherLake_WL_Monthly$Year <- substring(OtherLake_WL_Monthly$Month,1,4)
  
  # Step e. Add the result to a list
  WL_OtherLakes[[length(WL_OtherLakes)+1]] <- OtherLake_WL_Monthly
  
}

# Step 3. Create full monthly date for time series and merge it with the original one
SI_OtherLakes_full <- list()
for(i in 1:29) {
  
  #Step a. Create a new dataframe which contains full monthly date from 2000 to 2017
  Date_min <- min(as.Date(as.yearmon(SI_Mon_List[[i]]$Month)))
  Date_max <- max(as.Date(as.yearmon(SI_Mon_List[[i]]$Month)))
  Date_full <- seq(Date_min, Date_max, by="month")
  Date_full_DF <- data.frame(Date_full)
  
  #Step b. Add a month column for full time series dataframe
  Date_full_DF$Month <- substring(Date_full_DF$Date_full,1,7)
  names(Date_full_DF) <- c('Date','Month')
  
  #Step c. Merge full time series list and the original time series
  Full_Monthly_ts <- merge(Date_full_DF, SI_Mon_List[[i]], by = "Month",
                           all.x = TRUE)
  #Step d. Adjust the new dataframe (adjust year column)
  Full_Monthly_ts$Year <- substring(Full_Monthly_ts$Month,1,4)
  
  #Step e. Add the result to a list
  SI_OtherLakes_full[[length(SI_OtherLakes_full)+1]] <- Full_Monthly_ts
  
}

# Step 4. Check continuous years and the frequency for each year
OtherLakes_SIWL <- list()
for(i in 1:29) {
  
  # Step a. extract data by common year between water level and spectral index
  SI_OtherLakes_full[[i]]$Lake_WL <- WL_OtherLakes[[i]][,2][match(SI_OtherLakes_full[[i]]$Month,WL_OtherLakes[[i]]$Month)]
  Year_common <- intersect(SI_OtherLakes_full[[i]]$Year, WL_OtherLakes[[i]]$Year)
  Lake_SI_F1 <- SI_OtherLakes_full[[i]] %>% filter(Year %in% Year_common)
  
  # Step b. check whether the frequency of time series is same for each year
  Lake_fre_SI <- aggregate(cbind(count = Year) ~ Year,
                           data = Lake_SI_F1,
                           FUN = function(x){NROW(x)})
  
  # Step c. check the difference for year in the frequency dataframe
  diff_check_year <- diff(as.numeric(Lake_fre_SI$Year))
  
  # Step d. compute the lengths and values of runs of equal values in difference dataframe
  rle_year <- rle(diff_check_year)
  
  # Step e. find the continuous vector with more years
  index_year <- which(rle_year$lengths == max(rle_year$lengths[rle_year$values==1])):(max(rle_year$lengths[rle_year$values==1])+1)
  
  # Step f. create a new dataframe based on the year frequency
  Lake_count_new <- as.numeric(Lake_fre_SI$count[index_year]) 
  Lake_year_new <- as.numeric(Lake_fre_SI$Year[index_year])
  Lake_fre_new <- data.frame(Lake_year_new,Lake_count_new)
  names(Lake_fre_new) <- c('Year', 'count')
  
  # Step g. Check difference for the count component in the update frequency dataframe
  diff_check_count <- diff(as.numeric(Lake_fre_new$count))
  rle_count <- rle(diff_check_count)
  index_count <- which(rle_count$lengths == max(rle_count$lengths[rle_count$values==0])):(max(rle_count$lengths[rle_count$values==0])+1)
  
  # Step h. extract the year from frequency dataframe
  Year_Fre12 <- Lake_fre_new$Year[index_count]
  
  # Step i. filter dataframe based on the frequency of time series
  Lake_DataFrame <- Lake_SI_F1 %>% filter(Year %in% Year_Fre12)
  Lake_DataFrame$Year <- as.numeric(as.character(Lake_DataFrame$Year))
  
  # Step j. Add the result to a list
  OtherLakes_SIWL[[length(OtherLakes_SIWL)+1]] <- Lake_DataFrame
}

# Step 5. decompose the original time series for all lakes
Seasonal_OtherLakes <- list()
for(i in 1:29) {
  # Step a: Normalize dataframe 
  Norm_SIWL <- normalize(OtherLakes_SIWL[[i]])
  
  # Step b: convert dataframe into time series data
  TS_SI <- ts(Norm_SIWL[,3], frequency=12, start=c(min(OtherLakes_SIWL[[i]]$Year),1))
  TS_WL <- ts(Norm_SIWL$Lake_WL, frequency=12, start=c(min(OtherLakes_SIWL[[i]]$Year),1))
  
  # Step c. using stl to detrend the time series data into seasonal, trend components
  SI_STL <- stlplus(TS_SI, s.window = "periodic")
  SI_Seasonal = SI_STL$data[,2]
  WL_STL <- stlplus(TS_WL, s.window = "periodic")
  WL_Seasonal = WL_STL$data[,2]
  
  # Step d. combine seasonal component of spectral index and water level time series
  Seasonal_df <- data.frame(SI_Seasonal,WL_Seasonal,OtherLakes_SIWL[[i]]$Date)
  names(Seasonal_df) <- c('SI_Seasonal','WL_Seasonal','Date')
  
  Seasonal_OtherLakes[[length(Seasonal_OtherLakes)+1]] <- Seasonal_df
}

# Step 6. Create the harmonic model for both sepctral index and water level
# data for all lakes
HM_mod_List_SI <- list()
HM_mod_List_WL <- list()
Harm_Seasonal_SI <- list()
Harm_Seasonal_WL <- list()
  for(i in 1:29) {
    # Step a. Calculate date difference
    date1 <- strptime(Seasonal_OtherLakes[[i]]$Date, format="%Y-%m-%d")
    date2 <- strptime("01.01.1970", format="%d.%m.%Y")
    Dis_Date <- difftime(date1,date2,units="days")
    t <- as.numeric(Dis_Date/365)

    # Step b. Select lake
    NDPI_harm_WL <- Seasonal_OtherLakes[[i]]$WL_Seasonal
    NDPI_harm_SI <- Seasonal_OtherLakes[[i]]$SI_Seasonal
    
    # Step c. Convert ts data into numeric
    Harm_SI <- as.numeric(NDPI_harm_SI)
    Harm_WL <- as.numeric(NDPI_harm_WL)
    
    # Step d. 3th order harmonic model
    co <- cos(2 * pi * t)
    si <- sin(2 * pi * t)
    co2 <- cos(2 * pi * t * 2)
    si2 <- sin(2 * pi * t * 2)
    co3 <- cos(2 * pi * t * 3)
    si3 <- sin(2 * pi * t * 3)

    # Step e. fit the seasonal model using linear regression
    fitm3_SI <- lm(Harm_SI ~ t + co + si + co2 + si2 + co3 + si3)
    fitm3_WL <- lm(Harm_WL ~ t + co + si + co2 + si2 + co3 + si3)
     
    # Step f. Add to list
    HM_mod_List_SI[[length(HM_mod_List_SI)+1]] <- fitm3_SI
    HM_mod_List_WL[[length(HM_mod_List_WL)+1]] <- fitm3_WL
    Harm_Seasonal_SI[[length(Harm_Seasonal_SI)+1]] <- Harm_SI
    Harm_Seasonal_WL[[length(Harm_Seasonal_WL)+1]] <- Harm_WL
  }    

# Step 7. Caluculate harmonic metrics from sepctral index for all lakes
HM_Metrics_List <- list()
  for(i in 1:29) {
    # Step a. Calculate phase and amplitude metrics from spectral index time series for order 3
    phase1_or3_SI <- atan2((HM_mod_List_SI[[i]]$coefficients[3]),HM_mod_List_SI[[i]]$coefficients[4])
    phase2_or3_SI <- atan2((HM_mod_List_SI[[i]]$coefficients[5]),HM_mod_List_SI[[i]]$coefficients[6])
    phase3_or3_SI <- atan2((HM_mod_List_SI[[i]]$coefficients[7]),HM_mod_List_SI[[i]]$coefficients[8])
    amplitude1_or3_SI <- sqrt((HM_mod_List_SI[[i]]$coefficients[3])^2 + (HM_mod_List_SI[[i]]$coefficients[4])^2)
    amplitude2_or3_SI <- sqrt((HM_mod_List_SI[[i]]$coefficients[5])^2 + (HM_mod_List_SI[[i]]$coefficients[6])^2)
    amplitude3_or3_SI <- sqrt((HM_mod_List_SI[[i]]$coefficients[7])^2 + (HM_mod_List_SI[[i]]$coefficients[8])^2)
    
    # Step b. Calculate phase and amplitude metrics from water level time seriesfor order 3
    phase1_or3_WL <- atan2((HM_mod_List_WL[[i]]$coefficients[3]),HM_mod_List_WL[[i]]$coefficients[4])
    phase2_or3_WL <- atan2((HM_mod_List_WL[[i]]$coefficients[5]),HM_mod_List_WL[[i]]$coefficients[6])
    phase3_or3_WL <- atan2((HM_mod_List_WL[[i]]$coefficients[7]),HM_mod_List_WL[[i]]$coefficients[8])
    amplitude1_or3_WL <- sqrt((HM_mod_List_WL[[i]]$coefficients[3])^2 + (HM_mod_List_WL[[i]]$coefficients[4])^2)
    amplitude2_or3_WL <- sqrt((HM_mod_List_WL[[i]]$coefficients[5])^2 + (HM_mod_List_WL[[i]]$coefficients[6])^2)
    amplitude3_or3_WL <- sqrt((HM_mod_List_WL[[i]]$coefficients[7])^2 + (HM_mod_List_WL[[i]]$coefficients[8])^2)
    combine_or3 <- cbind(phase1_or3_SI, phase2_or3_SI, phase3_or3_SI, 
                         amplitude1_or3_SI, amplitude2_or3_SI,amplitude3_or3_SI,
                         phase1_or3_WL, phase2_or3_WL, phase3_or3_WL, 
                         amplitude1_or3_WL, amplitude2_or3_WL,amplitude3_or3_WL)
    names(combine_or3) <- c('phase1_SI','phase2_SI','phase3_SI','amlitude1_SI','amlitude2_SI','amlitude3_SI',
                            'phase1_WL','phase2_WL','phase3_WL','amlitude1_WL','amlitude2_WL','amlitude3_WL')
    
    # Step c. Add to list
    HM_Metrics_List[[length(HM_Metrics_List)+1]] <- combine_or3
  }

# Name for each element in the list
List_name <- c( "AssadLake","AtaturkLake","BalatonLake","BeysehirLake","ChiemseeLake",
                "ConstanceLake","DoiraniLake","DukanLake","ErcekLake","IseoLake","KarakayaReservoir","KastoriaLake",
                "KoroniaLake","MosulLake","OzerosLake","PamvotidaLake","PuenteNuevoReservoir","QadisiyahLake",
                "QarunLake","RazazzaLake","SevanLake","ThartharLake","TshchikskoyeLake","TsimlyanskReservoir",
                "UrmiaLake","VanLake","VegoritidaLake","VolviLakes","YesaLake")
names(HM_Metrics_List) <- List_name

# Combine all dataframe of all lists into one dataframe and export it
lapply(HM_Metrics_List, function(x) write.table( data.frame(x), 'D:/college/Minor thesis/DataAnalysis/Harmonic/Result/HarmonicMetric0106_try.csv'  , append= T, sep=',' ))

######################################################################
#################### Step 5: Plot for other lakes ####################
######################################################################

# Plot for the fitted curve of both NDPI and water level data
for(i in 1:29){
  # Step a. Calculate date difference
  date1 <- strptime(Seasonal_OtherLakes[[i]]$Date, format="%Y-%m-%d")
  date2 <- strptime("01.01.1970", format="%d.%m.%Y")
  Dis_Date <- difftime(date1,date2,units="days")
  t <- as.numeric(Dis_Date/365)
  
  Period <- c(1:12)
  
  # Step b. Find the boundary of y axis
  ymin <- min(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])-0.5
  ymax <- max(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])+0.5
  
  # Step c. Plot for fitted curve
  setwd('D:/college/Minor thesis/DataAnalysis/Harmonic/Result/Figure/Figure0106')
  Name <-  paste0(List_name[i],"_Fitted.png")
  png(Name,width = 2900,height = 1500,res = 300)
  par(family = "serif", cex.axis = 1.4, cex.lab = 1.4)
  plot(Period,HM_mod_List_SI[[i]]$fitted.values[1:12], col="green",type = 'o', lty=1,lwd=2.5,ylab = 'Amplitude',xlab = 'Period (Month intervals)', ylim = c(ymin,ymax))
  lines(Period,HM_mod_List_WL[[i]]$fitted.values[1:12], col="red",type = 'o', lty=1,lwd=2.5)
  legend("bottomright", cex = 1,lwd = c(2,2),xpd = TRUE, bty = "n", inset = c(-0.05,-0.3), col = c('green','red'), legend=c("Fitted curve for MNDWI","Fitted curve for water level"), lty=c(1,1))
  dev.off()

}

#######################################################
# Plot for the 1st term of both MNDWI and water level data
for(i in 1:29){
  # Step a. Calculate date difference
  date1 <- strptime(Seasonal_OtherLakes[[i]]$Date, format="%Y-%m-%d")
  date2 <- strptime("01.01.1970", format="%d.%m.%Y")
  Dis_Date <- difftime(date1,date2,units="days")
  t <- as.numeric(Dis_Date/365)
  
  # Step b. 3th order harmonic model
  co <- cos(2 * pi * t)
  si <- sin(2 * pi * t)
  co2 <- cos(2 * pi * t * 2)
  si2 <- sin(2 * pi * t * 2)
  co3 <- cos(2 * pi * t * 3)
  si3 <- sin(2 * pi * t * 3)
  
  SI_Term1 <- lm(Harm_Seasonal_SI[[i]] ~ t + co + si)
  SI_Term2 <- lm(Harm_Seasonal_SI[[i]] ~  co2 + si2)
  SI_Term3 <- lm(Harm_Seasonal_SI[[i]] ~  co3 + si3)
  WL_Term1 <- lm(Harm_Seasonal_WL[[i]] ~ t + co + si)
  WL_Term2 <- lm(Harm_Seasonal_WL[[i]] ~  co2 + si2)
  WL_Term3 <- lm(Harm_Seasonal_WL[[i]] ~  co3 + si3)
  Period <- c(1:12)
  
  # Step c. Find the boundary of y axis
  ymin <- min(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])-0.5
  ymax <- max(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])+0.5
  
  # Step d. Plot for first harmonic curve
  setwd('D:/college/Minor thesis/DataAnalysis/Harmonic/Figure/Final_Plot')
  Name <-  paste0(List_name[i],"_1st.png")
  png(Name,width = 2900,height = 1500,res = 300)
  par(family = "serif", cex.axis = 1.4, cex.lab = 1.4)
  plot(Period,SI_Term1$fitted.values[1:12], col="green",type = 'o', lty=1,lwd=2.5,ylab = 'Amplitude',xlab = 'Period (Month intervals)', ylim = c(ymin,ymax))
  lines(Period,WL_Term1$fitted.values[1:12], col="red",type = 'o', lty=1,lwd=2.5)
  legend("bottomright", cex = 1,lwd = c(2,2), xpd = TRUE,bty = "n", inset = c(-0.05,-0.3), col = c('green','red'), legend=c("1st term of harmonic model for MNDWI","1st term of harmonic model for water level"), lty=c(1,1))
  
  # Step e. Close the file
  dev.off()
  
}

#######################################################
# Plot for the 2nd term of both NDPI and water level data
for(i in 1:29){
  # Step a. Calculate date difference
  date1 <- strptime(Seasonal_OtherLakes[[i]]$Date, format="%Y-%m-%d")
  date2 <- strptime("01.01.1970", format="%d.%m.%Y")
  Dis_Date <- difftime(date1,date2,units="days")
  t <- as.numeric(Dis_Date/365)
  
  # Step b. 3th order harmonic model
  co <- cos(2 * pi * t)
  si <- sin(2 * pi * t)
  co2 <- cos(2 * pi * t * 2)
  si2 <- sin(2 * pi * t * 2)
  co3 <- cos(2 * pi * t * 3)
  si3 <- sin(2 * pi * t * 3)
  
  SI_Term1 <- lm(Harm_Seasonal_SI[[i]] ~ t + co + si)
  SI_Term2 <- lm(Harm_Seasonal_SI[[i]] ~  co2 + si2)
  SI_Term3 <- lm(Harm_Seasonal_SI[[i]] ~  co3 + si3)
  WL_Term1 <- lm(Harm_Seasonal_WL[[i]] ~ t + co + si)
  WL_Term2 <- lm(Harm_Seasonal_WL[[i]] ~  co2 + si2)
  WL_Term3 <- lm(Harm_Seasonal_WL[[i]] ~  co3 + si3)
  Period <- c(1:12)
  
  # Step c. Find the boundary of y axis
  ymin <- min(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])-0.5
  ymax <- max(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])+0.5
  
  # Step d. Plot for second harmonic curve
  setwd('D:/college/Minor thesis/DataAnalysis/Harmonic/Figure/Final_Plot')
  Name <-  paste0(List_name[i],"_2nd.png")
  png(Name,width = 2900,height = 1500,res = 300)
  par(family = "serif", cex.axis = 1.4, cex.lab = 1.4)
  plot(Period,SI_Term2$fitted.values[1:12], col="green",type = 'o', lty=1,lwd=2.5,ylab = 'Amplitude',xlab = 'Period (Month intervals)', ylim = c(ymin,ymax))
  lines(Period,WL_Term2$fitted.values[1:12], col="red",type = 'o', lty=1,lwd=2.5)
  legend("bottomright", cex = 1,lwd = c(2,2), xpd = TRUE,bty = "n", inset = c(-0.05,-0.3),col = c('green','red'), legend=c("2nd term of harmonic model for MNDWI","2nd term of harmonic model for water level"), lty=c(1,1))
  
  # Step e. Close the file
  dev.off()
  
}
#######################################################
# Plot for the 3rd term of both NDPI and water level data
for(i in 1:29){
  # Step a. Calculate date difference
  date1 <- strptime(Seasonal_OtherLakes[[i]]$Date, format="%Y-%m-%d")
  date2 <- strptime("01.01.1970", format="%d.%m.%Y")
  Dis_Date <- difftime(date1,date2,units="days")
  t <- as.numeric(Dis_Date/365)
  
  # Step b. 3th order harmonic model
  co <- cos(2 * pi * t)
  si <- sin(2 * pi * t)
  co2 <- cos(2 * pi * t * 2)
  si2 <- sin(2 * pi * t * 2)
  co3 <- cos(2 * pi * t * 3)
  si3 <- sin(2 * pi * t * 3)
  
  SI_Term1 <- lm(Harm_Seasonal_SI[[i]] ~ t + co + si)
  SI_Term2 <- lm(Harm_Seasonal_SI[[i]] ~  co2 + si2)
  SI_Term3 <- lm(Harm_Seasonal_SI[[i]] ~  co3 + si3)
  WL_Term1 <- lm(Harm_Seasonal_WL[[i]] ~ t + co + si)
  WL_Term2 <- lm(Harm_Seasonal_WL[[i]] ~  co2 + si2)
  WL_Term3 <- lm(Harm_Seasonal_WL[[i]] ~  co3 + si3)
  Period <- c(1:12)
  
  # Step c. Find the boundary of y axis
  ymin <- min(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])-0.5
  ymax <- max(HM_mod_List_SI[[i]]$fitted.values[1:12],HM_mod_List_WL[[i]]$fitted.values[1:12])+0.5
  
  # Step d. Plot for the third harmonic curve
  setwd('D:/college/Minor thesis/DataAnalysis/Harmonic/Figure/Final_Plot')
  Name <-  paste0(List_name[i],"_3rd.png")
  png(Name,width = 2900,height = 1500,res = 300)
  par(family = "serif", cex.axis = 1.4, cex.lab = 1.4)
  plot(Period,SI_Term3$fitted.values[1:12], col="green",type = 'o', lty=1,lwd=2.5,ylab = 'Amplitude',xlab = 'Period (Month intervals)', ylim = c(ymin,ymax))
  lines(Period,WL_Term3$fitted.values[1:12], col="red",type = 'o', lty=1,lwd=2.5)
  legend("bottomright", cex = 1,lwd = c(2,2), xpd = TRUE,bty = "n", inset = c(-0.05,-0.3), col = c('green','red'), legend=c("3rd term of harmonic model for MNDWI","3rd term of harmonic model for water level"), lty=c(1,1))
  
  # Step e. Close the file
  dev.off()
  
}
######################################