#########################################################################
# Spectral analysis: periodogram ()
# Yuting Zou
# August 2020
#########################################################################

#########################################################################
################### Step 1. Load library and data #######################
#########################################################################

# Step a. Install packages and load libraries
library("stats") #spec.pgram()
library('outliers')
library('dplyr')
library('stlplus')
library('ggplot2')
library('Rfast')
library('lattice')
library('latticeExtra')
library('zoo') # as.yearmonth()
library('BBmisc')

# Step b. Set workspace
setwd('D:/college/Minor thesis/DataAnalysis/FastFourierTransform/')

# Step c. Load Spectral Index data
MNDWI <- read.csv('Data/SpectralIndex/MNDWI.csv', header = T)

# Step d. Load water level data
path_WL <- 'D:/college/Minor thesis/DataAnalysis/FastFourierTransform/Data/WaterLevel/'
filenames_WL <- list.files(pattern="*.csv",path = path_WL, full.names=TRUE)
List_WL_new <- lapply(filenames_WL, read.csv)
names(List_WL_new[[9]]) <- c('Date', 'WaterLevel') 
names(List_WL_new) <- filenames_WL
List_WL_new <- List_WL_new[-21]

###############################################################################
################### Step 2. Reorganize water level data #######################
###############################################################################

# Step a. Function to add year for water level data of each study site
List_new = list()
for(i in List_WL_new){
  i = na.omit(i)
  y = substring(i$Date,1,4)
  i$Year <- y
  i$Date<- as.Date(i$Date, "%Y/%m/%d")
  # Export <- write.csv(i,"D:/college/Minor thesis/DataAnalysis/FastFourierTransform/Data/AWEI.csv", row.names = TRUE)
  List_new[[length(List_new)+1]] <- i
  
}
List_name <- c( "AssadLake","AtaturkLake","BalatonLake","BeysehirLake","ChiemseeLake",
                "ConstanceLake","DoiraniLake","DukanLake","ErcekLake","IseoLake","KarakayaReservoir","KastoriaLake",
                "KoroniaLake","MosulLake","OzerosLake","PamvotidaLake","PuenteNuevoReservoir","QadisiyahLake",
                "QarunLake","RazazzaLake","SevanLake","ThartharLake","TshchikskoyeLake","TsimlyanskReservoir",
                "UrmiaLake","VanLake","VegoritidaLake","VolviLakes","YesaLake")
names(List_new) <- List_name

####################################################################################################
####### Step 3: Reorganize water level and spectral index data for Iseo lake (8-day water level) ###
####################################################################################################

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

#Step o: Plot spectral index seasonal component
SI_Season_df <- data.frame(SI_Seasonal,WL_Seasonal,Iseo_DataFrame$Date)
names(SI_Season_df) <- c('SI_Fit','WL_Fit','Date')
d1 <- as.Date(SI_Season_df$Date)
d2 <- as.numeric(SI_Season_df$SI_Fit)
d3 <- as.numeric(SI_Season_df$WL_Fit)
obj1 <- xyplot(d2 ~ d1,
                mgp = c(3, 1, 0),
                col = c("red", "black"),
                xlab=list("Date", fontsize = 13,fontfamily="serif",font = 2,cex = 1.8),
                ylab = list("Seasonal component of 
spectral index time series",cex = 1.8, fontsize = 13,fontfamily="serif",font = 2, col = 'red'),
                ylab.right = list("Seasonal component of 
water level time series",cex = 1.8,fontsize = 13,fontfamily="serif",font = 2, col = 'blue',rot=270),
                par.settings = simpleTheme(col = 1),
                type = c("b","g"),
                layout.widths=list(ylab.axis.padding = 12.5),
                lty=1,
                lwd=2.5,
                pch = 8,
                scales=list(x=list(rot=45,tick.number=15,
                                  cex=1.2,fontfamily="serif",font = 2),
                           y=list(tick.number=5,fontfamily="serif",rot = 90,cex=10))
                 
  )
  
obj2 <- xyplot(d3 ~ d1,type = "o",col="blue",
               lty=1, lwd = 2.5,pch = 16,
               scales=list(x=list(rot=45,tick.number=15,
                                  cex=1.2,fontfamily="serif",font = 2),
                           y=list(tick.number=5,
                                  rot = 90,fontfamily="serif", cex=1.2))
 )
trellis.par.set(axis.text=list(cex=1.3, fontfamily = "serif",font = 2))
final_plot <- doubleYScale(obj1, obj2)

# Step p. Calculate frequency of time series using FFT

# For Spectral Index
specOut_SI<-spec.pgram(SI_Seasonal,taper=0,log='no')
maxSpecFreq_SI<-specOut_SI$freq[which.max(specOut_SI$spec)]
abline(v=maxSpecFreq_SI,lty=2,col='red')  # lty 指定线条类型
period<-1/maxSpecFreq_SI
period  # 48
# Create a vector of periods to label on the graph, units are in years
yrs.period <- rev(c(1/6, 1/5, 1/4, 1/3, 0.5, 1))
yrs.labels <- rev(c("1/6", "1/5", "1/4", "1/3", "1/2", "1"))
#Convert annual period to annual freq, and then to monthly freq
yrs.freqs <- 1/yrs.period * 1/46  
spec_df_SI <- data.frame(specOut_SI$spec,specOut_SI$freq)
names(spec_df_SI) <- c('spec','freq')
spec_df_SI$period <- 1/spec_df_SI$freq
# Plot the periodogram for NDPI 
PLOT_SI <- ggplot(data = subset(spec_df_SI[1:3,]), aes(x = freq, y = spec),lwd = 1.2,col = 'red') + 
  # scale_x_continuous("Period (years)", breaks = yrs.freqs, labels = yrs.labels) + 
  # scale_y_continuous('Spectrum') +
  theme(
    # legend.title=element_text(size=14, family="serif",face="bold"),
    text=element_text(size=20,family="serif",face="bold")) +
  geom_bar(stat="identity")

# For Water Level 
specOut_WL<-spec.pgram(WL_Seasonal,taper=0,log='no')
maxSpecFreq_WL<-specOut_WL$freq[which.max(specOut_WL$spec)]
abline(v=maxSpecFreq_WL,lty=2,col='red')  # lty 指定线条类型
period<-1/maxSpecFreq_WL
period  # 22.73684
# Create a vector of periods to label on the graph, units are in years
yrs.period <- rev(c( 1/6, 1/5, 1/4, 1/3, 0.5, 1))
yrs.labels <- rev(c("1/6", "1/5", "1/4", "1/3", "1/2", "1"))
#Convert annual period to annual freq, and then to monthly freq
yrs.freqs <- 1/yrs.period * 1/46  
spec_df_WL <- data.frame(specOut_WL$spec,specOut_WL$freq)
names(spec_df_WL) <- c('spec','freq')
spec_df_WL$period <- 1/spec_df_WL$freq
# Plot the periodogram for NDPI 
PLOT_WL <- ggplot(data = subset(spec_df_WL)) + 
  geom_line(aes(x = freq, y = spec),lwd = 1.2,col = 'red') + 
  scale_x_continuous("Period (years)", breaks = yrs.freqs, labels = yrs.labels) + 
  scale_y_continuous('Spectrum') +
  theme(
    # legend.title=element_text(size=14, family="serif",face="bold"),
    text=element_text(size=20,family="serif",face="bold"))

# # Find top three highest value
# Index_SI <- order(specOut_SI$spec, decreasing = T)[1:3]
# Index_WL <- order(specOut_WL$spec, decreasing = T)[1:3]
# SpecAmp_SI <- specOut_SI$spec[Index_SI]
# SpecAmp_WL <- specOut_WL$spec[Index_WL]
# SpecFreq_SI <- specOut_SI$freq[Index_SI]
# SpecFreq_WL <- specOut_WL$freq[Index_WL]
# Period_SI <- 1/SpecFreq_SI
# Period_WL <- 1/SpecFreq_WL

#######################################################################
##### Step 3: Calculate FFT for other lakes (Monthly frequency) #######
#######################################################################

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
    
  # Step a. Create a new dataframe which contains full monthly date from 2000 to 2017
  Date_min <- min(as.Date(as.yearmon(SI_Mon_List[[i]]$Month)))
  Date_max <- max(as.Date(as.yearmon(SI_Mon_List[[i]]$Month)))
  Date_full <- seq(Date_min, Date_max, by="month")
  Date_full_DF <- data.frame(Date_full)
  
  # Step b. Add a month column for full time series dataframe
  Date_full_DF$Month <- substring(Date_full_DF$Date_full,1,7)
  names(Date_full_DF) <- c('Date','Month')
    
  # Step c. Merge full time series list and the original time series
  Full_Monthly_ts <- merge(Date_full_DF, SI_Mon_List[[i]], by = "Month",
                           all.x = TRUE)
  # Step d. Adjust the new dataframe (adjust year column)
  Full_Monthly_ts$Year <- substring(Full_Monthly_ts$Month,1,4)
  
  # Step e. Add the result to a list
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

# Step 5. Caluculate Peason correlation coefficients for all lakes
Seasonal_OtherLakes <- list()
for(i in 1:29) {
  # Step a: Normalized dataframe
  Iseo_Norm <- normalize(OtherLakes_SIWL[[i]]) 
  
  # Step b: convert dataframe into time series data
  TS_SI <- ts(Iseo_Norm[,3], frequency=12, start=c(min(OtherLakes_SIWL[[i]]$Year),1))
  TS_WL <- ts(Iseo_Norm$Lake_WL, frequency=12, start=c(min(OtherLakes_SIWL[[i]]$Year),1))
  
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

# Step 6. Calculate FFT results for the seasonal component of all lakes
FFT_OtherLakes <- list()
Specgram_SI <- list()
Specgram_WL <- list()
Spec_plot_SI <- list()
Spec_plot_WL <- list()
for(i in 1:29) {
  
  # For Spectral Index
  specOut_SI<-spec.pgram(Seasonal_OtherLakes[[1]]$SI_Seasonal,taper=0,log='no')
  maxSpecFreq_SI<-specOut_SI$freq[which.max(specOut_SI$spec)]
  abline(v=maxSpecFreq_SI,lty=2,col='red')  # lty 指定线条类型
  
  # Create a vector of periods to label on the graph, units are in years
  yrs.period <- rev(c(1/6, 1/5, 1/4, 1/3, 0.5, 1, 3, 5, 10))
  yrs.labels <- rev(c("1/6", "1/5", "1/4", "1/3", "1/2", "1", "3", "5", "10"))
  #Convert annual period to annual freq, and then to monthly freq
  yrs.freqs <- 1/yrs.period * 1/12  
  spec_df_SI <- data.frame(specOut_SI$spec,specOut_SI$freq)
  names(spec_df_SI) <- c('spec','freq')
  spec_df_SI$period <- 1/spec_df_SI$freq
  
  # Plot the periodogram for NDPI 
  PLOT_SI <- ggplot(data = subset(spec_df_SI[1:3,])) + 
    geom_line(aes(x = freq, y = spec),lwd = 1.2,col = 'red') + 
    scale_x_continuous("Period (years)", breaks = yrs.freqs, labels = yrs.labels) + 
    scale_y_continuous('Spectrum') +
    theme(
      # legend.title=element_text(size=14, family="serif",face="bold"),
      text=element_text(size=20,family="serif",face="bold"))
  
 # For Water Level
  specOut_WL<-spec.pgram(Seasonal_OtherLakes[[i]]$WL_Seasonal,taper=0,log='no')
  maxSpecFreq_WL<-specOut_WL$freq[which.max(specOut_WL$spec)]
  abline(v=maxSpecFreq_WL,lty=2,col='red')  # lty 指定线条类型
  spec_df_WL <- data.frame(specOut_WL$spec,specOut_WL$freq)
  names(spec_df_WL) <- c('spec','freq')
  spec_df_WL$period <- 1/spec_df_WL$freq
  # Plot the periodogram for water level data
  PLOT_WL <- ggplot(data = subset(spec_df_WL)) +
             geom_line(aes(x = freq, y = spec),lwd = 1.2,col = 'red') +
             scale_x_continuous("Period (years)", breaks = yrs.freqs, labels = yrs.labels) +
            scale_y_continuous('Spectrum') +
             theme(
            # legend.title=element_text(size=14, family="serif",face="bold"),
            text=element_text(size=20,family="serif",face="bold"))


  # Find top three highest value
  Index_SI <- order(specOut_SI$spec, decreasing = T)[1:3]
  Index_WL <- order(specOut_WL$spec, decreasing = T)[1:3]
  SpecFreq_SI <- specOut_SI$freq[Index_SI]
  SpecFreq_WL <- specOut_WL$freq[Index_WL]
  Spec_SI <- specOut_SI$spec[Index_SI]
  Spec_WL <- specOut_WL$spec[Index_WL]
  period_SI <- 1/SpecFreq_SI
  period_WL <- 1/SpecFreq_WL
  
  # Combine all results into a dataframe
  Period_all <- data.frame(period_SI,period_WL,Spec_SI,Spec_WL)
  names(Period_all) <- c('Period_SI','Period_WL','Amp_SI','Amp_WL')
  
  #Step5. Add results to the new list
  FFT_OtherLakes[[length(FFT_OtherLakes)+1]] <- Period_all
  Specgram_SI[[length(Specgram_SI)+1]] <- specOut_SI
  Specgram_WL[[length(Specgram_WL)+1]] <- specOut_WL
  Spec_plot_SI[[length(Spec_plot_SI)+1]] <- PLOT_SI
  Spec_plot_WL[[length(Spec_plot_WL)+1]] <- PLOT_WL
  
}

FFT_df <- bind_rows(FFT_OtherLakes)
Name_df <- data.frame(List_name)
Full_df <- cbind(Name_df,FFT_df)

# Export the results into csv file
write.csv(Full_df,"D:/college/Minor thesis/DataAnalysis/FastFourierTransform/Result/0106/FFT_List_0106.csv", row.names = TRUE)

# Save all plots for water level periodogram of all lakes
for(i in 1:length(Spec_plot_SI)) {
  Name <-  paste0(List_name[i],"_SI.png")
  path_save <- 'D:/college/Minor thesis/DataAnalysis/FastFourierTransform/Result/0106'
  plots <- Spec_plot_SI[[i]]
  ggsave(filename=Name,plot=plots,dpi = 400, path = path_save)
}

# Save all plots for MNDWI periodogram of all lakes
for(i in 1:length(Spec_plot_WL)) {
  Name <-  paste0(List_name[i],"_WL.png")
  path_save <- 'D:/college/Minor thesis/DataAnalysis/FastFourierTransform/Result/0106'
  plots <- Spec_plot_WL[[i]]
  ggsave(filename=Name,plot=plots,dpi = 400, path = path_save)
}
