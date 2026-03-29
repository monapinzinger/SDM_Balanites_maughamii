#Species Distribution Model for the African coastal tree species "Balanites maughamii"
#Mona Pinzinger
#pinzin02@ads.uni-passau.de


#first of all: install and load helpful packages to download, visualize data and so on: 
install.packages("ggspatial")
install.packages("ggthemes")
library(ggthemes)
library(here)       
library(terra)     
library(tidyverse)  
library(rgbif)      
library(geodata)    
library(ggspatial)  
library(ggplot2)


dir.create(path = "data_SDM") #create folder to store data I will be analyzing 
dir.create(path = "figures")

# in order to get started: define extent 
Ba_mau_ext <- c(17,52,-35,15)     #extent for the coastal area of Eastern and Southern Africa with a little bit of mainlland buffer

# create and check main map 
Africa_map <- crop(x=world(resolution = 2, level = 0), y = Ba_mau_ext) 
plot(Africa_map, main = "Study Area: Eastern and Southern Africa")


# 1. Environmental variables  
# 1a. elevation 
  #download and save in data folder with resolution = 5 as it's good for a continental scale and a managable file size 
elev <- worldclim_global(var = "elev", res = 5, path = "data_SDM")

  # crop to my extent 
elev_crop <- crop(elev, Ba_mau_ext)

  # check with classic colours 
plot(elev_crop, col = terrain.colors(100),
     main = "Elevation in Eastern and Southern Africa")

plot(Africa_map, add = TRUE, border = "black", lwd = 0.8)

# 1b. download climate data and crop to my extent 

#1d. load climate data and crop to map extent 
bioclim_data <- worldclim_global(var = "bio",
                                 res = 5,
                                 path= "data_SDM")
bioclim_data_crop <- crop(bioclim_data, 
                          Ba_mau_ext) 

#check
plot(bioclim_data_crop)

#1e. add soil data as tree mainly prefers sandy soils

soil_data <- soil_world(var = "sand", depth = 5, res = 5,
                        path = "data_SDM")
soil_crop <- crop(soil_data, Ba_mau_ext)

#check 
plot(soil_crop)

#1f. build predictors stack for model --> everything needs to be on same grid 
env_data <- c(bioclim_data_crop, elev_crop, soil_crop)

## check what's different 
res(bioclim_data_crop)
res(elev_crop)
res(soil_crop)

## resample to bioclim grid 
soil_crop_resampled <- resample(soil_crop, bioclim_data_crop)
elev_crop_resampled <- resample(elev_crop, bioclim_data_crop)

## combination of predictors 
env_data <- c(bioclim_data_crop, elev_crop_resampled, soil_crop_resampled)


## check content to see if stack worked 
names(env_data)


#2. Data on presences/absences 
#2a. species occurence data for Balanites maughamii 

gbif_baMau_download <- occ_data(scientificName = "Balanites maughamii",
                                hasCoordinate = TRUE,
                                limit = 1000,
                                year="1970,2000")

#1970-2000, weil sich das mit Vorhandensein der worldclim Daten überschneidet, aber brauche mehr als nur human observation, sonst keine Daten 

### !!!!!! fact that i left out only human observation needs to be in discussion 


## extract data
gbif_baMau <- gbif_baMau_download$data
head(gbif_baMau) #shows 6 rows with 108 variables but no coordinate Uncertainty in Meters (probably bc. of no Human observation)

## only chose relevant columns 
gbif_baMau <- gbif_baMau[,c("key",
                            "decimalLatitude",
                            "decimalLongitude",
                            "occurrenceStatus")]

##check 
summary(gbif_baMau)
## check map with presences 

plot(Africa_map,
     axes = TRUE, 
     col = "grey95",
     main = "Occurences in Africa")
points(x = gbif_baMau$decimalLongitude, 
       y = gbif_baMau$decimalLatitude, 
       col = "orange", 
       pch = 16, 
       cex = 0.5)


#2b. Create absences 

background <- spatSample(x = env_data,
                         size = nrow(gbif_baMau),
                         values = FALSE,           
                         na.rm = TRUE,             
                         xy = TRUE)  
##check 
head(background)

## visualisation with both presences and absences 
plot(Africa_map,
     axes = TRUE, 
     col = "grey95",
     main = "Presences and Absences in East & South Africa")

points(background,
       col = "darkgrey",
       pch = 1,
       cex = 0.75)

points(x = gbif_baMau$decimalLongitude, 
       y = gbif_baMau$decimalLatitude, 
       col = "orange", 
       pch = 16, 
       cex = 0.5)

## put absences and presences in one dataframe 
### coordinates presences 

baMau_presences <- gbif_baMau[, c("decimalLongitude", "decimalLatitude")]
colnames(baMau_presences) <- c("longitude", "latitude")

### add column
baMau_presences$pa <- 1

### convert to dataframe
baMau_absences <- as.data.frame(background)
colnames(baMau_absences) <- c("longitude", "latitude")

### add column for absences
baMau_absences$pa <- 0

### putting it together
baMau_PA <- rbind(baMau_presences, baMau_absences)

#2c. extract environmental data for testpoints 

environment_df <- terra::extract(x = env_data,
                                 y = baMau_PA[, c("longitude", "latitude")],
                                 ID = FALSE)

## put environmental data in one dataframe with species data
baMau_complete_df <- cbind(baMau_PA, environment_df)

##check
head(baMau_complete_df)


#3. Model building

##to avoid overfitting: GLM with only selected variables
###clean all variable names
names(baMau_complete_df) <- make.names(names(baMau_complete_df))

glm_baMau <- glm(pa ~ wc2.1_5m_bio_1+ wc2.1_5m_bio_5+wc2.1_5m_bio_6+wc2.1_5m_bio_12+wc2.1_5m_bio_15+wc2.1_5m_elev+sand_0.5cm,
                 data = baMau_complete_df[,-c(1,2)],
                 family = binomial())
## second GLM for the climate only model for future prediction 
glm_baMau2 <- glm(pa ~ wc2.1_5m_bio_1+ wc2.1_5m_bio_5+wc2.1_5m_bio_6+wc2.1_5m_bio_12+wc2.1_5m_bio_15,
                  data = baMau_complete_df[,-c(1,2)],
                  family = binomial())

#4. prediction of possible distribution
## 4a. spatial prediction of GLM1
names(env_data) <- make.names(names(env_data))
head(env_data)
predict_baMau <- predict(env_data, glm_baMau, type = "response")

## visualize 
plot(predict_baMau, main = "Predicted Current Habitat Suitability") ## note for me: potential distribution matches with current presences points

##spatial prediction of GLM2 
predict_baMau2 <- predict(bioclim_data_crop, glm_baMau2, type = "response")

## visualize 
plot(predict_baMau2, main = "Predicted Current Habitat Suitability Climate Only")

## 4b. download climate forecast data

forecast_data <- cmip6_world(model = "MPI-ESM1-2-HR",
                             ssp = "245",
                             time = "2061-2080",
                             var = "bioc",
                             res = 5)


names(forecast_data) <- names(bioclim_data_crop)


forecast_data <- crop(x = forecast_data, y = Ba_mau_ext)

## 4c. prediction where species will occur in the future 
forecast_presence <- predict(forecast_data, glm_baMau2, type = "response")

##visualization
plot(Africa_map, 
     axes = TRUE, 
     col = "grey95")

##model-probability
plot(forecast_presence, add = TRUE)

##with Africa borders
plot(Africa_map, add = TRUE, border = "grey5")

##and actual occurences 
points(x = gbif_baMau$decimalLongitude, 
       y = gbif_baMau$decimalLatitude, 
       col = "orange", 
       pch = "+", 
       cex = 0.75)

#5. final SDM map 
## transform raster to dataframe

predict_baMau_df <- as.data.frame(predict_baMau2, xy = TRUE)

ggplot() +
  geom_raster(data = predict_baMau_df, 
              aes(x = x, y = y, 
                  fill = lyr1)) +
  annotation_spatial(Africa_map, fill = NA, colour = "black")+
  scale_fill_viridis_c(name = "Probability") +
  labs(x = "Latitude", 
       y = "Longitude") +
  ggtitle("Realised Niche for Balanites maughamii") +
  theme_bw()+
  coord_sf(expand = FALSE) +  # keeps north arrow outside if you adjust margins
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0, "cm"),   # moves arrow outside
    pad_y = unit(2.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(
    location = "tr",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(8.5, "cm"),
  )


