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
plot(Africa_map, main = "Study Area: south-east Africa")


# 1. Environmental variables that influence the distribution of Balanites maughamii 
# 1a. elevation 
  #download and save in data folder with resolution = 5 as it's good for a continental scale and a managable file size 
elev <- worldclim_global(var = "elev", res = 5, path = "data_SDM")

  # crop to my extent 
elev_crop <- crop(elev, Ba_mau_ext)

  # check with classic colours 
plot(elev_crop, col = terrain.colors(100),
     main = "Study Area: south-east Africa with Elevation")

plot(Africa_map, add = TRUE, border = "black", lwd = 0.8)

# 1b. download climate data and crop to my extent 
bioclim_data <- worldclim_global(var = "bio",
                                 res = 5,
                                 path= "data_SDM")
bioclim_data_crop <- crop(bioclim_data, 
                          Ba_mau_ext) 

  # check
plot(bioclim_data_crop)

# 1c. add soil data as tree mainly prefers sandy soils
  # download from soil_world and crop too extent 
soil_data <- soil_world(var = "sand", depth = 5, res = 5,
                        path = "data_SDM")
soil_crop <- crop(soil_data, Ba_mau_ext)

  # check 
plot(soil_crop)

# 1d. build predictors stack for model because everything needs to be on same grid 

  # check what's different 
res(bioclim_data_crop)
res(elev_crop)
res(soil_crop) #different than the other two 

  # resample to bioclim grid 
soil_crop_resampled <- resample(soil_crop, bioclim_data_crop)


  # combination of predictors 
env_data <- c(bioclim_data_crop, elev_crop, soil_crop_resampled)


  # check content to see if stack worked 
names(env_data)
  # save cropped environmental data in different folder 
dir.create(path = "data_output")
saveRDS(env_data, "data_output/env_data.rds")

#2. Data on presences/absences 
#2a. species occurence data for Balanites maughamii from gbif 

gbif_baMau_download <- occ_data(scientificName = "Balanites maughamii",
                                hasCoordinate = TRUE,
                                limit = 1000,
                                year="1970,2000") # no human observances limitation because then almost no data was available 
                                                  # timeframe because it's the same one as the climate data

 
  # extract data
gbif_baMau <- gbif_baMau_download$data
  # look at data 
head(gbif_baMau) # shows no coordinate Uncertainty in Meters

  # only chose relevant columns 
gbif_baMau <- gbif_baMau[,c("key",
                            "decimalLatitude",
                            "decimalLongitude",
                            "occurrenceStatus")]

  # check 
summary(gbif_baMau)
  # check map with presences 

plot(Africa_map,
     axes = TRUE, 
     col = "grey95",
     main = "Occurences in south-east Africa (GBIF)")
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
  # check 
head(background)

  # visualisation with both presences and absences 
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

  # put absences and presences in one dataframe 
  ## coordinates of the presences 

baMau_presences <- gbif_baMau[, c("decimalLongitude", "decimalLatitude")]
colnames(baMau_presences) <- c("longitude", "latitude")

  ## add column
baMau_presences$pa <- 1

  ## convert to dataframe
baMau_absences <- as.data.frame(background)
colnames(baMau_absences) <- c("longitude", "latitude")

  ## add column for absences
baMau_absences$pa <- 0

  ## putting it together
baMau_PA <- rbind(baMau_presences, baMau_absences)


#2c. extract environmental data for testpoints 

environment_df <- terra::extract(x = env_data,
                                 y = baMau_PA[, c("longitude", "latitude")],
                                 ID = FALSE)

## put environmental data in one dataframe with species data
baMau_complete_df <- cbind(baMau_PA, environment_df)

## check
head(baMau_complete_df)



# 3. Model building: Generalized linear model with the bioclimate data as explanatory variable 

  # to avoid overfitting: GLM 1 with only selected variables (Annual Mean Temperature, Max./Min. Temperature of warmest/coldest month, Annual Precipitation, Precipitation Seasonality, elevation, and sand)
  ## clean all variable names
names(baMau_complete_df) <- make.names(names(baMau_complete_df))

glm_baMau <- glm(pa ~ wc2.1_5m_bio_1+ wc2.1_5m_bio_5+wc2.1_5m_bio_6+wc2.1_5m_bio_12+wc2.1_5m_bio_15+wc2.1_5m_elev+sand_0.5cm,
                 data = baMau_complete_df[,-c(1,2)],
                 family = binomial())
                #family=binomial because the GLM predicts the probability of presence (between 0 and 1) --> binary presence/absence

  # GLM 2 for future prediction with climate only variables because no available data for the future for the other variables 
glm_baMau2 <- glm(pa ~ wc2.1_5m_bio_1+ wc2.1_5m_bio_5+wc2.1_5m_bio_6+wc2.1_5m_bio_12+wc2.1_5m_bio_15,
                  data = baMau_complete_df[,-c(1,2)],
                  family = binomial())



# 4. prediction of possible distribution
  # 4a. spatial prediction of GLM1

names(env_data) <- make.names(names(env_data))
head(env_data)

  ## prediction based on climate, sand, and elevation data 
predict_baMau <- predict(env_data, glm_baMau, type = "response") 

  ## save it to output folder 
saveRDS(env_data, "data_output/predict_baMau.rds")

  # visualize 
plot(predict_baMau, main = "Predicted Current Habitat Suitability") 

  # spatial prediction of GLM2 in order to predict future 
predict_baMau2 <- predict(bioclim_data_crop, glm_baMau2, type = "response")

  ## save it to output folder 
saveRDS(env_data, "data_output/predict_baMau2.rds")

  # visualize 
plot(predict_baMau2, main = "Predicted Current Habitat Suitability Climate Only")


  # 4b. download climate forecast data

forecast_data <- cmip6_world(model = "MPI-ESM1-2-HR",
                             ssp = "245",
                             time = "2061-2080",
                             var = "bioc",
                             res = 5)     #same resolution as other data  

  # use same names 
names(forecast_data) <- names(bioclim_data_crop)

  # crop
forecast_data <- crop(x = forecast_data, y = Ba_mau_ext)


  # 4c. prediction where species will occur in the future 
forecast_presence <- predict(forecast_data, glm_baMau2, type = "response")

  ## save it to output folder 
saveRDS(env_data, "data_output/forecast_presence.rds")

  # visualization step by step 
plot(Africa_map, 
     axes = TRUE, 
     col = "grey95")

  ## model-probability
plot(forecast_presence, add = TRUE)

  ## with Africa borders
plot(Africa_map, add = TRUE, border = "grey")

  
## and actual occurences 
points(x = gbif_baMau$decimalLongitude, 
       y = gbif_baMau$decimalLatitude, 
       col = "orange", 
       pch = "+", 
       cex = 0.75)


#5. Final SDM maps
  ## transform raster to dataframe in order for ggplot to work 

predict_baMau_climate_df <- as.data.frame(predict_baMau2, xy = TRUE) #dataframe for current potential distribution climate only 

predict_baMau_df <- as.data.frame(predict_baMau, xy = TRUE) #dataframe for current potential distribution with whole environmental data 

predict_baMau_future_df <- as.data.frame(forecast_presence, xy = TRUE) #dataframe for future potential distribution under changing climate conditions 


  #5a. final map: Current Habitat Suitability of B. maughamii (Full Environmental Model)

plotSDM_full_env <-
ggplot() +
  geom_raster(data = predict_baMau_df, 
              aes(x = x, y = y, 
                  fill = lyr1)) +
  annotation_spatial(Africa_map, fill = NA, colour = "black")+
  scale_fill_viridis_c(name = "Probability") +
  labs(x = "Longitude", 
       y = "Latitude") +
  ggtitle(expression(bold("Current Habitat Suitability of ") ~ bolditalic(B.~maughamii) ~ bold("(Full Environmental Model)"))) +
  theme_bw()+
  theme(
    plot.title = element_text(
      hjust = 0.5,        
      size = 11))+
  coord_sf(expand = FALSE) +  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0, "cm"),   
    pad_y = unit(2.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(
    location = "br"
  )


  # save under figures 
ggsave(filename = here("figures", "SDM_full_environment.png"), plot = plotSDM_full_env)


  #5b. environmental map with occurence points 
  # Occurrence records were overlaid on current habitat suitability predictions to visually assess model performance
  # occurrence points dataframe
gbif_points_df <- gbif_baMau %>%
  dplyr::select(decimalLongitude, decimalLatitude) %>%
  dplyr::rename(longitude = decimalLongitude,
                latitude = decimalLatitude)

  # final map with points 


install.packages("ggnewscale")  #in order to have probability and points in legend 
library(ggnewscale)

plotSDM_points <-
ggplot() +
  geom_raster(data = predict_baMau_df,
              aes(x = x, y = y, fill = lyr1)) +
  
  scale_fill_viridis_c(name = "Probability") +
  
  # reset fill scale
  ggnewscale::new_scale_fill() +
  geom_point(data = gbif_points_df,
             aes(x = longitude, y = latitude, fill = "Occurrence"),
             shape = 21,
             size = 1.8,
             stroke = 0.3,
             colour = "black",
             alpha = 0.85) +
  
  scale_fill_manual(
    name = "",
    values = c("Occurrence" = "orange"),
    guide = guide_legend(override.aes = list(shape = 21, size = 3))
  ) +
  
  annotation_spatial(Africa_map, fill = NA, colour = "black", linewidth = 0.3) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  ggtitle(expression(bold("Current Habitat Suitability of ") ~
                       bolditalic(B.~maughamii) ~
                       bold(" (Full Model) with Occurence Points"))) +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  
  coord_sf(expand = FALSE) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0, "cm"),
    pad_y = unit(2.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(location = "br")

  # save under figures 
ggsave(filename = here("figures", "SDM_Occurence_points.png"), plot = plotSDM_points)

  #5c. final map: Current Habitat Suitability of B. maughamii (Climate-Only Model)
plotSDM_climate  <- 
  ggplot() +
  geom_raster(data = predict_baMau_climate_df, 
              aes(x = x, y = y, 
                  fill = lyr1)) +
  annotation_spatial(Africa_map, fill = NA, colour = "black")+
  scale_fill_viridis_c(name = "Probability") +
  labs(x = "Longitude", 
       y = "Latitude") +
    ggtitle(expression(bold("Current Habitat Suitability of ") ~ bolditalic(B.~maughamii) ~ bold("(Climate-Only Model)"))) +
  theme_bw()+
  theme(
    plot.title = element_text(
      hjust = 0.5,        
      size = 11))+
  coord_sf(expand = FALSE) +  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0, "cm"),   
    pad_y = unit(2.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(
    location = "br"
  )

  # save in figures 
ggsave(filename = here("figures", "SDM_climate_only.png"), plot = plotSDM_climate)  


  #5d. final map: Projected Future Habitat Suitability of B. maughamii (SSP245, 2061–2080)
plotSDM_future <- 
ggplot() +
  geom_raster(data = predict_baMau_future_df, 
              aes(x = x, y = y, 
                  fill = lyr1)) +
  annotation_spatial(Africa_map, fill = NA, colour = "black")+
  scale_fill_viridis_c(name = "Probability") +
  labs(x = "Longitude", 
       y = "Latitude") +
  ggtitle(expression(bold("Projected Future Habitat Suitability of ") ~ bolditalic(B.~maughamii) ~ bold("(2061-2080)"))) +
  theme_bw()+
  theme(
    plot.title = element_text(
      hjust = 0.5,        
      size = 11))+
  coord_sf(expand = FALSE) +  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0, "cm"),   
    pad_y = unit(2.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(
    location = "br"
  )

  # save in figures 
ggsave(filename = here("figures", "SDM_future.png"), plot = plotSDM_future) 


