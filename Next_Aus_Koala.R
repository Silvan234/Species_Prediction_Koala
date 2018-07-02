# Loop Koala

setwd("E:/Uni_Wue_Master/Second_Semester/Spatial Modelling and Prediction")

#packages
library(Rcpp)
library(rlang)
library(dplyr)
library(purrr)


library(rms)
library(raster)
library(mgcv) # gam
source("varImpBiomod.R")
library(randomForest)
library(dismo)
library(rgdal)
library(ellipse)
library(randomForest)
library(rJava)
library(XML)
library(RStoolbox)
library(stringr)
library(sp)
library(rworldmap)

study_area <- readOGR("./Shapefile_Aus/Shapefile_Aus.shp")

biocrop <- brick("bioclim.tif")

#' Collinearity
cm <- cor(getValues(biocrop), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))

if (file.exists("./Species_occurence_data/Koala_species_gt1970.shp")) {
  species <- readOGR("./Species_occurence_data/Koala_species_gt1970.shp")
} else {
  
  species0 <- gbif('Phascolarctos', 'cinereus')
  species <- subset(species0,select=c("lat","lon","year"))
  species  <- species[species$year >= 1970,] #exclude data before 1970
  species <- na.omit(species)
  coordinates(species) <- c("lon", "lat")  # set spatial coordinates
  
  # Add projection information
  proj4string(species) <- CRS("+proj=longlat +datum=WGS84")
  
  # Save species records in mif-format (preserves full column names)
  writeOGR(species, "Species_occurence_data", 
           "Koala_species_gt1970", driver="ESRI Shapefile")
}

env <- subset(biocrop, c("bioclim.12", "bioclim.17"))

set.seed(2)
background <- randomPoints(env, 2000, species) 
#' Select only one presence record in each cell of the environmental layer
presence <- gridSample(species, env, n = 1)

#' 
#' Now we combine the presence and background points, adding a 
#' column "species" that contains the information about presence (1)
#' and background (0)
fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("species" = rep(c(1,0), 
                                                                     c(nrow(presence), nrow(background)))),
                                   match.ID = FALSE,
                                   proj4string = CRS(projection(env)))

fulldata@data
#' Add information of environmental conditions at point locations
fulldata@data <- cbind(fulldata@data, extract(env, fulldata))
fulldata@data

#' 
# Split data set into a training and test data set
set.seed(2)
fold <- kfold(fulldata, k = 5)
traindata <- fulldata[fold != 1, ]
testdata <- fulldata[fold == 1, ]

traindata
testdata

########### GAM #####################################################################################
#dir.create("Koala_results")

varnames <- c("bioclim.12",  "bioclim.17")

## Generalized additive models
gammodel <- gam(species ~ s(bioclim.12)  + s(bioclim.17),
                family="binomial", data=traindata)
summary(gammodel)

plot(gammodel, pages=1,ylim=c(-100,100))

# Evaluate model on test data
# a) Predict to test data
gamtest <- predict(gammodel, newdata = testdata, type = "response")


# b) Calculate performance indices

png(filename= paste("./Koala_results/Performance_Indices_Koala_Clim.png"))
val.prob(gamtest, testdata[["species"]])
dev.off() 

# Variable importance
gamimp <- varImpBiomod(gammodel, varnames,
                       traindata)
gamimp

png(filename = paste("./Koala_results/Variable_Importance_Koala_Clim.png"))
barplot(100 * gamimp/sum(gamimp), ylab = "Variable importance (%)",
        main = "Variable Importance Koala")
dev.off()





# Prediction map
gammap <- predict(env, gammodel, type = "response")


plot(gammap)

#writeRaster(gammap, paste("./Koala_results/Pred_Map_Koala_Clim"),format="GTiff",overwrite=T)

Gammap_Clim <- print(ggR(gammap,geom_raster = T,stretch="hist",quantiles=c(0.05,0.95))+
             scale_fill_gradientn(colours = rev(terrain.colors(5)),limits=c(0,1),  
                                  space = "Lab",name=paste("Probability \n"),na.value = NA)+
             scale_x_continuous("Lon", sec.axis = dup_axis())+
             scale_y_continuous("Lat", sec.axis = dup_axis())+
             theme(plot.title=element_text(color="red",hjust =0.5))+
             labs(title= paste("Prediction map Koala Climate ")))


ggsave(Gammap_Clim$plot, filename= paste("./Koala_results/gg_Pred_Map_Koala_Climate.png"))




##########################################################

# Loop Koala Climate and Eucalyptus

##########################################################

All_Eucal <- list.files("./Eucalyptus_results", full.names = T, pattern=".tif")
All_Eucal
All_Eucal_stack <- stack(All_Eucal)

names(All_Eucal_stack) <- c("E.camaldulensis","E.microcorys","E.ovata", "E.propinqua", "E.punctata",
                            "E.tereticornis", "E.viminalis")

names(All_Eucal_stack)

Euc_Biocrop <- stack(biocrop,All_Eucal_stack)
Euc_Biocrop
names(Euc_Biocrop)


env <- subset(Euc_Biocrop, c("bioclim.12", "bioclim.17", "E.camaldulensis",
                         "E.microcorys","E.ovata", "E.propinqua", "E.punctata",
                         "E.tereticornis", "E.viminalis"))

if (file.exists("./Species_occurence_data/Koala_species_gt1970.shp")) {
  species <- readOGR("./Species_occurence_data/Koala_species_gt1970.shp")
} else {
  
  species0 <- gbif('Phascolarctos', 'cinereus')
  species <- subset(species0,select=c("lat","lon","year"))
  species  <- species[species$year >= 1970,] #exclude data before 1970
  species <- na.omit(species)
  coordinates(species) <- c("lon", "lat")  # set spatial coordinates
  
  # Add projection information
  proj4string(species) <- CRS("+proj=longlat +datum=WGS84")}


set.seed(2)
background <- randomPoints(env, 2000, species) 
#' Select only one presence record in each cell of the environmental layer
presence <- gridSample(species, env, n = 1)

#' 
#' Now we combine the presence and background points, adding a 
#' column "species" that contains the information about presence (1)
#' and background (0)
fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("species" = rep(c(1,0), 
                                                                     c(nrow(presence), nrow(background)))),
                                   match.ID = FALSE,
                                   proj4string = CRS(projection(env)))

fulldata@data
#' Add information of environmental conditions at point locations
fulldata@data <- cbind(fulldata@data, extract(env, fulldata))
fulldata@data

#' 
# Split data set into a training and test data set
set.seed(2)
fold <- kfold(fulldata, k = 5)
traindata <- fulldata[fold != 1, ]
testdata <- fulldata[fold == 1, ]

traindata
testdata

########### GAM #####################################################################################

varnames <- c("bioclim.12", "bioclim.17", "E.camaldulensis",
              "E.microcorys","E.ovata", "E.propinqua", "E.punctata",
              "E.tereticornis", "E.viminalis")

## Generalized additive models
gammodel <- gam(species ~ s(bioclim.12)  + s(bioclim.17) + s(E.camaldulensis) + s(E.microcorys)
                + s(E.ovata)+ s(E.propinqua)+ s(E.punctata)+ s(E.tereticornis)+ s(E.viminalis),
                family="binomial", data=traindata)
                
summary(gammodel)

plot(gammodel, pages=1,ylim=c(-100,100))

# Evaluate model on test data
# a) Predict to test data
gamtest <- predict(gammodel, newdata = testdata, type = "response")


# b) Calculate performance indices

png(filename= paste("./Koala_results/Performance_Indices_Koala_Clim_Eucal.png"))
val.prob(gamtest, testdata[["species"]])
dev.off() 

# Variable importance
gamimp <- varImpBiomod(gammodel, varnames,
                       traindata)
gamimp

png(filename = paste("./Koala_results/Variable_Importance_Koala_Clim.png"))
barplot(100 * gamimp/sum(gamimp), ylab = "Variable importance (%)",
        main = "Variable Importance Koala")
dev.off()





# Prediction map
gammap <- predict(env, gammodel, type = "response")


plot(gammap)

#writeRaster(gammap, paste("./Koala_results/Pred_Map_Koala_Clim_Eucal"),format="GTiff",overwrite=T)

Gammap_ClimEuc <- print(ggR(gammap,geom_raster = T,stretch="hist",quantiles=c(0.05,0.95))+
             scale_fill_gradientn(colours = rev(terrain.colors(5)),limits=c(0,1),  
                                  space = "Lab",name=paste("Probability \n"),na.value = NA)+
             scale_x_continuous("Lon", sec.axis = dup_axis())+
             scale_y_continuous("Lat", sec.axis = dup_axis())+
             theme(plot.title=element_text(color="red",hjust =0.5))+
             labs(title= paste("Prediction map Koala Climate and Eucalyptus ")))

ggsave(Gammap_ClimEuc$plot, filename= paste("./Koala_results/gg_Pred_Map_Koala_Climate_Eucal.png"))








# corLocal

Koala_PM_Clim    <- raster("./Koala_results/Pred_Map_Koala_Clim.tif")
Koala_PM_ClimEuc <- raster("./Koala_results/Pred_Map_Koala_Clim_Eucal.tif") 

Koala_PM_corLocal <- corLocal(Koala_PM_Clim,Koala_PM_ClimEuc,method="spearman")
warnings()

plot(Koala_PM_corLocal)

p <- print(ggR(gammap,geom_raster = T,stretch="hist",quantiles=c(0.05,0.95))+
             scale_fill_gradientn(colours = rev(terrain.colors(5)),limits=c(-1,1),  
                                  space = "Lab",name=paste("Correlation \n coefficient \n"),na.value = NA)+
             scale_x_continuous("Lon", sec.axis = dup_axis())+
             scale_y_continuous("Lat", sec.axis = dup_axis())+
             theme(plot.title=element_text(color="red",hjust =0.5))+
             labs(title= paste("Map CorLocal Koala")))

ggsave(p$plot, filename= paste("./Koala_results/gg_CorLocal_Koala.png"))
