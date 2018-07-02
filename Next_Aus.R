# Loop Eucalyptus

setwd("E:/Uni_Wue_Master/Second_Semester/Spatial Modelling and Prediction")

#packages
library(Rcpp)
library(rlang)
library(dplyr)
library(purrr)
library(ggplot2)


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


# Get the shape file of Australia
data(wrld_simpl)# wrld_simpl is from maptools
wrld_simpl$NAME
Aus =wrld_simpl[wrld_simpl$NAME== "Australia",]
study_area <- Aus

writeOGR(Aus,"./Shapefile_Aus","Shapefile_Aus",driver = "ESRI Shapefile") 

study_area <- readOGR("./Shapefile_Aus/Shapefile_Aus.shp")


bio <- raster::getData("worldclim", var = "bio", res = 10)
biocrop <- crop(bio, extent(study_area) + 10)
writeRaster(biocrop,"bioclim.tif",format="GTiff")


biocrop <- brick("bioclim.tif")
biocrop

#' Collinearity
cm <- cor(getValues(biocrop), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))

#' ### Select an uncorrelated subset of environmental variables
env <- subset(biocrop, c("bioclim.7", "bioclim.13", "bioclim.14"))

Euc_species <- c("microcorys","viminalis","ovata", "propinqua","punctata","tereticornis","camaldulensis")
Euc_species <- c("microcorys")


for (i in 1:length(Euc_species)){
  
  species0 <- gbif("Eucalyptus",Euc_species[i])
  species <- subset(species0,select=c("lat","lon","year"))
  species  <- species[species$year >= 1970,] #exclude data before 1970
  species <- na.omit(species)
  coordinates(species) <- c("lon", "lat")  # set spatial coordinates
  
  
  
  # Add projection information
  proj4string(species) <- CRS("+proj=longlat +datum=WGS84")
 
  
  study_area <- spTransform(study_area,CRS=crs(species))
  
  crs(study_area)
  print(paste("CRS study_area is:",crs(study_area)))
  crs(species)
  print(paste("CRS species is:",crs(species)))
  
  
  Aus_E <- species[complete.cases(extract(biocrop, species)), ]
  print(paste("CRS Aus_E is:",crs(Aus_E)))
  
  writeOGR(Aus_E, "Eucalyptus_occurence_data", 
           layer =paste("Eucalyptus",sep="_",Euc_species[i],"gt1970"), driver="ESRI Shapefile",
           overwrite_layer = T) 
  
  
  
  Eucalyptus <- Aus_E
  
  

  #' Selecting 2000 random background points, excluding cells where
  #' the species is present
  set.seed(2)
  background <- randomPoints(env, 2000, Eucalyptus) 
  #' Select only one presence record in each cell of the environmental layer
  presence <- gridSample(Eucalyptus, env, n = 1)
 

  #' 
  #' Now we combine the presence and background points, adding a 
  #' column "species" that contains the information about presence (1)
  #' and background (0)
  fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                     data = data.frame("Eucalyptus" = rep(c(1,0), 
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
  

  
  
  ########### GAM #####################################################################################
  
  #' We can now use a range of statistical methods to estimate the
  #' probability of species occurrence.
  #' Unfortunately, there are often subtle differences in how the models
  #' are specified and in which data formats are useable
  
  
  
  varnames <- c("bioclim.7", "bioclim.13", "bioclim.14")
  
  ## Generalized Linear Model
  
  ## Generalized additive models
  gammodel <- gam(Eucalyptus ~ s(bioclim.7) + s(bioclim.13) + s(bioclim.14),
                  family="binomial", data=traindata)
  summary(gammodel)
  gammodel
  
  
  plot(gammodel, pages=1,ylim=c(-100,100))
 
  # Evaluate model on test data
  # a) Predict to test data
  gamtest <- predict(gammodel, newdata = testdata, type = "response")
  
  
  # b) Calculate performance indices
  #png(filename= paste("./Eucalyptus_results/Performance_Indices_Eucalyptus_",sep="",Euc_species[i],".png"))
  val.prob(gamtest, testdata[["Eucalyptus"]])
  #dev.off() 
  
  
  # Variable importance
  gamimp <- varImpBiomod(gammodel, varnames,
                         traindata)
  
  #png(filename = paste("./Eucalyptus_results/Variable_Importance_Eucalyptus_",sep="",Euc_species[i],".png"))
  barplot(100 * gamimp/sum(gamimp), ylab = "Variable importance (%)",
          main = paste("Variable Importance Eucalyptus", Euc_species[i]))
  #dev.off()
  
  
  # Prediction map
  gammap <- predict(env, gammodel, type = "response")
  
  
 
  
  

  # writeRaster(gammap, paste("./Eucalyptus_results/Pred_Map_Eucalyptus",
  #                   sep="_",Euc_species[i],"gt1970"),format="GTiff",overwrite=T)
  
  
  
  p <- print(ggR(gammap,geom_raster = T,stretch="hist",quantiles=c(0.05,0.95))+
               scale_fill_gradientn(colours = rev(terrain.colors(5)),limits=c(0,1),  
                                    space = "Lab",name=paste("Probability \n"),na.value = NA)+
               scale_x_continuous("Lon", sec.axis = dup_axis())+
               scale_y_continuous("Lat", sec.axis = dup_axis())+
               theme(plot.title=element_text(color="red",hjust =0.5))+
               labs(title= paste("Prediction map Eucalyptus",Euc_species[[i]])))
  
  ggsave(p$plot, filename= paste("./Eucalyptus_results/gg_Pred_Map_Eucalyptus_",sep ="",Euc_species[[i]],".png"))
  
  
  print(paste("Processing done for Eucalyptus", Euc_species[i])) 
  
}





