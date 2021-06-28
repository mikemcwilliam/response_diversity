


library("ggplot2")
library("plyr")
library("mapdata")
library("maps")
library("PBSmapping")
library("cowplot")

setwd("~/Documents/PhD/OneDrive - James Cook University/04_GBR/mapping")
load("mappingFunctions.RData")

# source of this code
# http://www.rpubs.com/spoonerf/countrymapggplot2


library(maptools)
library(raster)
library(rgdal)
library(ggplot2)
library(plyr)


# 1 choose country using the code: 
getData('ISO3')
aus<-getData("GADM", country="AUS", level=1)


# 2 set projection 
# choose UTM square from: http://www.dmap.co.uk/utmworld.htm 
# UTM code is longutudinal number followed by North (N) or South (S)
# for Lizard - UTM 55S
# search this UTM code here https://spatialreference.org 
# Extract EPSG number for projection
aus_UTM<-spTransform(aus, CRS("+init=EPSG:32355"))  #32555 #32755

# 3 localise the area
aus_UTM@data$NAME_1
qld<-aus_UTM[aus_UTM@data$NAME_1 == "Queensland",]
qld = spTransform(qld, "+init=epsg:4326")
qld_df<-fortify(qld)
          
          
# 4 select exact location and set limits
lizard_lat<--14.6645        
lizard_long<-145.4651
xlim<-c(lizard_long-0.04, lizard_long+0.02) 
ylim<-c(lizard_lat-0.04, lizard_lat+0.02)    
  
  
  
  
# 5 plot               
lizzie<- ggplot() + 
 geom_polygon(data=qld_df, aes(long,lat,group=group), fill="black")+
 geom_point(aes(x=145.451, y=-14.644),fill="red",shape=25,size=2, stroke=0.3)+
 theme(aspect.ratio=1)+xlim(min(qld_df$long),max(qld_df$long))+
 ylim(min(qld_df$lat),max(qld_df$lat))+
coord_cartesian(xlim,ylim)#+theme_opts
lizzie
