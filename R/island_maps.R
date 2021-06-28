

#library(maptools)
#library(raster)
#library(rgdal)
#library("mapdata")
#library("maps")
#library("PBSmapping")
#load("figs/mapping/mappingFunctions.RData")
#source("figs/mapping/regroup_close.R")

# PACIFIC OCEAN WITH NO LINES
ylim=c(-35, 30)
xlim=c(70, 350)
center=160

worldmap <- map_data ("world")
worldmap$long <- ifelse(worldmap$long < center - 180 , worldmap$long + 360, worldmap$long)
worldmap <- ddply(worldmap, .(group), RegroupElements, "long", "group") # now regroup
worldmap <- ddply(worldmap, .(group), ClosePolygons, "long", "order") # close polys 
colnames(worldmap)<-c("X","Y","PID","POS","region","subregion", "PID2")
worldmap <- clipPolys(worldmap, xlim=xlim,ylim=ylim, keepExtra=TRUE)
head(worldmap)

sites<-data.frame(site=c("rio","tia","liz"), 
lat=c(18.4745,-17.4928, -14.6645),
long=c(-77.4673, -149.8977, 145.4651))
sites$long <- ifelse(sites$long < center - 180 , sites$long + 360, sites$long)
sites

map_theme<-theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.line=element_blank())

pacific<-ggplot()+ 
coord_map(ylim=ylim, xlim=xlim) +
geom_polygon(data=worldmap,aes(X,Y,group=PID2),fill="white",color="black",size=0.035) +
geom_point(data=sites, aes(x=long, y=lat, col=site), shape=22, size=2, stroke=1.2)+
geom_point(data=sites, aes(x=long, y=lat, col=site), shape=15, size=2, alpha=0.2)+
guides(fill=F,size=F, col=F)+map_theme+ theme(panel.background=element_rect(fill="whitesmoke"))+
scale_colour_manual(values=c(gbrcol, caribcol, polycol))
pacific

# ISLANDS 
aus<-getData("GADM", country="AUS", path="figs/mapping", level=1)
pol<-getData("GADM", country="PYF", path="figs/mapping", level=1)
jam<-getData("GADM", country="JAM", path="figs/mapping", level=0)


aus_UTM<-spTransform(aus, CRS("+init=EPSG:32355")) 
qld<-aus_UTM[aus_UTM@data$NAME_1 == "Queensland",]
qld = spTransform(qld, "+init=epsg:4326")
qld_df<-fortify(qld)
liz_lat<--14.6645        
liz_long<-145.4651
L_xlim<-c(liz_long-0.04, liz_long+0.02) 
L_ylim<-c(liz_lat-0.04, liz_lat+0.024)              
          
pol_UTM<-spTransform(pol, CRS("+init=EPSG:2977"))  
mor<-pol_UTM[pol_UTM@data$NAME_1 == "Îles du Vent",]
mor = spTransform(mor, "+init=epsg:4326")
mor_df<-fortify(mor)  
mor_lat<--17.5388    
mor_long<--149.8295
M_xlim<-c(mor_long-0.1, mor_long+0.1) 
M_ylim<-c(mor_lat-0.1, mor_lat+0.1)   
 
jam_UTM<-spTransform(jam, CRS("+init=EPSG:3450")) 
#jam<-jam_UTM[jam_UTM@data$NAME_1 == "Saint Ann",]
jam = spTransform(jam, "+init=epsg:4326")
jam_df<-fortify(jam)
jam_lat<-18.1096 
jam_long<--77.2975
J_xlim<-c(jam_long-1, jam_long+1) 
J_ylim<-c(jam_lat-1, jam_lat+1)  

# 1 degree = 111 km
# 0.01 degree = 1.11 km
# 0.0045 degree = 0.5 km LIZ
# 0.018 degree = 2km MOR
# 0.18 degree = 20km JAM

pliz<-ggplot() + geom_polygon(data=qld_df, aes(long,lat,group=group), fill="whitesmoke", col='black', size=0.1)+
 geom_point(aes(x=145.453, y=-14.641),fill="red",shape=25,size=1, stroke=0.1)+
 #geom_segment(aes(x=145.453+0.00225, xend=145.453-0.00225,
 #y=-14.643+0.004,yend=-14.643+0.004), size=0.5)+
 theme(aspect.ratio=1)+xlim(min(qld_df$long),max(qld_df$long))+
 ylim(min(qld_df$lat),max(qld_df$lat)+0.05)+
coord_cartesian(L_xlim,L_ylim)+map_theme


pmor<-ggplot() + geom_polygon(data=mor_df, aes(long,lat,group=group), fill="whitesmoke", col='black', size=0.1)+
 geom_point(aes(x=-149.8977, y=-17.473),fill="red",shape=25,size=1, stroke=0.1)+
 #geom_segment(aes(x=-149.8977-0.009, xend=-149.8977+0.009,
 #y=-17.475+0.015,yend=-17.475+0.015), size=0.5)+
 theme(aspect.ratio=1)+xlim(min(mor_df$long),max(mor_df$long))+
 ylim(min(mor_df$lat),max(mor_df$lat))+
coord_cartesian(M_xlim,M_ylim)+map_theme

pjam<-ggplot() + geom_polygon(data=jam_df, aes(long,lat,group=group), fill="whitesmoke", col='black', size=0.1)+
 geom_point(aes(x=-77.44, y=18.64),fill="red",shape=25,size=1, stroke=0.1)+
 #geom_segment(aes(x=-77.44-0.09, xend=-77.44+0.09,
# y=18.62+0.15,yend=18.62+0.15), size=0.5)+
 theme(aspect.ratio=1)+xlim(min(jam_df$long),max(jam_df$long))+
 ylim(min(jam_df$lat),max(jam_df$lat)+1)+
coord_cartesian(J_xlim,J_ylim)+map_theme



plot_grid(pacific, plot_grid(pliz, pmor, pjam, nrow=1), nrow=2)

