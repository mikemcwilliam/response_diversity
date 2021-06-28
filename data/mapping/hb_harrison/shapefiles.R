# Loading in the shapefile
#Queensland and GBR reefs
GBR_feat <- readOGR(dsn = "/Users/jc171846/Dropbox/Project Data/DECRA/Map/3dgbr_geomorph/shape", 
                    layer = "gbr_features",
                    verbose = F)
GBR_feat <- spTransform(GBR_feat, CRS("+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#Outline of the Coral Sea Marine Sanctuary
CS_feat <- readOGR(dsn = "/Users/jc171846/Dropbox/Project Data/DECRA/Map/3dgbr_geomorph/shape", 
                   layer = "qld_gbrwha_cscz",
                   verbose = F)
CS_feat <- spTransform(CS_feat, CRS("+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#Add dryreefs
CS_dryreef <- readOGR(dsn = "/Users/jc171846/Dropbox/Project Data/DECRA/Map/3dgbr_geomorph/shape", 
                      layer = "coralsea_dryreef",
                      verbose = F)
CS_dryreef <- spTransform(CS_dryreef, CRS("+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#Add reefs
CS_reef <- readOGR(dsn = "/Users/jc171846/Dropbox/Project Data/DECRA/Map/3dgbr_geomorph/shape", 
                   layer = "coralsea_reef",
                   verbose = F)
CS_reef <- spTransform(CS_reef, CRS("+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))