

#sitedata[,c("1995 Lizard 1-7m","2002 Lizard 1-7m","2011 Lizard 1-7m")]
#----sitedata[,c("1995 Lizard 1-7m")]<-sitedata[,c("1996 Lizard 1-7m")]
#----sitedata[c("Galaxea","Astreopora","Mussidae"),c("1995 Lizard 1-7m","2002 Lizard 1-7m","2011 Lizard 1-7m")]<-0
sitedata[is.na(sitedata)]<-0
#sitedata[,"1995 Lizard 1-7m"]<-rowMeans(cbind(data.frame(sitedata[,"1995 Lizard 1-7m"]),data.frame(sitedata[,"1996 Lizard 1-7m"])))
#sitedata<-sitedata[apply(sitedata, 1, function(x) !all(abs(x)<1)),] # RARE



# coral abundances at 3 times in 3 places

library("reshape2")

substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x))}

cols<-c("Region", "Site", "Year", "Zone", "Replicate", "Taxon","x")

taxa <- read.csv("data/align_taxa.csv")
zones <- read.csv("data/align_zones.csv")

# rio bueno, jamaica

jam<-read.csv("data/original/jamaica.csv")
jam$Taxon <- taxa$jam_new[match(jam$Taxon, taxa$jam_old)]
jam <- jam[!is.na(jam$Taxon),]
jam <- aggregate(.~Taxon, jam, sum)

jam <- melt(jam, id.vars="Taxon", value.name="x", variable.name="Site")
jam$Year <- substrRight(as.character(jam$Site), 4)
jam$Zone <- zones$jam_new[match(jam$Site, zones$jam_old)]
jam$Replicate <- 1
jam$Site <- "Rio Bueno"
jam$Region <- "Jamaica"
jam <- jam[,cols]
head(jam)

# moorea, french polynesia

fp<-read.csv("data/original/Benthic_data_Moorea.csv")
head(fp)
colnames(fp)
fp$Year <- substrRight(as.character(fp$Date), 4)
fp$Replicate <- fp$Trans
fp$Site <- fp$Location
fp<-melt(fp, id.vars=c("Year","Site", "Zone","Replicate"), value.name="x", variable.name="Taxon")

fp$Taxon <- taxa$fp_new[match(fp$Taxon, taxa$fp_old)]
fp$Zone <- zones$fp_new[match(fp$Zone, zones$fp_old)]
fp <- fp[!is.na(fp$Taxon),]
fp$x <- as.numeric(fp$x)
fp$Region <- "Polynesia"
fp <- fp[,cols]
fp <- aggregate(x~., fp, sum)
fp <- fp[,cols]
head(fp)

# Literature FP data (Bouchon 1981)
# Bouchon data is presented as % relative abundance
# multiply by total coral cover to get absolute abundance  
# cover values not presented in bouchon
# cover presented without millepora in berumen
# berumencover/bouchoncover = no_millepora_proportion
# bouchoncover = bermumen cover/no_millepora_proportion

fp_b<-read.csv("data/literature/bouchon_Moorea.csv")
berumen <- read.csv("data/literature/berumen.csv")
cover <- data.frame(ID = unique(fp_b$ID))
cover$berumen<-berumen$Total[match(cover$ID, berumen$ID)]
cover$millepora<-aggregate(Cover~ID, fp_b[fp_b$Genus=="Millepora",], sum)$Cover
cover$bouchon <- cover$berumen/(1-(cover$millepora/100))
cover

fp_b$Total <- cover$bouchon[match(fp_b$ID, cover$ID)]
fp_b$Cover <- (fp_b$Cover/100)*fp_b$Total

fp_b$Replicate <- 1
fp_b$Site <- "Tiahura"
fp_b$x <- fp_b$Cover
fp_b$Taxon <- taxa$fpBouch_new[match(fp_b$Taxa, taxa$fpBouch_old)]
fp_b$Zone <- zones$fpBouch_new[match(fp_b$Zone, zones$fpBouch_old)]
fp_b$Region <- "Polynesia"
fp_b <- fp_b[!is.na(fp_b$Taxon),]
fp_b <- fp_b[,cols]
fp_b <- aggregate(x~., fp_b, sum)

# lizard island, gbr

gbr<-read.csv("data/original/Coral Data Lizard Island 95-17.csv",stringsAsFactors=FALSE)
gbr<-melt(gbr, id.vars=c("Year","Site", "Zoneno","Replicate"), value.name="x", variable.name="Taxon")

gbr$Taxon <- taxa$gbr_new[match(gbr$Taxon, taxa$gbr_old)]
gbr <- gbr[!is.na(gbr$Taxon),]
gbr$Zone <- zones$gbr_new[match(gbr$Zoneno, zones$gbr_old)]
gbr$Region <- "GBR"
gbr <- gbr[,cols]
gbr$x <- as.numeric(gbr$x)
gbr <- aggregate(x~., gbr, sum)



# Additional GBR data (Lizard 1996)

gbr_b <- read.csv("data/original/Coral surveys Lizard 1995-99.csv",stringsAsFactors=FALSE)
gbr_b<-melt(gbr_b, id.vars=c("Month","Site", "Zone","Replicate"), value.name="x", variable.name="Taxon")
gbr_b$Site<-ifelse(gbr_b$Site=="NR", "1. North Reef", ifelse(gbr_b$Site=="WM", "2. Washing Machine",ifelse(gbr_b$Site=="LH", "3. Lizard Head", ifelse(gbr_b$Site=="SI", "4. South Island", NA))))
gbr_b$Year<-paste("19",substrRight(gbr_b$Month, 2),
"b", sep="")

gbr_b$Taxon <- taxa$gbr_new[match(gbr_b$Taxon, taxa$gbr_old)]
gbr_b <- gbr_b[!is.na(gbr_b$Taxon),]
gbr_b$Zone <- zones$gbr96_new[match(gbr_b$Zone, zones$gbr96_old)]
gbr_b$Region <- "GBR"
gbr_b <- gbr_b[,cols]

# merge 
abun <- rbind(jam, fp, fp_b, gbr)
abun$ID <- apply(abun[ ,c("Year","Site","Zone")], 1, paste, collapse=".")
abun$IDsite<-apply(abun[ ,c("Site","Zone")], 1, paste, collapse=".")
head(abun)

# subset
sites <- c("1. North Reef", "Rio Bueno", "Tiahura")
abun <- abun[abun$Site %in% sites,]
depths <- c("1-7m", "7-15m", "15-30m")
abun <- abun[abun$Zone %in% depths,]

#moorea<-subset(moorea, Taxon!="Echinopora" & Taxon!="Stylocoeniella") 

# aggregate replicates.... 
avabun <- aggregate(x~., subset(abun, select=-Replicate), mean)

# Define three points
avabun$Points<-
ifelse(avabun$Region=="Jamaica" & avabun$Year==1977, 1, 
ifelse(avabun$Region=="Jamaica" & avabun$Year==1993, 2, 
ifelse(avabun$Region=="Jamaica" & avabun$Year==2013, 3, 
ifelse(avabun$Region=="GBR" & avabun$Year==1995, 1, 
#ifelse(avabun$Region=="GBR" & avabun$Year==1996, 1.5, 
ifelse(avabun$Region=="GBR" & avabun$Year==2002, 2,
ifelse(avabun$Region=="GBR" & avabun$Year==2011, 3,
ifelse(avabun$Region=="Polynesia" & avabun$Year==1979, 1, 
ifelse(avabun$Region=="Polynesia" & avabun$Year==1982, 2,
ifelse(avabun$Region=="Polynesia"& avabun$Year==2003, 3, 0)))))))))

head(avabun)

write.csv(avabun, file="data/abundance.csv")


##### TRAIT DATA    

traits<-read.csv("data/traits/traits.csv")
rownames(traits)<-traits$Group.1

traits$X<-NULL
traits$Group.1<-NULL
traitlist<-traits[,c("height","sa_vol","spaces", "polypcat","growthcat","skelcat","sizecat")]
traitdata<-data.frame(sapply(traitlist, function(x) as.numeric(as.character(x))))
sapply(traitdata, class)

# Export data
head(traits)
colnames(traits)
export<-traits[,c(1,2,3,4,9,10,11,12,13,14,15,16)]
colnames(export)
colnames(export)<-c("CS_raw","CW_raw","GR_raw","SD_raw","Morphology", "CH_cat","SV_cat","IB_cat","CW_cat","GR_cat","CS_cat","SD_cat")

write.csv(export, file="data/traits.csv")
