
####### dip your feet in
rm(list=ls())
library("ggplot2")
library("cowplot")
library("reshape2")
library("FD") 
library("ggalt")
library("png")
library("grid")
library("FNN")


# data
#source("R/data_prep.R")
abun <- read.csv("data/abundance.csv")
traits <- read.csv("data/traits.csv")

# colours
cols <- c("#fbb4ae", "#b3cde3", "#bee2b7")
names(cols) <- c("Jamaica","GBR", "Polynesia")

################################
# coral cover relative to t0

sites <- aggregate(x~., subset(abun, select=-c(X, Taxon)), sum)

t0 <- NULL
for(x in unique(sites$IDsite)){
sub <- sites[sites$IDsite==x,]
cov <- sub[which.min(sub$Year),"x"]
t0 <- rbind(t0, data.frame(IDsite=x, t0=cov))}
sites$t0 <- t0$t0[match(sites$IDsite, t0$IDsite)]
sites$rel_cover <- (sites$x-sites$t0)/sites$t0*100

ggplot(sites, aes(as.numeric(Year), rel_cover))+
geom_line(aes(col=Zone))+
stat_summary(geom="line", fun="mean")+
facet_wrap(~Region)

################################
# trait diversity

tlist <- c("CW_cat", "GR_cat",  "SD_cat", "CH_cat", "SV_cat", "IB_cat", "CS_cat")
gower<-gowdis(traits[,tlist], w=c(1,1,1,1,1,1,0.5))
pca<-pcoa(gower)
traits$pca1<--pca$vectors[,1]
traits$pca2<--pca$vectors[,2]
traits$pca3<-pca$vectors[,3]
traits$pca4<-pca$vectors[,4]

a <- acast(abun, Taxon~ID, value.var="x")
a[is.na(a)]<-0
a <- a[!as.vector(rowSums(a)==0),]

t <- traits[,c("pca1","pca2", "pca3")]
rownames(t) <- traits$X
t <- t[rownames(a),]
rownames(t) == rownames(a)

fd<-as.data.frame(dbFD(a=t(a), x=dist(t,method = "euclidean")))
cwms<-functcomp(a=as.matrix(t(a)), x=t)
fd$ID<-rownames(fd)
sitefd<-merge(x=fd, y=sites, by="ID", all.x=TRUE, all.y=TRUE)
head(sitefd)

ggplot(sitefd[sitefd$Points>0,], aes(as.numeric(Year), FDis, col=Region))+stat_summary(fun="median", geom="line")+stat_summary(fun="median", geom="point")

################################
# relative FD change

fd0 <- NULL
for(x in unique(sitefd$IDsite)){
sub <- sitefd[sitefd$IDsite==x,]
fd <- sub[which.min(sub$Year),"FDis"]
fd0 <- rbind(fd0, data.frame(IDsite=x, fd0=fd))}
sitefd$fd0 <- fd0$fd0[match(sitefd$IDsite, fd0$IDsite)]
sitefd$rel_fd <- (sitefd$FDis-sitefd$fd0)/sitefd$fd0*100

recov <- subset(sitefd, Points==3)
recov <- melt(recov[,c("Region", "Zone","rel_fd","rel_cover")])
ggplot(recov, aes(variable, value, group=Zone))+
geom_bar(stat="Identity", position=position_dodge(), col="white")+facet_wrap(~Region, ncol=1, scales="free_y")

################################
# figure 2
#source("R/island_maps.R")
source("R/fig_2.R")
fig2
#ggsave("figs/fig2.png", fig2, width=130, height=145, units="mm", dpi=200)

################################
# trait space

cmd<-cmdscale(gower, add=TRUE, k=3)
fit<-envfit(cmd, traits[,tlist], perm=1000, na.rm=TRUE)
vects<-as.data.frame(scores(fit, "vectors"))
vects$lab<-c("CW","GR","SD","CH","SV","IB","CS")
vects

ggplot()+
stat_density_2d(data=traits, aes(x=pca1, y=pca2), size=0.15, alpha=0.75, bins=6)+
geom_point(data=traits, aes(x=pca1, y=pca2), col="red", size=0.75, alpha=0.5)+
geom_text(data=vects, aes(-Dim1/3.5, -Dim2/3.5, label=lab))+
ylim(c(-0.27,0.3))+xlim(-0.42, 0.5)+
labs(x="PCoA axis 1", y="PCoA axis 2")

################################
# trait shifts (slope)

slope <- subset(abun, Zone=="7-15m" & Points>0)
slope[,c("pca1","pca2")] <- traits[,c("pca1","pca2")][match(slope$Taxon, traits$X),]
slope[,c("cwm1","cwm2")]<- cwms[,c("pca1","pca2")][match(slope$ID, rownames(cwms)),]
head(slope)

ggplot()+
geom_encircle(data=subset(slope, x>0), aes(x=pca1, y=pca2,fill=Region), alpha=0.5)+
geom_segment(data=subset(slope, x>0), aes(x=cwm1, y=cwm2, xend=pca1, yend=pca2, col=Region))+
geom_point(data=subset(slope, x>0), aes(x=pca1, y=pca2, size=x, fill=Region), shape=21)+
facet_grid(Points~Region, switch="y")

################################
# response diversity 1

rd<-dcast(Region+Taxon+pca1+pca2~Points, data=slope, value.var=("x"))
rd<-data.frame(rd)
rd$Delta<-rd[,"X3"]-rd[,"X1"]
rd$WinLose<-ifelse(rd$Delta>0, "Gain","Loss")
rd$Delta<-abs(rd$Delta)

rd_a<-acast(Taxon~Region+WinLose, data=rd, value.var="Delta")
rd_a <- rd_a[rownames(t),]
rownames(rd_a)==rownames(t)
rd_cwm<-functcomp(a=as.matrix(t(rd_a)), x=t)
rd[,c("cwm1","cwm2")]<- rd_cwm[,c("pca1","pca2")][match(paste(rd$Region, rd$WinLose, sep="_"), rownames(rd_cwm)),]
head(rd)

ggplot()+
geom_encircle(data=rd, aes(x=pca1, y=pca2,fill=WinLose), alpha=0.5)+
geom_segment(data=rd,aes(x=cwm1, y=cwm2, xend=pca1, yend=pca2, col=WinLose))+
geom_point(data=rd, aes(x=pca1, y=pca2, size=Delta, fill=WinLose),shape=21)+ 
facet_wrap(~Region)

RDcalc <- function(R){
#R="Polynesia"
df<-subset(rd, Region==R)
df<-subset(df, Delta>0)
delta<-acast(WinLose~Taxon, data=df, value.var="Delta")
delta[is.na(delta)]<-0
fd<-data.frame(FDis=as.matrix(fdisp(a=delta, dist(data.frame(df[,c("pca1","pca2")], row.names=df$Taxon)))$FDis))
fd$WinLose<-rownames(fd)
cwm<-acast(WinLose~variable, data=unique(melt(df[,c("WinLose","cwm1","cwm2")], id.var="WinLose")), value.var="value")
a2<-(cwm["Gain","cwm1"]-cwm["Loss","cwm1"])^2
b2<-(cwm["Gain","cwm2"]-cwm["Loss","cwm2"])^2
cwmdist<-sqrt(a2+b2)
fd$cwmdist<-cwmdist
fd$Region<-R
fd}

rddat<-rbind(RDcalc("GBR"),RDcalc("Polynesia"),RDcalc("Jamaica"))
rddat$Region<-factor(rddat$Region, levels=c("GBR","Polynesia","Jamaica"))
rddat

################################
# figure 3
source("R/fig_3.R")
Fig3
#ggsave("figs/fig3.png", Fig3, width=155, height=130, units="mm", dpi=300)

################################
# response diversity 2

rd$Mort<-(rd[,"X1"]-rd[,"X2"])
rd$Recov<-rd[,"X3"]-rd[,"X2"]
rd$Replace<-rd[,"X3"]-rd[,"X1"]

tX1<-aggregate(X1~Region, rd, sum)
tX2<-aggregate(X2~Region, rd, sum)
tX3<-aggregate(X3~Region, rd, sum)
rd$tX1<-tX1$X1[match(rd$Region, tX1$Region)]
rd$tX2<-tX2$X2[match(rd$Region, tX2$Region)]
rd$tX3<-tX3$X3[match(rd$Region, tX3$Region)]

rd$relReplace<-(rd$X3/rd$tX3)-(rd$X1/rd$tX1)

maxReplace<-aggregate(rd$Replace, by=list(Region=rd$Region), FUN=max)
rd$relReplace2<-rd$Replace/(maxReplace$x[match(rd$Region, maxReplace$Region)])

################################
# rd3 

# DISTINCTIVENESS/CHANGE IN COVER

nnsum<-function(R, n){
#R="Polynesia"
#n=2
dat<-subset(rd, Region==R)
dat$ID<-c(1:nrow(dat))
dat<-dat[!is.na(dat$pca1),]
dist<-dist(dat[,c("pca1","pca2")])
nearest=get.knn(dist, k=n)$nn.index
rownames(nearest)<-dat$Taxon
nn<-data.frame(nearest)
names<-data.frame(t(apply(nn, 1, function(x) dat$Taxon[match(x, dat$ID)])))
colnames(names)<-letters[seq( from = 1, to = n )]
shift<-data.frame(t(apply(nn, 1, function(x) dat$Replace[match(x, dat$ID)])))
shift$sum<-rowMeans(shift, na.rm=T)
shift$r_sum=shift$sum/max(abs(shift$sum),na.rm=T)
dat$r_rep=dat$Replace/max(abs(dat$Replace), na.rm=T)
df<-cbind(dat, shift[,c("sum","r_sum")], n=n)
#fit<-lm(r_rep~r_sum, df)
#df$Rsq<-summary(fit)$r.squared
#plot(fit)
fit<-cor.test(df$r_sum,df$r_rep,method="spearman")
df$r<-fit$estimate
df$p<-fit$p.value
df}

nndat8<-rbind(nnsum("GBR", 8),nnsum("Polynesia", 8),nnsum("Jamaica", 8))
nndat4<-rbind(nnsum("GBR", 4),nnsum("Polynesia", 4),nnsum("Jamaica", 4))
nndat2<-rbind(nnsum("GBR", 2),nnsum("Polynesia", 2),nnsum("Jamaica", 2))
nndat<-rbind(nndat8, nndat4, nndat2)
head(nndat)

stat<-unique(nndat[,c("r","p","Region","n")])
stat$p<-ifelse(stat$p>0.01,as.character(round(stat$p, 2)), "<0.01")
stat$r<-round(stat$r, 2)

ggplot()+geom_hline(yintercept=0, size=0.2)+geom_vline(xintercept=0, size=0.2)+
scale_radius(range=c(1,4))+
geom_text(data=subset(stat, n==4), aes(x=0.65, y=1, label=paste("r =",r)), size=1.5)+
geom_text(data=subset(stat, n==4), aes(x=0.65, y=0.8, label=paste("p =",p)), size=1.5)+
geom_point(data=subset(nndat, n==4), aes(x=r_rep, y=r_sum, fill=Region), shape=21, stroke=0.2)+
ggtitle("Neighbour responses")+
facet_wrap(~Region, nrow=3, strip.position="top")+
xlab(expression(Delta~"abundance of taxa"))+ylab(expression(Delta~"abundance of neighbours"))


################################
# mort vs recovery

Tot_Loss<-aggregate(Mort~Region, rd, sum)
rd$Tot_Loss<-Tot_Loss$Mort[match(rd$Region, Tot_Loss$Region)]
rd$Rel_mort<-rd$Mort/rd$Tot_Loss

ggplot()+geom_hline(yintercept=0)+
ggtitle("Mortality vs regeneration")+
geom_segment(data=rd, aes(x=Rel_mort, xend=Rel_mort,y=Replace, yend=0), size=0.1)+
geom_point(data=rd, aes(x=Rel_mort, y=Replace, fill=Region), shape=21, stroke=0.2, size=2)+
geom_text(data=subset(rd, Replace>2), aes(x=Rel_mort, y=Replace+1, label=Taxon), size=1.5)+
scale_radius(range=c(1,4))+
labs(x="Initial mortality", y="Change in % cover after recovery")


################################
# figure 4

source("R/fig_4.R")
Fig4
#ggsave("figs/Fig4.png", Fig4, width=125, height=85, units="mm", dpi=200)


