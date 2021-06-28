



################################
# sensitivity analysis


traitnum<-function(N){
sub<-traitdata[,sample(ncol(traitdata),N)]
#weights<-data.frame(t=colnames(traitdata), w=c(1,1,1,1,1,1,0.5))
#colnames(sub)
#weights$w[match(weights$t,colnames(sub))]
##############################
gower<-gowdis(sub)
pca<-pcoa(gower)
axes<-data.frame(pca1=pca$vectors[,1],pca2=pca$vectors[,2],pca3=pca$vectors[,3], pca4=pca$vectors[,4],pca5=pca$vectors[,5])
rownames(axes)<-rownames(traits)
##############################
merged<-merge(x=sitedata, y=axes, by="row.names", all.x=TRUE, all.y=FALSE)
t<-merged[,c("pca1","pca2", "pca3")]
a<-merged[,colnames(sitedata)]
rownames(t)<-merged$Row.names
rownames(a)<-merged$Row.names
fd<-melt(fdisp(a=t(a), d=dist(t,method = "euclidean"))$FDis, value.name="FDis")
cwms<-functcomp(a=as.matrix(t(a)), x=t)
fd$ID<-rownames(fd)
subdat<-merge(x=fd, y=covchange, by="ID", all.x=TRUE, all.y=TRUE)
##############################
# find original/change in cover
year_by_site<-acast(subdat, Year~IDsite.x, value.var="FDis")
whichfirst<-function(col){NonNA <- col[which(!is.na(col))]; (col*0)+NonNA[1]} 
firstval_wide<-apply(year_by_site, 2, whichfirst)
firstval<-na.omit(melt(firstval_wide, varnames=c("Yr","IDsite"), value.name="OriginalFD"))
firstval$ID<-apply(firstval[ ,c("Yr","IDsite")], 1, paste, collapse=" ")
sens<-subset(merge(x=subdat, y=firstval, by="ID", all.x=T, all.y=T), select=-c(Yr))
sens$fdchange<-(sens$FDis-sens$OriginalFD)/sens$OriginalFD*100
##############################
sensloss<-melt(subset(sens, Points==3)[,c("covchange","fdchange", "FDis","Region", "Zone")])
sensloss$Ntraits<-N
subset(sensloss, variable=="fdchange")}

reps<-100
sens1<-rbind(
do.call("rbind", replicate(reps, traitnum(3), simplify=F)),
do.call("rbind", replicate(reps, traitnum(4), simplify=F)),
do.call("rbind", replicate(reps, traitnum(5), simplify=F)),
do.call("rbind", replicate(reps, traitnum(6), simplify=F)),
do.call("rbind", replicate(reps, traitnum(7), simplify=F)))
sens1agg<-aggregate(sens1$value, by=c(Region=list(sens1$Region), Zone=list(sens1$Zone), Ntraits=list(sens1$Ntraits)), FUN=mean)
sens1agg$se<-aggregate(sens1$value, by=c(Region=list(sens1$Region), Zone=list(sens1$Zone), Ntraits=list(sens1$Ntraits)), FUN=se)$x
L_quant<-function(z){quantile(z, probs = 0.025)}
U_quant<-function(z){quantile(z, probs = 0.975)}
sens1agg$L_quant<-aggregate(sens1$value, by=c(Region=list(sens1$Region), Zone=list(sens1$Zone), Ntraits=list(sens1$Ntraits)), FUN=L_quant)$x
sens1agg$U_quant<-aggregate(sens1$value, by=c(Region=list(sens1$Region), Zone=list(sens1$Zone), Ntraits=list(sens1$Ntraits)), FUN=U_quant)$x
sens1agg$sig<-ifelse(sens1agg$U_quant>0, "NS", "S")
head(sens1agg)

head(sens1)
ntraits<-ggplot(sens1, aes(x=as.factor(Ntraits), y=value, fill=Region))+
geom_hline(yintercept=0, colour="grey")+
geom_boxplot(size=0.2, outlier.size=0.05)+
scale_fill_manual(values=c(gbrcol, polycol, caribcol))+
facet_wrap(~Zone, nrow=3)+
scale_y_continuous(limits=c(-100,25), breaks=c(-80, -40,0))+
labs(x="Number of Traits", y="Recovery deficit")+
theme(legend.title=element_blank(),legend.text=element_text(size=8), axis.text.x=element_text(size=8), axis.title.x=element_text(size=8), 
axis.text.y=element_text(size=8), axis.title.y=element_blank(), 
strip.text=element_text(size=8, face="bold"), 
strip.background=element_rect(fill="white"))
ntraits

##################################################
##################################################
##################################################
##################################################

# quality

source("maire/Maire_Quality.R")

maire<-quality_funct_space(traitdata, nbdim=7, metric="Gower", plot="quality_funct_space")
maire$meanSD

axisnum<-function(N){
	#N<-3
gower<-gowdis(traitdata, w=c(1,1,1,1,1,1,0.5))
pca<-pcoa(gower)
axes<-data.frame(pca1=pca$vectors[,1],pca2=pca$vectors[,2],pca3=pca$vectors[,3], pca4=pca$vectors[,4],pca5=pca$vectors[,5],pca6=pca$vectors[,6],pca7=pca$vectors[,7])
rownames(axes)<-rownames(traits)
##############################
merged<-merge(x=sitedata, y=axes, by="row.names", all.x=TRUE, all.y=FALSE)
t<-merged[,colnames(axes)]
t<-t[,c(1:N)]
a<-merged[,colnames(sitedata)]
rownames(t)<-merged$Row.names
rownames(a)<-merged$Row.names
fd<-melt(fdisp(a=t(a), d=dist(t,method = "euclidean"))$FDis, value.name="FDis")
cwms<-functcomp(a=as.matrix(t(a)), x=t)
fd$ID<-rownames(fd)
subdat<-merge(x=fd, y=covchange, by="ID", all.x=TRUE, all.y=TRUE)
##############################
# find original/change in cover
year_by_site<-acast(subdat, Year~IDsite.x, value.var="FDis")
whichfirst<-function(col){NonNA <- col[which(!is.na(col))]; (col*0)+NonNA[1]} 
firstval_wide<-apply(year_by_site, 2, whichfirst)
firstval<-na.omit(melt(firstval_wide, varnames=c("Yr","IDsite"), value.name="OriginalFD"))
firstval$ID<-apply(firstval[ ,c("Yr","IDsite")], 1, paste, collapse=" ")
sens<-subset(merge(x=subdat, y=firstval, by="ID", all.x=T, all.y=T), select=-c(Yr))
sens$fdchange<-(sens$FDis-sens$OriginalFD)/sens$OriginalFD*100
##############################
sensloss<-melt(subset(sens, Points==3)[,c("covchange","fdchange","Region", "Zone")])
sensloss$Naxes<-N
subset(sensloss, variable=="fdchange")}


sens2<-rbind(axisnum(2),axisnum(3),axisnum(4),axisnum(5),axisnum(6),axisnum(7))
sens2$MSD<-sens2$Naxes
sens2$MSD<-ifelse(sens2$MSD==2,maire$meanSD["m_2D"],ifelse(sens2$MSD==3,maire$meanSD["m_3D"], ifelse(sens2$MSD==4,maire$meanSD["m_4D"],ifelse(sens2$MSD==5,maire$meanSD["m_5D"], ifelse(sens2$MSD==6, maire$meanSD["m_6D"],ifelse(sens2$MSD==7, maire$meanSD["m_7D"],sens2$MSD))))))
head(sens2)


qual<-ggplotGrob(ggplot(data=subset(sens2, Region=="Polynesia" & Zone=="1-7m"))+
geom_bar(aes(x=as.factor(Naxes), y=MSD), stat="identity", fill="grey80", col="black", size=0.1, width=0.5)+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
axis.text.y=element_text(size=6), axis.title.y=element_text(size=10))+
scale_y_continuous(expand=c(0,0), position="left"))


axes<-ggplot(sens2, aes(x=Naxes, y=value, colour=Region, fill=Region, shape=Zone))+
geom_hline(yintercept=0, col="grey")+
geom_line(aes(group=interaction(Region, Zone)), col="black", size=0.2)+geom_point(colour="black", stroke=0.2)+
guides(colour="none", fill="none")+
scale_colour_manual(values=c(gbrcol, polycol, caribcol))+
scale_fill_manual(values=c(gbrcol, polycol, caribcol))+
#annotation_custom(qual, xmin=3, xmax=7.5, ymin=-110, ymax=-60)+
#facet_wrap(~Zone, nrow=3)+
scale_shape_manual(values=c(24, 21, 25))+
labs(x="Number of Axes", y="% loss of FDis after recovery")+
lims(y=c(-100,25))+
theme(axis.text=element_text(size=8), axis.title=element_text(size=8), 
strip.text=element_text(size=8, face="bold"), 
strip.background=element_rect(fill="white"), legend.position=c(0.65, 0.9), legend.title=element_blank(), legend.text=element_text(size=8),legend.key.size=unit(3, "mm"))

naxes<-plot_grid(qual, axes, nrow=2, rel_heights=c(0.25,1), align="v")


FigS1<-plot_grid(naxes, ntraits, rel_widths=c(1,1.3), labels=c("(a)","(b)"),label_size=10, nrow=1, hjust=c(0,0.7))
FigS1

#jpeg('FigS1.jpg', width=130, height=115, units="mm", res=300)
#FigS1
#dev.off()



#dis1<-gowdis(traitlist, w=c(1,1,1,1,1,1,0.5))
#df.dis<-as.matrix(dis1)
#subrd<-subset(rd, Region=="Jamaica")
#df.dis2<-df.dis[as.character(subrd$Taxon),as.character(subrd$Taxon)]
#rownames(df.dis2)==subrd$Taxon
#groups<-subrd$WinLose
#dis2<-as.dist(df.dis2)
#mod <- betadisper(dis2, groups)
#anova(mod)
#permutest(mod)
#TukeyHSD(mod,which="group",ordered=FALSE,conf.level = 0.95)
#boxplot(mod, ylab = "Distance to centroid")
#scores(mod, display = c("sites", "centroids"), choices = c(1,2))
#plot(mod)


