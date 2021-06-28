
random<-function(R){
	#R="Polynesia"
 df<-subset(rd, Region==R)
df<-subset(df, Delta>0)
df<-data.frame(Taxon=df$Taxon, pca1=df$pca1, pca2=df$pca2, 
df[sample(1:nrow(df)),c("Delta","WinLose")])
delta<-acast(WinLose~Taxon, data=df, value.var="Delta")
delta[is.na(delta)]<-0
t<-data.frame(df[,c("pca1","pca2")], row.names=df$Taxon)
fd<-data.frame(FDis=as.matrix(fdisp(a=delta, dist(t))$FDis))
fd$WinLose<-rownames(fd)
cwm<-functcomp(a=delta, x=t)
a2<-(cwm["Gain","pca1"]-cwm["Loss","pca1"])^2
b2<-(cwm["Gain","pca2"]-cwm["Loss","pca2"])^2
cwmdist<-sqrt(a2+b2)
fd$cwm<-c(cwmdist, NA)
fd$Region<-R
fd}

n.iter=1000
null<-rbind(
do.call(rbind, replicate(n.iter, random("GBR"), simplify=F)),
do.call(rbind, replicate(n.iter, random("Polynesia"), simplify=F)),
do.call(rbind, replicate(n.iter, random("Jamaica"), simplify=F)))
null$Region<-factor(null$Region, levels=c("GBR","Polynesia","Jamaica"))
nrow(null)
head(null)

nullmean<-aggregate(null$FDis, by=c(WinLose=list(null$WinLose), Region=list(null$Region)), FUN=mean)
nullse<-aggregate(null$FDis, by=c(WinLose=list(null$WinLose), Region=list(null$Region)), FUN=se)
nullmean$se<-nullse$x


dod1<-0.5
dod2<-0.5
disp<-ggplot()+
#geom_errorbar(data=nullmean, aes(x=WinLose, ymax=x+se, ymin=x-se, col=Region), width=0, position=position_dodge(dod1))+
#geom_path(data=nullmean, aes(x=WinLose, y=x,col=Region, group=Region), linetype="dotted", position=position_dodge(dod))+
geom_boxplot(data=null, aes(x=WinLose, y=FDis, fill=Region), fatten=0, position=position_dodge(dod1), size=0.2, width=0.2, outlier.size=0, outlier.color="white")+
geom_path(data=rddat, aes(x=WinLose, y=FDis, col=Region, group=Region), position=position_dodge(dod2), linetype="dotted", size=0.75)+
geom_point(data=rddat, aes(x=WinLose, y=FDis, fill=Region), size=1.2, shape=4, stroke=1, position=position_dodge(dod2))+
guides(col=F, fill=F)+
ylab("Trait Diversity (FDis)")+
#facet_wrap(~Region, nrow=3)+
scale_fill_manual(values=c(gbrcol,polycol,caribcol))+
scale_colour_manual(values=c(gbrcol,polycol,caribcol))+
theme(axis.title.x=element_text(size=8, colour="white"), 
axis.title.y=element_text(size=8), axis.text=element_text(size=8))

head(null)
head(rddat)
cwmplot<-ggplot()+
geom_boxplot(data=null, aes(y=cwm, x="cwm", fill=Region), fatten=0, position=position_dodge(dod1), size=0.2, width=0.2, outlier.size=0, outlier.color="white", na.rm=T)+
geom_point(data=rddat, aes(x="cwm", y=cwmdist, fill=Region), size=1.2, shape=4, stroke=1, position=position_dodge(dod2))+
guides(col=F, fill=F)+
ggtitle(paste("Null models of \n winner/loser similarity"))+
labs(y="Community weighted distance")+
scale_fill_manual(values=c(gbrcol,polycol,caribcol))+
scale_colour_manual(values=c(gbrcol,polycol,caribcol))+
theme( axis.title.x=element_text(size=2, colour="white"), axis.text.x=element_text(size=2, colour="white"),
axis.title.y=element_text(size=8), axis.text.y=element_text(size=8), 
plot.title=element_text(size=8))



bars<-plot_grid(cwmplot, disp, nrow=2)
bars

(nrow(subset(null, Region=="GBR" & cwm < 2.401))/n.iter)*100
(nrow(subset(null, Region=="Polynesia" & cwm < 0.303))/n.iter)*100
(nrow(subset(null, Region=="Jamaica" & cwm < 2.2618))/n.iter)*100



rd2<-subset(rd, Taxon!="Mussidae") # SCALE ~200m
rd2<-subset(rd2, Original > 1) # SCALE ~200m
#rd2<-rd

head(rd2)
rd2$logRatio<-log10((rd2$Recovering+1)/(rd2$Disturbed+1))
#rd2$total.fecun<-log(10^(rd2$polypden)*10^(rd2$fecun))
rd2$total.fecun<-(rd2$polypden)*(rd2$fecun)


ggplot(data=rd2, aes(x=(total.fecun), y=Replace))+
#geom_text(aes(label=Taxon, col=Region))+
geom_point(aes(fill=Region), shape=21)+guides(col="none")+
scale_colour_manual(values=c(gbrcol, polycol, caribcol))+
scale_fill_manual(values=c(gbrcol, polycol, caribcol))+
geom_smooth(method="lm", formula=y~poly(x,1))

lm1<-lm(Replace~total.fecun, data=rd2)
summary(lm1)

# quantify rd somehow... 

library("ade4")

disims<-function(R){
	#R="Polynesia"
dat<-subset(rd, Region==R)
rownames(dat)<-dat$Taxon
# response difference
rdist<-dist(data.frame(x=dat$relReplace,row.names=dat$Taxon))
rdmat<-as.matrix(rdist)
rdmat[upper.tri(rdmat, diag=T)]<-NA
rlong<-melt(rdmat)
hist(sqrt(rlong$value))
# trait dissimilarity
tdist<-dist(data.frame(pca1=dat$pca1,pca2=dat$pca2,row.names=dat$Taxon))
tmat<-as.matrix(tdist)
tmat[upper.tri(tmat, diag=T)]<-NA
tlong<-melt(tmat, value.name="tdist")
hist(tlong$tdist)
# check
rlong$Var1==tlong$Var1
rlong$Var2==tlong$Var2
# combine
both<-cbind(tlong, rd=rlong$value)
both$rd<-(both$rd)  /max(both$rd, na.rm=T)
both$tdist<-both$tdist/max(both$tdist, na.rm=T)
both$Region<-R
both<-both[!is.na(both$tdist),]
# max/min models
both$t_max<-both$tdist[order(both$tdist)]
both$r_max<-both$rd[order(both$rd)]
both$t_min<-both$tdist[order(both$tdist)]
both$r_min<-both$rd[order(both$rd, decreasing=T)]
# mantel test
mant<-mantel.rtest(tdist, (rdist), nrepet=999)
both$mantR<-round(mant[1]$obs,2)
both$mantP<-round(mant[5]$pvalue,2)
both}


new<-rbind(disims("GBR"),disims("Polynesia"),disims("Jamaica"))
new$Region<-ifelse(new$Region=="GBR","Lizard Island",ifelse(new$Region=="Polynesia","Mo'orea","Jamaica"))
new$Region<-factor(new$Region, levels=c("Lizard Island","Mo'orea","Jamaica"))
head(new)

random<-function(R){
	#R="Polynesia"
 n.iter=10
dat<-disims(R)
df<-dat[,c("rd","tdist")]
rdf<- melt(replicate(n.iter, df[sample(1:nrow(df)),"rd"]), value.name="rd")
rdf$tdist<-rep(df$tdist, n.iter)
rdf$Region<-R
rdf}

null<-rbind(random("GBR"),random("Polynesia"),random("Jamaica"))
null$Region<-factor(null$Region, levels=c("GBR","Polynesia","Jamaica"))
head(null)


redun1<-ggplot(data=new, aes(x=tdist, y=rd))+
geom_point(data=new,aes(col=Region),shape=4, size=0.2)+
#geom_text(data=new,aes(col=Region, label=Var2), size=3)+
#geom_smooth(method="lm", formula=y~poly(x,2), se=F)+
#stat_density2d(bins=7, geom="polygon", aes(fill=Region), size=3)+
stat_density2d(bins=7,  aes(col=Region), size=0.5)+
scale_colour_manual(values=c(gbrcol,polycol,caribcol))+
scale_fill_manual(values=c(gbrcol,polycol,caribcol))+
facet_wrap(~Region, nrow=1, scale="free_y")+
#stat_density2d(data=null, aes(x=tdist, y=rd), bins=7, col="grey")+
geom_text(data=unique(new[,c("Region","mantR", "mantP")]), aes(x=0.5, y=1.1, label=paste("Mantel r =", mantR),group=Region), size=3)+
#coord_cartesian(xlim=c(0,1), ylim=c(0,1.1))+xlim(c(-1,2))+ylim(c(-1,2))+
labs(x="Trait disimilarity", y="Response difference")+
guides(col="none", fill="none")+theme(legend.title=element_blank(), legend.position=c(-0.03, 0.15), legend.text=element_text(size=8), axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10), 
strip.background=element_blank())
redun1


nullplot<-ggplot(data=new, aes(tdist, y=rd, col=Region))+
stat_density2d(bins=7, size=0.5)+
stat_density2d(data=null, aes(x=tdist, y=rd), bins=7, col="grey")+
facet_wrap(~Region)
nullplot







