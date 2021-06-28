
staghorn<-readPNG("data/silhouettes/STAGHORN1.png")
staghorn<-rasterGrob(staghorn, interpolate=TRUE)
digitate<-readPNG("data/silhouettes/digitate.png")
digitate<-rasterGrob(digitate, interpolate=TRUE)
encrusting<-readPNG("data/silhouettes/encrusting.png")
encrusting<-rasterGrob(encrusting, interpolate=TRUE)
laminar<-readPNG("data/silhouettes/laminar.png")
laminar<-rasterGrob(laminar, interpolate=TRUE)
massive<-readPNG("data/silhouettes/massive.png")
massive<-rasterGrob(massive, interpolate=TRUE)
solitary<-readPNG("data/silhouettes/solitary.png")
solitary<-rasterGrob(solitary, interpolate=TRUE)
submassive<-readPNG("data/silhouettes/submassive.png")
submassive<-rasterGrob(submassive, interpolate=TRUE)
tabular<-readPNG("data/silhouettes/tabular.png")
tabular<-rasterGrob(tabular, interpolate=TRUE)

m <- aggregate(pca1~Morphology, traits, mean)
m$pca2 <- aggregate(pca2~Morphology, traits, mean)$pca2
m$lab<-c(1,2,3,4, 5,6,7,8,9,10,11,12)

morphspace<-ggplot()+
annotation_custom(staghorn, xmin=m[2,2]-0.37, xmax=m[2,2]+0.1, ymin=m[2,3]-0.2, ymax=m[2,3]+0.1)+
annotation_custom(digitate, xmin=m[4,2]-0.1, xmax=m[4,2]+0.27, ymin=m[4,3]-0.15, ymax=m[4,3]+0.07)+
annotation_custom(tabular, xmin=m[12,2]-0.3, xmax=m[12,2]+0.17, ymin=m[12,3]-0.12, ymax=m[12,3]+0.18)+
annotation_custom(massive, xmin=m[9,2]-0.3, xmax=m[9,2]+0.18, ymin=m[9,3]-0.34, ymax=m[9,3]+0.17)+
annotation_custom(submassive, xmin=m[11,2]-0.11, xmax=m[11,2]+0.14, ymin=m[11,3]-0.15, ymax=m[11,3]+0.03)+
annotation_custom(solitary, xmin=m[10,2]-0.17, xmax=m[10,2]+0.08, ymin=m[10,3]-0.23, ymax=m[10,3]+0.08)+
annotation_custom(laminar, xmin=m[8,2]-0.23, xmax=m[8,2]+0.15, ymin=m[8,3]-0.11, ymax=m[8,3]+0.2)+
annotation_custom(encrusting, xmin=m[6,2]-0.28, xmax=m[6,2]+0.08, ymin=m[6,3]-0.05, ymax=m[6,3]+0.15)+
geom_point(data=m, aes(x=pca1, y=pca2),  fill="bisque", size=4, stroke=0.2, na.rm=T, shape=21)+
geom_text(data=m, aes(x=pca1, y=pca2, label=lab),size=2.5, colour="black", na.rm=T)+
geom_point(data=vects, aes(x=Dim1*-0.26, y=Dim2*-0.26), size=4)+
geom_text(data=vects, aes(x=Dim1*-0.26, y=Dim2*-0.26, label=lab), size=2, col="white", fontface="bold")+
labs(title="Coral trait space", x="PCoA axis 1", y="PCoA axis 2")+
theme_cowplot()+
theme(plot.title=element_text(size=8, hjust=0.5), axis.title=element_text(size=8), axis.text=element_text(size=8))+
xlim(-0.4, 0.5)+ylim(-0.3, 0.3)+
guides(color=guide_legend("none"),shape=guide_legend("none"),fill=guide_legend("none"), size=guide_legend("none"))
morphspace

labels <- read.csv("data/traitlabels.csv")
head(labels)
labels$pca1 <- traits$pca1[match(labels$X, traits$X)]
labels$pca2 <- traits$pca2[match(labels$X, traits$X)]

tx<-1.2
ty<-1.2
traitspace2<-ggplot()+
stat_density_2d(data=traits, aes(x=pca1, y=pca2), size=0.3, alpha=0.75, bins=6)+
geom_segment(data=labels, aes(x=pca1, y=pca2, xend=as.numeric(lab1)*tx, yend=as.numeric(lab2)*ty), size=0.15, na.rm=T, col="grey60")+
geom_point(data=traits, aes(x=pca1, y=pca2), col="red", size=0.75, alpha=0.5)+
geom_label(data=labels, aes(x=as.numeric(lab1)*tx, y=as.numeric(lab2)*ty, label=labs, size=as.numeric(size)), na.rm=T, label.padding=unit(0.1, "mm"), label.size=0, label.r=unit(0.1,"mm"))+
scale_y_continuous(breaks=c(-0.2,0,0.2), limits=c(-0.3,0.38))+
scale_x_continuous(breaks=c(-0.25, 0, 0.25, 0.5), limits=c(-0.6, 0.72))+
theme_cowplot()+
theme(axis.text=element_text(size=8), axis.title=element_text(size=8), plot.title=element_text(size=8, hjust=0.5))+
labs(x="PCoA axis 1", y="PCoA axis 2")+ggtitle('Coral trait space')+
scale_size(range=c(1.4,1.8))+guides(size="none")
traitspace2

space<-plot_grid(traitspace2+theme(axis.text.x=element_blank(), axis.title.x=element_blank()), morphspace+theme(plot.title=element_blank()), nrow=2, labels=c("(a)","(b)"),label_size=10, hjust=-1, rel_heights=c(0.88,1))
space

#############################################

slope$Time<-ifelse(slope$Points==1, "Original", ifelse(slope$Points==2, "Disturbed",ifelse(slope$Points==3, "Recovering",NA)))
slope$Time<-factor(slope$Time, levels=c("Original","Disturbed","Recovering"))
head(slope)

shift<-function(R){
dat<-subset(slope, Region==R)
dat<-dat[order(-dat$x),]
ggplot()+
geom_encircle(data=subset(dat, x>0), aes(x=pca1, y=pca2,fill=Region), alpha=0.5, col=NA)+
geom_segment(data=subset(dat, x>0), aes(x=cwm1, y=cwm2, xend=pca1, yend=pca2, col=Region))+
geom_point(data=subset(dat, x>0), aes(x=pca1, y=pca2, size=x, fill=Region), shape=21)+
scale_radius(range=c(0.01,4), limits=c(0, max(dat$x)), breaks=c(5), labels=c("= 5% cover"))+
scale_y_continuous(breaks=c(-0.2,0,0.2), limits=c(-0.25, 0.3))+
scale_x_continuous(breaks=c(-0.25, 0,0.25, 0.5), limits=c(-0.43, 0.53))+
guides(fill="none", col="none")+
scale_colour_manual(values=cols)+
scale_fill_manual(values=cols)+
panel_border(colour = "gray80", size = 0.5, linetype = 1,remove = FALSE)+
theme_cowplot()+
theme(axis.text=element_blank(), 
panel.border=element_rect(colour="black", fill=NA), plot.title=element_text(size=8), 
axis.line=element_blank(),
axis.title=element_blank(), 
 axis.ticks.length=unit(0.5, "mm"), 
 strip.text=element_text(size=7, face="bold"), strip.background=element_blank(), 
 strip.placement = "outside", 
 legend.position=c(-0.03,-0.04), 
 legend.title=element_blank(), 
 legend.text=element_text(size=6), 
 legend.background=element_blank(), 
 plot.margin = unit(c(0,0.05,0,0), "cm"))+
facet_grid(Time~Region, switch="y")}

shifts<-plot_grid(shift("GBR"),
shift("Polynesia")+theme(axis.ticks.y=element_blank(), strip.text.y=element_blank()),
shift("Jamaica")+theme(axis.ticks.y=element_blank(), strip.text.y=element_blank()),
nrow=1, rel_widths=c(1.1,1,1))
shifts

#############################################

rdspace<-function(R){
	#R="GBR"
dat<-subset(rd, Region==R)
dat$Title<-"Response diversity"
calc<-subset(rddat, Region==R)
cwmdist<-as.numeric(round(calc[1,"cwmdist"],1))
dat<-dat[order(-dat$Delta),]
dat<-subset(dat, "X1">0 | "X3">0)
maxcov<-max(dat[,c("X1","X2","X3")])
ggplot()+
geom_point(data=traits, aes(x=pca1*1.3, y=pca2*1.3), col="white",size=0,  na.rm=T)+
geom_encircle(data=subset(dat, WinLose=="Gain"),aes(x=pca1, y=pca2, fill=WinLose), col=NA,show.legend=F, alpha=0.5)+
geom_encircle(data=subset(dat, WinLose=="Loss"),aes(x=pca1, y=pca2), col="grey", fill=NA, show.legend=F,  linetype="dotted")+
geom_segment(data=subset(dat, WinLose=="Loss"),aes(x=cwm1, y=cwm2, xend=pca1, yend=pca2), size=0.2, show.legend=F, col="grey40")+
geom_segment(data=subset(dat, WinLose=="Gain"),aes(x=cwm1, y=cwm2, xend=pca1, yend=pca2, col=WinLose), size=0.6, show.legend=F, col=cols[R])+
facet_wrap(~Title, strip.position="left")+
annotate("text", x = -0.25, y = -0.4, label = cwmdist, size=2)+
geom_point(data=dat, aes(x=pca1, y=pca2, size=Delta, fill=WinLose),shape=21)+ 
scale_fill_manual(values=c(cols[R],"grey"), labels=c("= Gain", "= Loss"))+ #scale_colour_manual(values=c(col,"white"))+
labs(y="Response diversity")+
scale_y_continuous(breaks=c(-0.2,0,0.2), limits=c(-0.25, 0.3))+
scale_x_continuous(breaks=c(-0.25, 0,0.25, 0.5), limits=c(-0.43, 0.53))+
guides(size = "none", col="none", fill="none")+
scale_radius(range=c(0.01,5), limits=c(0,maxcov))+
panel_border(colour = "gray80", size = 0.5, linetype = 1,remove = FALSE)+
theme_cowplot()+
theme(panel.border=element_rect(colour="black", fill=NA), axis.text=element_blank(), 
plot.title=element_text(size=7, face="bold", hjust=0.5), 
#axis.title=element_blank(), 
axis.ticks.length=unit(0.5, "mm"), 
strip.text=element_blank(), 
axis.title=element_blank(),
strip.background=element_rect(fill="white"), 
legend.title=element_blank(), 
legend.text=element_text(size=6),
 legend.position=c(-0.05, 1.03), 
 legend.key.size=unit(3,"mm"), 
 legend.direction="horizontal",
  strip.placement = "outside",
  legend.spacing.x=unit(1,"mm"), 
axis.line=element_blank(),
legend.background=element_blank(),
 plot.margin = unit(c(0,0.05,0.15,0), "cm"))}

rdplot<-plot_grid(rdspace("GBR")+theme(strip.text=element_text(size=7, face="bold", colour="white"))+guides(fill="none"),
rdspace("Polynesia")+theme(axis.ticks.y=element_blank())+ggtitle("Response diversity")+guides(fill="none"),
rdspace("Jamaica")+guides(fill="none")+theme(axis.ticks.y=element_blank(),),
nrow=1, align="h", rel_widths=c(1.1,1,1))
rdplot

#############################################
arrow1<-readPNG("data/silhouettes/up.png")
arrow1 <-rasterGrob(arrow1, interpolate=TRUE)
arrow2<-readPNG("data/silhouettes/down.png")
arrow2 <-rasterGrob(arrow2, interpolate=TRUE)

rdtxt<-function(R){
	#R="GBR"
dat<-subset(rd, Region==R)
calc<-subset(rddat, Region==R)
cwmdist<-as.numeric(round(calc[1,"cwmdist"],2))
ggplot(data=data.frame(x=c(1:3), y=c(1:3)), aes(x=x,y=y))+
annotation_custom(arrow1, xmin=0.7, xmax=2.4, ymin=2.3, ymax=2.8)+
annotation_custom(arrow2, xmin=0.65, xmax=2.45, ymin=1.7, ymax=2.2)+
annotate("text", x=2, y=2.5, label ="  FDis=", size=2)+
annotate("label", x=3, y=2.5, label=paste(round(calc["Gain","FDis"],2),""), size=2, fill=cols[R], label.size=0.1, label.padding = unit(0.1, "lines"))+
annotate("text", x = 2, y = 2, label ="  FDis=", size=2)+
annotate("label", x=3, y=2, label=paste(round(calc["Loss","FDis"],2),""), size=2, fill="grey", label.size=0.1, label.padding = unit(0.1, "lines"))+
annotate("text", x = 2, y = 1.5, label = expression(paste(~Delta,"cwm=")), size=2)+
annotate("label", x = 3, y = 1.5, label = paste(cwmdist,""), size=2, fill="white", label.size=0.1, label.padding = unit(0.1, "lines"))+
lims(x=c(1,4), y=c(1,2.7))+
theme_void()
}

txt<-plot_grid(rdtxt("GBR"),rdtxt("Polynesia"),rdtxt("Jamaica"), nrow=1)
txt

blobs<-plot_grid(NULL, shifts, NULL, rdplot,txt, NULL,  rel_heights=c(0.05,1,0.08,0.41,0.15,0.05), nrow=6, labels=c("","(c)","","(d)","",""), label_size=10, hjust=0, vjust=0.3)



Fig3<-plot_grid(space, blobs, nrow=1, rel_heights=c(1,0.5))
Fig3



