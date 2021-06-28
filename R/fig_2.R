
# add Trapon
trapon <- read.csv("data/literature/trapon_Moorea.csv", row.names=1)
trapon$t0 <- cover$t0[match(trapon$IDsite, cover$IDsite)]
trapon$rel_cover <- (trapon$x-trapon$t0)/trapon$t0*100

colnames(trapon) == colnames(cover)

sitesplus<-rbind(cover, trapon)
head(sitesplus)

#traponDeep<-subset(covchange, IDsite.x=="Tiahura 7-15m" & Year %in% traponSlope$Year)
#traponDeep$ID<-"traponSlope"
#traponDeep$Zone<-"15-30m"
#traponDeep$IDsite.x<-"Tiahura 15-30m"
#covchange<-rbind(covchange, traponDeep)

events=data.frame(
Event=c("S", "A", "D", "S", "B", "A", "A","B", "S","B"),
Year=c(1980, 1981, 1985, 1992-1, 1992+1,1996, 2007-1, 2007+1, 2015-1, 2015+1),
Col=c("A", "C", "D", "A", "B", "C", "C","B", "A","B"),
Region=c("Jamaica","Polynesia", "Jamaica","Polynesia","Polynesia","GBR", "Polynesia","Polynesia","GBR","GBR"))
events

sitesplus$Region <- factor(sitesplus$Region, levels=c("GBR","Polynesia", "Jamaica"))
events$Region <- factor(events$Region, levels=c("GBR","Polynesia", "Jamaica"))

points <- sitesplus[sitesplus$Points > 0,c("Points","rel_cover","Region","Year")]
points <- aggregate(rel_cover~., points, mean)
points$y <- ifelse(points$Region=="Polynesia", 25, 20)

covplot<-ggplot(sitesplus, aes(as.numeric(Year), rel_cover))+
geom_hline(yintercept=c(0,-109))+
geom_segment(data=points, aes(x=Year, xend=Year, y=rel_cover, yend=y), size=0.25, linetype="dotted")+
geom_line(aes(group=Zone),  size=0.2, col="grey")+
stat_summary(geom="line", fun="mean", size=0.75)+
stat_summary(aes(col=Region), geom="line", fun="mean", size=0.7)+
geom_point(aes(fill=Region),  size=0.15, col="black")+
facet_wrap(~Region, ncol=1)+
scale_colour_manual(values=cols)+
xlim(c(1977, 2020))+
scale_y_continuous(limits=c(-110,30), expand=c(0,0), breaks=c(0,-40,-80))+
labs(y="% change in cover", x="Year")+
geom_point(data=events, aes(x=Year, y=-100), shape=22, size=3, stroke=0.1)+
geom_point(data=points, aes(x=Year, y=y), size=3, shape=21, fill="white", stroke=0.25)+
geom_text(data=points, aes(x=Year, y=y, label=Points), size=2)+
ggtitle("Coral cover")+
geom_text(data=events, aes(x=Year, y=-100, label=Event),size=2, fontface="bold")+
guides(col="none", fill="none")+
theme_cowplot()+
theme(axis.line.x=element_blank(), axis.title=element_text(size=8),axis.text=element_text(size=8), strip.text=element_text(size=8), strip.background=element_blank(), plot.title=element_text(size=10, hjust=0.5))
covplot

maxfd<-aggregate(FDis~Region+Year+Points,subset(fd, Points>0), max)

fdplot <- ggplot(fd[fd$Points>0,], aes(as.numeric(Year), FDis))+
geom_path(aes(group=paste(Year,Region)),col="grey")+
stat_summary(aes(group=Region), fun="median", geom="line", size=0.8)+
stat_summary(fun="median", geom="line", size=0.75, aes(col=Region))+
geom_point(aes(shape=Zone, fill=Region),stroke=0.2, size=1.2)+
#stat_summary(fun="median", geom="point")+
scale_colour_manual(values=cols)+
scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
labs(y="Trait diversity (FDis)", x="Time Point")+ggtitle("Trait diversity")+
scale_shape_manual(values=c(24, 21, 25))+
geom_text(data=maxfd, aes(x=Year, y=FDis+0.012, label=Points), size=2.5)+
theme_cowplot()+
theme(axis.text=element_text(size=8), axis.title=element_text(size=8),legend.key.size=unit(2.5, "mm"), legend.title=element_blank(),  legend.text=element_text(size=7), legend.position=c(0.03, 0.14),plot.title=element_text(size=8, hjust=0.5))
fdplot

recov$variable <- factor(recov$variable, levels=c("rel_cover", "rel_fd"))
labs<-subset(recov, Region=="GBR")

recov$Zone <- factor(recov$Zone, levels=c("1-7m","7-15m","15-30m"))


rplot <- ggplot(data=recov, aes(x=variable, y=value, fill=Region, group=Zone))+
geom_hline(yintercept=0)+
geom_bar(stat="identity", position="dodge", col="black", size=0.25, width=0.6)+
guides(fill="none")+
geom_text(data=labs, aes(x=variable, y=15, label=rep(c("1-7 ","| 7-15 |",
"     15-30"),2)), position = position_dodge(width = 0.7), size=2.2)+
scale_fill_manual(values=cols)+
labs(title="Recovery deficit", y="% difference", x="")+
geom_vline(xintercept=1.5, linetype="dotted", size=0.25)+
facet_wrap(~Region, nrow=3, scales="free_y")+
scale_x_discrete(labels = c(expression(Delta~"Cover"),expression(Delta~"FDis")))+
theme_cowplot()+
theme(axis.title=element_text(size=8),legend.title=element_blank(),
strip.text=element_blank(),strip.background=element_blank(), axis.text.y=element_text(size=7), axis.text.x=element_text(size=8), legend.key.height=unit(5,"mm"), legend.key.size=unit(3, "mm"), legend.text=element_text(size=8), legend.background=element_rect(size=10, fill="grey91"), plot.title=element_text(size=8, hjust=0.5))
rplot


fig2<-plot_grid(covplot, plot_grid( fdplot, rplot, nrow=2, rel_heights=c(1,1), align="v", labels=c("(b)","(c)"), label_size=12), nrow=1, labels=c("(a)",""), rel_widths=c(1, 0.9), label_size=12)
fig2




#annotate("text", label=paste("Lizard Island, \n GBR"), x=1986, y=-60, size=1.8, col="grey30", fontface="bold"),
#annotation_custom(ggplotGrob(pmor), xmin=2006.5, xmax=2023.5, ymin=-70, ymax=-10)+
 # annotate("text", label=paste("Mo'orea, \n French Polynesia"), x=2015, y=-60, size=1.8, col="grey30", fontface="bold"),
#annotation_custom(ggplotGrob(pjam), xmin=1992, xmax=2012, ymin=-60, ymax=0)+
  #annotate("text", label=paste("Discovery Bay, \n Jamaica"), x=2002, y=-50, size=1.8, col="grey30", fontface="bold")
  
