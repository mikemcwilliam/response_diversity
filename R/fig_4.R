

mort<-ggplot()+geom_hline(yintercept=0)+
ggtitle("Mortality vs regeneration")+
geom_segment(data=rd, aes(x=Rel_mort, xend=Rel_mort,y=Replace, yend=0), size=0.1)+
geom_point(data=rd, aes(x=Rel_mort, y=Replace, fill=Region), shape=21, stroke=0.2, size=2)+
geom_text(data=subset(rd, Replace>2), aes(x=Rel_mort, y=Replace+1, label=Taxon), size=1.5)+
scale_radius(range=c(1,4))+
labs(x="Initial mortality", y="Change in % cover after recovery")+
scale_colour_manual(values=cols)+
scale_fill_manual(values=cols)+
guides(col="none", size="none")+
theme_cowplot()+
theme(legend.title=element_blank(), legend.position=c(0.03, 0.15), legend.text=element_text(size=8), axis.text=element_text(size=8), axis.title=element_text(size=8), plot.title=element_text(size=8, hjust=0.5))
mort


df_poly <- data.frame(
    x=c(-Inf, -Inf, Inf),
    y=c(Inf, -Inf, -Inf))


neighbs<-ggplot()+geom_hline(yintercept=0, size=0.2)+geom_vline(xintercept=0, size=0.2)+
scale_radius(range=c(1,4))+
geom_text(data=subset(stat, n==4), aes(x=0.65, y=1, label=paste("r =",r)), size=1.5)+
geom_text(data=subset(stat, n==4), aes(x=0.65, y=0.8, label=paste("p =",p)), size=1.5)+
geom_point(data=subset(nndat, n==4), aes(x=r_rep, y=r_sum, fill=Region), shape=21, stroke=0.2)+
ggtitle("Neighbour responses")+
facet_wrap(~Region, nrow=3, strip.position="top")+
xlab(expression(Delta~"abundance of taxa"))+ylab(expression(Delta~"abundance of neighbours"))+
scale_fill_manual(values=cols)+
xlim(c(-1,1))+ylim(c(-1,1))+
guides(fill="none", size="none")+
theme_cowplot()+
theme(legend.title=element_blank(), legend.text=element_text(size=5), axis.text=element_text(size=6), axis.title=element_text(size=7), plot.title=element_text(size=8), legend.position=c(0.3, 0.6),
strip.text=element_text(size=8, colour="white"), 
 axis.ticks.length=unit(0.75, "mm"),
 strip.background=element_blank(), axis.ticks.x=element_line(),
legend.key.size=unit(3, "mm"))
neighbs


#Fig4<-plot_grid(plot_grid(bars, neighbs, nrow=1,labels=c("(a)","(b)"), label_size=12), mort, rel_heights=c(1, 1), nrow=2, labels=c("", "(c)"), label_size=12)


#Fig4<-plot_grid( neighbs, mort, nrow=2,labels=c("(a)","(b)"), label_size=10, rel_heights=c(1,2), hjust=-2)
#Fig4


Fig4<-plot_grid( neighbs, mort, nrow=1,labels=c("(a)","(b)"), label_size=10, rel_widths=c(1,2), hjust=c(0,-1))
Fig4