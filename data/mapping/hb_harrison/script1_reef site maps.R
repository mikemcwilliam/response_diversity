#Produce map of reefs with survey sites
zoom = data.frame(reef = unique(site.info$Reef),
                  aspect = if_else(unique(site.info$Reef) %in% c("Ashmore","Flinders", "Lihou", "Marion", "Willis"), 1.02, 
                                   if_else(unique(site.info$Reef) %in% c("Osprey", "Holmes", "Kenn", "Saumarez", "Wreck"),1.01, 1.005))) 


for( i in 1:length(unique(site.info$Reef))){
reef = zoom$reef[i]
aspect = zoom$aspect[i]

pdf(paste("figures/site_maps/Sampling sites_", reef, ".pdf", sep = ""))
data <- site.info[site.info$Reef == reef,]
print(ggplot() +
  geom_polygon(data = CS_reef, 
                     aes(group = group, 
                         x = long, 
                         y = lat),
                     fill = "grey95") +
  geom_polygon(data = CS_dryreef, 
               aes(group = group, 
                   x = long, 
                   y = lat),
               fill = "coral")  + 
  geom_polygon(data = GBR_feat, 
               aes(group = group, 
                   x = long, 
                   y = lat),
               fill = "coral")  + 
    geom_point(data = data, aes(x = X, y = Y), alpha = 1, col = "cadetblue4") +
  geom_text_repel(data = data, aes(x = X, y = Y, label = Site), alpha = 1, col = "cadetblue3") +
  coord_cartesian(ylim = c(max(data$Y)*(1/aspect), min(data$Y)*(1*aspect)),
            xlim = c(min(data$X)*(1/aspect), max(data$X)*(1*aspect))) + 
  theme_void() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_rect(
            fill = "white",
            colour = "black",
            size = 1))
) #end of print
dev.off()
}




# Reef map

data <- site.info %>%
  group_by(Reef) %>%
  summarise(sites = n()) %>%
  right_join(reef.info, by = "Reef")



pdf("figures/Map_Surveyed reefs.pdf")
print(ggplot() +
        geom_polygon(data = GBR_feat, 
                     aes(group = group, 
                         x = long, 
                         y = lat),
                     fill = "grey40") +
        geom_polygon(data = CS_reef, 
                     aes(group = group, 
                         x = long, 
                         y = lat),
                     fill = "grey40") +
        geom_point(data = data, aes(x = X, y = Y, size = sites, col = SectorRegion), alpha = .7) +
        scale_size(breaks = c(3,6,12), range=c(2,10), name = "Number of\nSurvey Sites") +
        geom_label_repel(data = data, aes(x = X, y = Y, label = Reef, col = SectorRegion), alpha = 1, size = 4,
                         point.padding = unit(.5, "lines"), 
                         segment.size = .7, position = "identity") +
        scale_color_manual(values = SectorRegionCol, guide = FALSE) +
        coord_cartesian(ylim = c(2600000, 4450000),
                      xlim = c(1900000, 3400000)) +
        theme_void() +
        theme(legend.position = c(.85,.85), legend.title=element_text(size=18), legend.text=element_text(size=16),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.background = element_rect(fill = "white", colour = "black", size = 1))
)

dev.off()


