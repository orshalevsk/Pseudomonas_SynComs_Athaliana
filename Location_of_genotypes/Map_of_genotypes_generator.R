'''This on is to map the different biological material to its original sampling location'''
library(ggplot2)
library(ggmap)
library(ggsn)

local_plants_and_microbes <- read.csv2("~/ownCloud/documents/Experiments/Synthetic Pseudomonas communities on different genotypes/Full scale experiment/full scale 2018/Location_of_genotypes/plant_genotypes_and_microbes_coordinates.csv")
local_plants_and_microbes$Longitude <- as.numeric(as.character(local_plants_and_microbes$Longitude))
local_plants_and_microbes$Latitude <- as.numeric(as.character(local_plants_and_microbes$Latitude ))

# you need this API key to download the map, here it is nicely explained how to do it: https://www.littlemissdata.com/blog/maps 
ggmap::register_google(key = "AIzaSyC1tlxPEaTIgz4_80UiSm68m22TEMGCvpo")

# there are many different styles, just google it and see what you like
basicMap = get_stamenmap(bbox = c(left = 8.5, bottom = 48.25, right = 9.5, top = 48.75), 
                         extent ="device", maptype = "terrain-background", zoom = 11, 
                         color = "bw", force = FALSE)

local_plants_and_microbes$long <- local_plants_and_microbes$Longitude
local_plants_and_microbes$lat <- local_plants_and_microbes$Latitude


pdf("~/ownCloud/My papers/Syncoms_paper/Figures/Location_of_genotypes/map_genotypes_bacteria_and_Tuebingen.pdf",useDingbats=FALSE, width = 10, height = 10)
ggmap(basicMap) + 
  geom_point(data = local_plants_and_microbes, 
             aes(x = Longitude, y = Latitude, color=Category, shape=Category), 
             size =4 )+
  geom_text(data = local_plants_and_microbes[!local_plants_and_microbes$Category%in%c("bacteria_det","bacteria_ey","city") & local_plants_and_microbes$name_clear!="Tu-Wal",], aes(x = Longitude, y = Latitude, label = name), 
            size = 5, vjust = 0, hjust = -0.25,color="green4")+
  geom_text(data = local_plants_and_microbes[!local_plants_and_microbes$Category%in%c("bacteria_det","bacteria_ey","city") & local_plants_and_microbes$name_clear=="Tu-Wal",], aes(x = Longitude, y = Latitude, label = name), 
            size = 5, vjust = +0.35, hjust = +1.15,color="green4")+
  geom_text(data = local_plants_and_microbes[!local_plants_and_microbes$Category%in%c("bacteria_det","bacteria_ey","plant"),], aes(x = Longitude, y = Latitude, label = name), 
            size = 5, vjust = +2, hjust = +1,color="black")+
  geom_text(data = local_plants_and_microbes[local_plants_and_microbes$Category%in%c("bacteria_ey"),], aes(x = Longitude, y = Latitude, label = name), 
            size = 5, vjust = +2, hjust = +1,color="purple")+
  geom_text(data = local_plants_and_microbes[local_plants_and_microbes$Category%in%c("bacteria_det"),], aes(x = Longitude, y = Latitude, label = name), 
            size = 5, vjust = -1, hjust = 1,color="purple")+
  scale_colour_manual(values = c("purple","purple","black","green4")) +
  scale_shape_manual(values=c(4,4, 15, 16))+
  scale_size_manual(values=c(9,3,4,5))+
  theme(legend.position="none")+
  ylab("Latitude")+
  xlab("Longitude") +
  scalebar( dist = 5, dist_unit = "km", ,x.min = 8.5,x.max = 9.4,y.min = 48.28,y.max = 48.75,
            transform = TRUE, model = "WGS84")
dev.off()
