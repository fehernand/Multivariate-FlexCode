
library("ggmap")
library(maptools)
library(maps)
require(slam)


###########
load("RDA/predtwitter.rda")
load("RDA/data_twitter.rda")

obs = 14

CDE = matrix(pred$CDE[[obs]], ncol = 1000, nrow = 1000)

redim = seq(1,1000,5)

CDE = CDE[redim,redim] 
z1 = pred$z[[1]][redim]
z2 = pred$z[[2]][redim]

library(reshape2)
  
data.plot = melt(CDE)

data.plot$Var1 = z1[data.plot$Var1]
data.plot$Var2 = z2[data.plot$Var2]

zoom = 50

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world) +
  geom_tile(data = data.plot %>% dplyr::filter(value > 0.0004), aes(Var2,Var1, fill=value),alpha = 1) +
  scale_fill_gradientn(colours=rev(rainbow(100, start=0, end=0.75)))+
  geom_point(aes(x = zTest[obs,2], y = zTest[obs,1]), size = 1.5) +
  scale_x_continuous(limits = c(zTest[obs,2] - zoom,zTest[obs,2]+ zoom)) + scale_y_continuous(limits = c(zTest[obs,1] - zoom,zTest[obs,1]+ zoom))

