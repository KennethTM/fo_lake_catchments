library(tidyverse);library(sf);library(terra);library(tmap)

#DEM figure
dem <- rast("data/dem_50m.tif")
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 45, 45)

fig_dem <- tm_shape(dem)+
  tm_raster(style = "cont", palette = terrain.colors(10), 
            legend.show = TRUE, title = "Elevation")+
  tm_scale_bar(position = c("left", "bottom"))

tmap_save(fig_dem, "figures/figure_dem.png")

fig_demhill <- tm_shape(hill) +
  tm_raster(palette="-Greys", style="cont", legend.show=FALSE)+
  tm_shape(dem, raster.downsample = FALSE) +
  tm_raster(alpha = 0.5, palette = "cividis",
            legend.show = TRUE, title = "Elevation", style="cont")+
  tm_scale_bar(position = c("left", "bottom"))

tmap_save(fig_demhill, "figures/figure_demhill.png")

#NDVI figure
ndvi <- rast("data/sentinel_ndvi.tif")
ndvi_50m <- aggregate(ndvi, fact=5)

fig_ndvi <- tm_shape(ndvi_50m)+
  tm_raster(style = "cont", palette = "PRGn", breaks=c(-1, -0.5, 0, 0.5, 1), 
            legend.show = TRUE, title="NDVI", legend.reverse = TRUE)+
  tm_scale_bar(position = c("left", "bottom"))+
  tm_legend()

tmap_save(fig_ndvi, "figures/figure_ndvi.png")

#TCI figure
tci <- rast("data/sentinel_tci.tif")
tci_50m <- aggregate(tci, fact=5)

fig_tci <- tm_shape(tci_50m)+
  tm_rgb()+
  tm_scale_bar(position = c("left", "bottom"))

tmap_save(fig_tci, "figures/figure_tci.png")
