library(tidyverse);library(sf);library(terra);library(whitebox)
library(mapview);library(gdalUtilities);library(osmdata)

#Data from https://www.foroyakort.fo/

#Digital surface model
gdalinfo("foroyakort/FO_DSM_2017_FOTM_2M.tif")
dem <- rast("foroyakort/FO_DSM_2017_FOTM_2M.tif")

#Slope map
slope <- terrain(dem, v="slope", unit="degrees", filename="data/slope.tif")

#Vector data
st_layers("foroyakort/lendiskort.gdb/")
st_layers("foroyakort/fyrisitingarlig.gdb/")

#Island polygons
islands <- st_read("foroyakort/lendiskort.gdb/", layer = "oyggjar") |> 
  st_cast("MULTIPOLYGON") |> 
  filter(!is.na(name)) |> 
  select(name)

mapview(islands)
st_write(islands, "data/fo_islands.sqlite", delete_dsn = TRUE)

#Water polygon (including lakes)
water <- st_read("foroyakort/lendiskort.gdb/", layer = "bruksoki") |> 
  filter(area_type=="Vatn") |> 
  st_cast("MULTIPOLYGON") |> 
  select(feature_id) |> 
  mutate(area_m2 = as.numeric(st_area(Shape)),
         shoreline_m = as.numeric(st_length(st_cast(Shape, "MULTILINESTRING")))) |> 
  filter(area_m2 > 10)

#Filter out polygons that are river segments, port basins, or other non-lake polygon
water_islands <- water |> 
  st_centroid() |> 
  st_join(islands) |> 
  st_drop_geometry() |> 
  select(feature_id, name)

river_feature_ids <- c(301, 34, 6923, 6061, 6923, 7092, 772, 714, 7091, 6170, 7090, 41, 6166)
other_feature_ids <- c(559, 561, 780, 809, 787, 835, 843) ###check for more near ports

lakes <- water |> 
  left_join(water_islands) |> 
  filter(!is.na(name),
         !(feature_id %in% river_feature_ids),
         !(feature_id %in% other_feature_ids))

#Extract lake elevation
lakes_elev <- exact_extract(dem, lakes, fun="mean")

lakes_elev <- lakes |> 
  add_column(elevation_m = lakes_elev) |> 
  arrange(feature_id)

mapview(lakes_elev)
st_write(lakes_elev, "data/fo_lakes.sqlite", delete_dsn = TRUE)

#ArcticDEM Mosaic tiles
#https://www.pgc.umn.edu/data/arcticdem/
st_layers("arcticdem/ArcticDEM_Mosaic_Index_v4_1_gdb.gdb/")
tiles <- st_read("arcticdem/ArcticDEM_Mosaic_Index_v4_1_gdb.gdb/", layer="ArcticDEM_Mosaic_Index_v4_1_2m") |> 
  st_zm()
mapview(tiles)

tiles_idx <- c(7373, 7372, 4321, 7374, 4322, 3286, 4517)
tiles_sub <- tiles[tiles_idx, ]
mapview(tiles_sub)

#Copy-paste fileurl into browser and download

raw_dem_files <- list.files("arcticdem", pattern = "*_dem.tif", 
                            recursive = TRUE, full.names = TRUE)

gdalbuildvrt(raw_dem_files, "arcticdem/dem.vrt")
gdalinfo("arcticdem/dem.vrt")

#5316 is used for vector data
#for arcticdem vertical reference is height above the WGS84 ellipsoid
#Must to adjusted to geoid (orthometric height prior to use)

#egm08_25.gtx downloaded from osgeo github and moved to sf proj
#Negative værdier forekommer - anden geoid model?
gdalwarp("arcticdem/dem.vrt",
         "data/dem_20m.tif",
         s_srs = "EPSG:3413",
         t_srs= "EPSG:5316+3855",
         tr=c(20, 20),
         cutline = "data/fo_islands.sqlite",
         crop_to_cutline = TRUE,
         r="bilinear",
         co="COMPRESS=LZW",
         overwrite = TRUE)

#Openstreetmap vector data
query <- opq(bbox = c(-7.756, 61.133, -6.053, 62.441)) %>%
  add_osm_feature(key = 'natural', value = 'water')

resp <- osmdata_sf(query)

osm_natural_water <- st_cast(resp$osm_polygons, "MULTIPOLYGON") |>
  bind_rows(resp$osm_multipolygons)

st_write(osm_natural_water, "data/osm_natural_water.sqlite", delete_dsn = TRUE)
