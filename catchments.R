library(tidyverse);library(sf);library(terra);library(whitebox)
library(mapview);library(gdalUtilities);library(osmdata)
library(rmapshaper);library(exactextractr);library(nngeo)

#Load data
lakes <- st_read("data/fo_lakes.sqlite")

island_names <- unique(lakes$name)

#Preproces dem for each island
for(i in island_names){
  
  if(!dir.exists(paste0("catchments/", i))){
    dir.create(paste0("catchments/", i))
  }
  
  gdalwarp("foroyakort/FO_DSM_2017_FOTM_2M.tif",
           paste0("catchments/", i, "/dem.tif"),
           tr=c(2, 2),
           cutline = "data/fo_islands.sqlite",
           crop_to_cutline = TRUE,
           cwhere=paste0("name='", i, "'"),
           co="COMPRESS=LZW",
           overwrite = TRUE,
           srcnodata = -9999,
           dstnodata = -9999)
  
  #Whitebox tools processing
  wbt_fill_single_cell_pits(dem=paste0("catchments/", i, "/dem.tif"),
                            output=paste0("catchments/", i, "/dem_single.tif"))
  
  wbt_breach_depressions_least_cost(dem = paste0("catchments/", i, "/dem_single.tif"),
                                    output = paste0("catchments/", i, "/dem_breach.tif"),
                                    dist=25,
                                    fill=TRUE)
  
  wbt_fill_depressions_wang_and_liu(dem = paste0("catchments/", i, "/dem_breach.tif"),
                                    output = paste0("catchments/", i, "/dem_breach_fill.tif"))
  
  wbt_d8_pointer(dem = paste0("catchments/", i, "/dem_breach_fill.tif"),
                 output = paste0("catchments/", i, "/dirs.tif"))
  
  wbt_d8_flow_accumulation(input = paste0("catchments/", i, "/dirs.tif"),
                           output = paste0("catchments/", i, "/accum.tif"),
                           pntr = TRUE,
                           log = TRUE)
  
  wbt_basins(paste0("catchments/", i, "/dirs.tif"),
             paste0("catchments/", i, "/basins.tif"))
  
}

#Delineate lake catchments for each island
for(i in island_names){
  
  island_lakes <- lakes |> 
    filter(name == i)
  
  if(!dir.exists(paste0("catchments/", i, "/tmp"))){
    dir.create(paste0("catchments/", i, "/tmp"))
  }
  
  print(i)
  
  raster_template <- rast(paste0("catchments/", i, "/dirs.tif"))
  
  for(j in 1:nrow(island_lakes)){
    
    print(paste0("Lake ", j, " of ", nrow(island_lakes)))
    
    row <- island_lakes[j, ]
    feature_id <- row$feature_id
    
    site_rast <- paste0("catchments/", i, "/tmp", "/lake_", feature_id, ".tif")
    watershed_rast <- paste0("catchments/", i, "/tmp", "/watershed_", feature_id, ".tif")
    watershed_vect <- paste0("catchments/", i, "/tmp", "/watershed_", feature_id, ".shp")
    
    if(file.exists(watershed_vect)){
      print("Skipping...")
      next
    }
      
    site_rasterize <- rasterize(vect(row), raster_template, field=1)
    
    writeRaster(site_rasterize, site_rast, overwrite=TRUE)
    
    wbt_watershed(paste0("catchments/", i, "/dirs.tif"),
                  site_rast,
                  watershed_rast)
    
    wbt_raster_to_vector_polygons(watershed_rast, watershed_vect)
    
  }
  
}

#Load catchments
catchment_paths <- list.files("catchments", pattern="*.shp", recursive = TRUE, full.names = TRUE)
lake_feature_ids <- as.integer(sub("*.shp", "", sub("watershed_", "", basename(catchment_paths))))
catchment_polys <- lapply(catchment_paths, function(x){
  st_read(x, quiet=TRUE) |> 
    st_union() |> 
    st_as_sf()})
 
lake_catchments <- do.call(rbind, catchment_polys) |>
  st_make_valid() |> 
  st_cast("MULTIPOLYGON") |>
  add_column(lake_feature_id = lake_feature_ids) |> 
  arrange(lake_feature_ids)

st_write(lake_catchments, "data/catchments.sqlite", delete_layer = TRUE)

#Simplify catchments
lake_catchments_simple <- ms_simplify(lake_catchments, keep=0.1, keep_shapes = TRUE)|> 
  st_make_valid() |> 
  st_cast("MULTIPOLYGON") 

lakes_nohole <- st_remove_holes(lakes)

lake_catch_diff <- mapply(function(c, l){st_difference(c, l)}, lake_catchments_simple$x, lakes_nohole$geometry)

lake_catchments_nolake <- lake_catchments_simple
lake_catchments_nolake$x <- st_sfc(lake_catch_diff)
lake_catchments_nolake_clean <- lake_catchments_nolake %>% 
  st_set_crs(st_crs(lake_catchments_simple)) %>% 
  st_cast("MULTIPOLYGON")

#Extract mean, median, min, max elevation and slope for catchments
dem <- rast("foroyakort/FO_DSM_2017_FOTM_2M.tif")
slope <- rast("data/slope.tif")
ndvi <- rast("data/sentinel_ndvi.tif")

stat_funs <- c("min", "mean", "median", "max")

dem_stats <- exact_extract(dem, lake_catchments_nolake_clean, fun=stat_funs)
names(dem_stats) <- paste0(stat_funs, "_elevation_m")

slope_stats <- exact_extract(slope, lake_catchments_nolake_clean, fun=stat_funs)
names(slope_stats) <- paste0(stat_funs, "_slope_degrees")

ndvi_stats <- exact_extract(ndvi, lake_catchments_nolake_clean, fun=stat_funs)
names(ndvi_stats) <- paste0(stat_funs, "_ndvi")

lake_catchments_feats <- lake_catchments_nolake_clean |> 
  mutate(area_m2 = as.numeric(st_area(x))) |> 
  bind_cols(dem_stats, 
            slope_stats,
            ndvi_stats) 

lake_catchments_feats |> 
  st_write("data/catchments_simple_nolake.sqlite", delete_layer = TRUE)

#Format label for popups and write to geojson
lakes |> 
  mutate(area_m2 = sprintf(area_m2, fmt = "%0.0f"),
         shoreline_m = sprintf(shoreline_m, fmt = "%0.0f"),
         elevation_m = sprintf(elevation_m, fmt = "%0.0f")) |> 
  select(feature_id, area_m2, shoreline_m, elevation_m) |> 
  st_transform(4326) |>
  st_write("fo_lakes.geojson", delete_dsn = TRUE, layer_options="COORDINATE_PRECISION=5")

lake_catchments_feats |> 
  mutate(area_m2 = sprintf(area_m2, fmt = "%0.0f"),
         min_elevation_m = sprintf(min_elevation_m, fmt = "%0.0f"),
         mean_elevation_m = sprintf(mean_elevation_m, fmt = "%0.0f"),
         max_elevation_m = sprintf(max_elevation_m, fmt = "%0.0f"),
         mean_slope_degrees = sprintf(mean_slope_degrees, fmt = "%0.1f"),
         mean_ndvi = sprintf(mean_ndvi, fmt = "%0.2f")) |> 
  select(lake_feature_id, area_m2, min_elevation_m, mean_elevation_m, max_elevation_m,
         mean_slope_degrees, mean_ndvi) |> 
  st_transform(4326) |> 
  st_write("fo_catchments.geojson", delete_dsn = TRUE, layer_options="COORDINATE_PRECISION=5")
