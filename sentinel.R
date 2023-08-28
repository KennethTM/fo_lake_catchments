library(tidyverse);library(sf);library(terra);library(tmap)

#Islands
islands <- st_read("foroyakort/lendiskort.gdb/", layer = "oyggjar") |> 
  st_cast("MULTIPOLYGON") |> 
  filter(!is.na(name)) |> 
  select(name) |> 
  st_transform(32629) |> 
  vect()

terraOptions(memfrac=0.6)

tiles <- list.files("sentinel/", full.names = TRUE, pattern="*.SAFE.zip")

#describe(tiles[1], sds=TRUE) |> View()

#True color image
tci_list <- lapply(tiles, \(x) rast(x, subds=4) |> crop(ext(islands)))
tci_merge <- merge(sprc(tci_list), na.rm=FALSE)
tci_mask <- mask(tci_merge, islands, updatevalue=255) #255 for white, 0 for black bg

writeRaster(tci_mask, "data/sentinel_tci.tif", datatype="INT1U", NAflag=NA, overwrite=TRUE)
writeRaster(tci_mask, "data/sentinel_tci.png", datatype="INT1U", NAflag=NA, overwrite=TRUE)
writeRaster(tci_mask, "data/sentinel_tci.jpeg", filetype="JPEG", datatype="INT1U", NAflag=NA, overwrite=TRUE)

#NDVI
QUANTIFICATION_VALUE <- 10000
ADD_BOA_OFFSET <- -1000

ndvi_list <- lapply(tiles, function(x){
  
  #Calculate ndvi for each tiles
  x_bands <- rast(x, subds=1, lyrs=c(3, 4)) |> 
    crop(ext(islands))
  names(x_bands) <- c("B4", "B8")
  NAflag(x_bands) <- 0
  x_reflec <- (x_bands+ADD_BOA_OFFSET)/QUANTIFICATION_VALUE
  x_reflec[x_reflec < 0] = 0
  x_ndvi <- (x_reflec[["B8"]]-x_reflec[["B4"]])/(x_reflec[["B8"]]+x_reflec[["B4"]])
  
  #Only include veg/non veg areas
  #levels(x_scl) for classes
  x_scl <- rast(x, subds=2, lyrs=9)
  x_scl_rs <- resample(x_scl, x_ndvi, method="near")
  x_ndvi_mask <- mask(x_ndvi, x_scl_rs, maskvalues=c(4, 5), inverse=TRUE)
  
  return(x_ndvi_mask)
})

ndvi_merge <- merge(sprc(ndvi_list))
ndvi_mask <- mask(ndvi_merge, islands)

writeRaster(ndvi_mask, "data/sentinel_ndvi.tif", datatype="FLT4S", NAflag=-9999, overwrite=TRUE)
