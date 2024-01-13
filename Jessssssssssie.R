#---
#title: "UpperMerced_DamHeights_ForHarshad"
#$output: html_document
#date: "2024-01-13"
---
  
library(tidyverse)
library(ggplot2)
library(dplyr)
library(dataRetrieval)
library(ggsn)
library(sf)
library(maps)
library(sp)
library(ncdf4)
library(sf)
library(sp)
library(lwgeom)
library(terra)
library(raster)
library(mapview)

#funky projection
wkt <- 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],
          AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,
            AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

#Beaver dam volumes by 1x1km grid across the entire domain, based on BRAT model
vol_by_grid <- readRDS("vol_by_grid.rds")
rvol <- vol_by_grid
mapview(rvol)

#Make it a raster
rvol <- rast(ncols = 174, nrows = 63, xmin = -120.985, xmax = -119.255, ymin = 37.325, ymax = 37.925)
rvol
values(rvol) <- vol_by_grid$total_vol_in_cell 
plot(rvol,1, colNA="grey70") #hmm upside down
mapview(raster(rvol))
#the raster is upside down so I flipped it. The final watershed should look like the one mapped below. 
rvol  <- flip(rvol, direction="vertical")
plot(rvol,1, colNA="grey70")
mapview(raster(rvol))

#For now, here's the netcdf file: 
inputs <- nc_open("lis_input.nc")
print(inputs)
#Look at the variables 
attributes(inputs$var)

#Look at the dimensions (lat and lon)
attributes(inputs$dim)

#In this chunk I am wrangling the netcdf file so I can pull a single variable into an array. 

#Get latitude and longitude with the ncvar_get function and store each into their own object:
lat <- ncvar_get(inputs, "lat")
nlat <- dim(lat) #to check it matches the metadata: 63
lon <- ncvar_get(inputs, "lon")
nlon <- dim(lon) #to check, should be 174
# Check your lat lon dimensions match the information in the metadata we explored before:
print(c(nlon, nlat))

#Lat and Lon are both stored as two dimensional data arrays. Let's select just the unique lats and longs from those arrays: 
lat <- lat[1,]
lon <- lon[,1]

HYMAP_grid_area <- ncvar_get(inputs, "HYMAP_grid_area")
dim(HYMAP_grid_area) #to check; this should give you 174 and 63

#Pull elevation for fun just to see what elevation looks like: 
elevation <- ncvar_get(inputs, "ELEVATION")
image(lon, lat, elevation)
summary(elevation)

HYMAP_river_height <- ncvar_get(inputs, "HYMAP_river_height")
image(lon,lat,HYMAP_river_height)


#The purpose of this code is to use the beaver dam volumes (the rvol raster from above) 
#along with riverbed geometry (river height, width, length, floodplain height) for each grid cell to 
#back-calculate the height a beaver dam would need to be to store the  appropriate volume of water. 
#Then we gotta add that height to elevation to get total elevation of dam above sea level. 

#First, I pull arrays for individual variables from the netcdf file: 
HYMAP_grid_area <- ncvar_get(inputs, "HYMAP_grid_area")
HYMAP_floodplain_height <- ncvar_get(inputs, "HYMAP_floodplain_height")
HYMAP_river_width <- ncvar_get(inputs, "HYMAP_river_width")
HYMAP_river_length <- ncvar_get(inputs, "HYMAP_river_length")
HYMAP_river_height <- ncvar_get(inputs, "HYMAP_river_height")
elevation <- ncvar_get(inputs, "ELEVATION")

dim(damvol.df)


#Next, I rename my raster just for fun: 
damvol <- rvol
mapview(raster(rvol))
damvol.df <- as.matrix(t(raster(damvol)), xy=T) #Converted your raster to a matrix. Using the raster as-is is problematic as
#it uses the whole thing for every iteration of the loop instead of the specific cell you're interested in
dim(damvol.df)
str(damvol.df)
damvol.df


damheightsasl_array <- matrix(NA_real_, 174, 63) #Using a matrix seems easier than an array for some reason
floodheight_array <- matrix(NA_real_, 174, 63)
floodvol_array <- matrix(NA_real_, 174, 63)
rivervol_array <- matrix(NA_real_, 174, 63)


for (i in 1:length(lon)){ #iterate over every longitude
  for (j in 1:length(lat)){ #iterate over every latitude
    damvol_val = damvol.df[i,j] #Using the matrix instead of the raster in your original code
    rivervol = HYMAP_river_width[i,j] * HYMAP_river_length[i,j] * HYMAP_river_height[i,j] #calculate the volume of the river in that particular cell (width*length*height)
    floodvol = damvol_val - rivervol #calculate the volume of water left in the floodplain                             
    area = HYMAP_grid_area[i,j] #pull in area
    height = HYMAP_floodplain_height[i,j,] #pull in floodplain height, which is a 3 dim variable 
    volume = area * height #floodplane volume 
    model = lm(volume ~ height) #make a lm
    m = coef(model)[2]
    b = coef(model)[1]
    floodheight = (floodvol-b)/m #back calculate floodplain height
    floodheight[floodheight<0] = 0
    damheightsasl_array[i,j] = floodheight + elevation[i,j] #add to the elevation variable 
    floodheight_array[i,j] = floodheight
    floodvol_array[i,j] = floodvol
    rivervol_array[i,j] = rivervol
    
 }} 



damheightsasl_array <- t(damheightsasl_array) #Fixing the flipping and rotation issue
damheight <- rast(damheightsasl_array)
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
damheight <- raster(damheight) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(damheight) <- e
crs(damheight) <- wkt
mapview(damheight)


floodheight_array <- t(floodheight_array) #Fixing the flipping and rotation issue
floodheight_r <- rast(floodheight_array)
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
floodheight_r <- raster(floodheight_r) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(floodheight_r) <- e
crs(floodheight_r) <- wkt
mapview(floodheight_r)


floodvol_array <- t(floodvol_array) #Fixing the flipping and rotation issue
floodvol_r <- rast(floodvol_array)
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
floodvol_r <- raster(floodvol_r) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(floodvol_r) <- e
crs(floodvol_r) <- wkt
mapview(floodvol_r)


rivervol_array <- t(rivervol_array) #Fixing the flipping and rotation issue
rivervol_r <- rast(rivervol_array)
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
rivervol_r <- raster(flip(rivervol_r, direction = 'vertical')) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(rivervol_r) <- e
crs(rivervol_r) <- wkt
mapview(rivervol_r)


HYMAP_river_length_r <- rast(t(HYMAP_river_length))
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
HYMAP_river_length_r <- raster(flip(HYMAP_river_length_r, direction = 'vertical')) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(HYMAP_river_length_r) <- e
crs(HYMAP_river_length_r) <- wkt
mapview(HYMAP_river_length_r)


HYMAP_river_width_r <- rast(t(HYMAP_river_width))
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
HYMAP_river_width_r <- raster(flip(HYMAP_river_width_r, direction = 'vertical')) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(HYMAP_river_width_r) <- e
crs(HYMAP_river_width_r) <- wkt
mapview(HYMAP_river_width_r)


HYMAP_river_height_r <- rast(t(HYMAP_river_height))
e <- extent(c(-120.985, -119.255, 37.325, 37.925))
HYMAP_river_height_r <- raster(flip(HYMAP_river_height_r, direction = 'vertical')) #Because I was lazy and didn't want to find out the equivalent for doing this with a spatraster :P
extent(HYMAP_river_height_r) <- e
crs(HYMAP_river_height_r) <- wkt
mapview(HYMAP_river_height_r)



