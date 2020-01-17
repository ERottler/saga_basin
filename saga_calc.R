###

#Watershed delineation using R
#Erwin Rottler, University of Potsdam

###

pacman::p_load(RSAGA, rgdal, raster, sp, gdalUtils, maptools, leaflet, readr, htmlwidgets)

base_dir <- "/home/rottler/ownCloud/RhineFlow/rhine_obs/saga_basin/"
data_dir <- "/media/rottler/data2/basin_data/"

#set path to SAGA
env <- rsaga.env()

#prep_dem----

#Paths DEM files
path_dem_E30N20 <- paste0(data_dir, "eu_dem/eu_dem_v11_E30N20/eu_dem_v11_E30N20.TIF")
path_dem_E30N30 <- paste0(data_dir, "eu_dem/eu_dem_v11_E30N30/eu_dem_v11_E30N30.TIF")
path_dem_E40N20 <- paste0(data_dir, "eu_dem/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF")
path_dem_E40N30 <- paste0(data_dir, "eu_dem/eu_dem_v11_E40N30/eu_dem_v11_E40N30.TIF")

#load DEM files
dem_E30N20 <- raster(path_dem_E30N20)
dem_E30N30 <- raster(path_dem_E30N30)
dem_E40N20 <- raster(path_dem_E40N20)
dem_E40N30 <- raster(path_dem_E40N30)

#make template raster to build onto
x_min_min <- min(c(dem_E30N20@extent[1], dem_E30N30@extent[1], dem_E40N30@extent[1], dem_E40N20@extent[1]))
x_max_max <- min(c(dem_E30N20@extent[2], dem_E30N30@extent[2], dem_E40N30@extent[2], dem_E40N20@extent[2]))
y_min_min <- min(c(dem_E30N20@extent[3], dem_E30N30@extent[3], dem_E40N30@extent[3], dem_E40N20@extent[3]))
y_max_max <- min(c(dem_E30N20@extent[4], dem_E30N30@extent[4], dem_E40N30@extent[4], dem_E40N20@extent[4]))

my_extent <- extent(x_min_min, x_max_max, y_min_min, y_max_max)

my_template <- raster(my_extent)
projection(my_template) <- dem_E30N20@crs

path_eu_dem_25 <- paste0(data_dir, "eu_dem/processed/eu_dem_25.tif")
writeRaster(my_template, filename = path_eu_dem_25, format = "GTiff", overwrite = T)

#Merge all dem into one big raster
all_dem <- c(path_dem_E30N20, path_dem_E30N30, path_dem_E40N20, path_dem_E40N30)
mosaic_rasters(gdalfile = all_dem, dst_dataset = path_eu_dem_25, of = "GTiff")

#Calculate different resolutions
path_eu_dem_100 <- paste0(data_dir, "eu_dem/processed/eu_dem_100.tif")
gdal_translate(src_dataset = path_eu_dem_25,
               dst_dataset = path_eu_dem_100,
               tr = c(100, 100))

path_eu_dem_500 <- paste0(data_dir, "eu_dem/processed/eu_dem_500.tif")
gdal_translate(src_dataset = path_eu_dem_25,
               dst_dataset = path_eu_dem_500,
               tr = c(500, 500))

path_eu_dem_1000 <- paste0(data_dir, "eu_dem/processed/eu_dem_1000.tif")
gdal_translate(src_dataset = path_eu_dem_25,
               dst_dataset = path_eu_dem_1000,
               tr = c(1000, 1000))

#Convert to .sgrd (SAGA GIS format)
path_eu_dem_100_saga <- paste0(data_dir, "eu_dem/processed/eu_dem_100.sgrd")
rsaga.import.gdal(in.grid = path_eu_dem_100, out.grid = path_eu_dem_100_saga)

path_eu_dem_500_saga <- paste0(data_dir, "eu_dem/processed/eu_dem_500.sgrd")
rsaga.import.gdal(in.grid = path_eu_dem_500, out.grid = path_eu_dem_500_saga)

#Fill sinks
path_eu_dem_100_saga_fill <- paste0(data_dir, "eu_dem/processed/eu_dem_100_fill.sgrd")
rsaga.fill.sinks(in.dem = path_eu_dem_100_saga,
                 out.dem = path_eu_dem_100_saga_fill,
                 method = "wang.liu.2006")

path_eu_dem_500_saga_fill <- paste0(data_dir, "eu_dem/processed/eu_dem_500_fill.sgrd")
rsaga.fill.sinks(in.dem = path_eu_dem_500_saga,
                 out.dem = path_eu_dem_500_saga_fill,
                 method = "wang.liu.2006")

#prep_gauges----

dir_grdc <- "/media/rottler/data2/GRDC_DAY/" #path to grdc discharge data

file_names <- list.files(path = dir_grdc, pattern = "*.Cmd", full.names = F)
file_paths <- list.files(path = dir_grdc, pattern = "*.Cmd", full.names = T)

#Read meta information from data files

for(i in 1:length(file_paths)){
  
  print(paste(i, "of", length(file_paths)))
  
  #get rows with meta information
  meta_rows <- read_lines(file_paths[i], n_max = 32)
  meta_rows <- iconv(meta_rows, "UTF-8", "ASCII", "")
  #Name
  row_name <- meta_rows[11]
  sta_name <- substr(row_name, 26, nchar(row_name))
  #Longitude
  row_long <- meta_rows[14]
  sta_long <- substr(row_long, 24, nchar(row_long))
  #Latitude
  row_lati <- meta_rows[13]
  sta_lati <- substr(row_lati, 24, nchar(row_lati))
  #Start/End time series 
  row_seri <- meta_rows[24]
  sta_seri <- substr(row_seri, 26, nchar(row_seri)-13)
  end_seri <- substr(row_seri, 36, nchar(row_seri)-3)
  
  meta_sing <- c(sta_name, sta_lati, sta_long, sta_seri, end_seri, file_names[i])
  
  
  if(i == 1){
    
    grdc_meta <- meta_sing
    
  }else{
    
    grdc_meta <- rbind(grdc_meta, meta_sing)
    
  }
  
  
}

colnames(grdc_meta) <- c("name", "latitude", "longitude", "start_series", "end_series", "file")
rownames(grdc_meta) <- NULL
grdc_meta <- as.data.frame(grdc_meta, stringsAsFactors=FALSE)
grdc_meta$latitude   <- as.numeric(grdc_meta$latitude)
grdc_meta$longitude  <- as.numeric(grdc_meta$longitude)
grdc_meta$start_series  <- as.numeric(grdc_meta$start_series)
grdc_meta$end_series  <- as.numeric(grdc_meta$end_series)

#Leaflet map with gauging stations
map <- leaflet(grdc_meta) %>%
  # Base groups
  addProviderTiles(providers$Stamen.TerrainBackground, group = "TerrainBackground") %>%
  addProviderTiles(providers$OpenTopoMap,              group = "Open Topo") %>%
  addProviderTiles(providers$Esri.WorldImagery,        group = "WorldImagery") %>%
  # Overlay groups
  addCircleMarkers(lng = grdc_meta$longitude, lat = grdc_meta$latitude, popup = ~paste0(grdc_meta$name, 
                                                                                        ", lat: ", grdc_meta$latitude,
                                                                                        ", lon: ", grdc_meta$longitude,
                                                                                        ", start: ", grdc_meta$start_series,
                                                                                        ", end: ", grdc_meta$end_series,
                                                                                        ", file: ", grdc_meta$file),
                   labelOptions = labelOptions(noHide = T, textOnly = TRUE, direction = "bottom"), stroke = F, group = "Stations", 
                   color = "darkblue", fillOpacity = 0.7) %>%
  # Layers control
  addLayersControl(
    baseGroups = c("TerrainBackground", "Open Topo", "WorldImagery"),
    overlayGroups = c("Stations"),
    options = layersControlOptions(collapsed = FALSE)
  )

map

saveWidget(map, file=paste0(base_dir, "calc_res/grdc_map.html"))

write.table(grdc_meta, file = paste0(base_dir, "calc_res/grdc_meta.csv"), sep = ";", row.names = T, quote = T)
grdc_meta <- read.table(file = paste0(base_dir, "calc_res/grdc_meta.csv"), sep = ";", header = T)

#Select river gauge
# names_sel <- c("MELLINGEN", "BRUGG", "DOMAT/EMS", "BASEL, RHEINHALLE", "KOELN")
names_sel <- c("BASEL, RHEINHALLE")

grdc_sel <- grdc_meta[grdc_meta$name %in% names_sel, ]             
          
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gauges_84 <- SpatialPoints(data.frame(lon = grdc_sel$longitude, lat = grdc_sel$latitude), proj4string =  crswgs84)
gauges    <- spTransform(gauges_84, dem_E30N20@crs)

eu_dem_500_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_500_fill.sdat"))
plot(eu_dem_500_fill)
points(gauges)

#delineation----

eu_dem_100_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_100_fill.sdat"))
eu_dem_500_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_500_fill.sdat"))

#Create target raster: gauges get values -1
target_raster <- rasterize(gauges, eu_dem_100_fill, -1)

#Write as .tif
path_target_raster <- paste0(data_dir, "eu_dem/processed/target_raster.tif")
writeRaster(target_raster, path_target_raster, overwrite = T)

#Transform to .sgrd
path_target_raster_saga <- paste0(data_dir, "eu_dem/processed/target_raster.sgrd")
rsaga.import.gdal(in.grid = path_target_raster, out.grid = path_target_raster_saga) #transform to .sgrd

#Calculate drainage area
path_target_catch <- paste0(data_dir, "eu_dem/processed/basins/target_catch.sgrd")
rsaga.geoprocessor(lib = "ta_hydrology", module = 4, env = env, 
                   param = list(TARGET =  path_target_raster_saga,
                                ELEVATION = path_eu_dem_500_saga_fill,
                                AREA = path_target_catch,
                                METHOD = 2,
                                CONVERGE = 0.001))

#Cells in this grid contain the probability of how likely a cell belongs to the 
#Reclassify the grid, so that all pixels with a probability greater than zero considered catchment
rsaga.geoprocessor(lib = "grid_tools", module = 15, env = env, 
                   param = list(INPUT =  path_target_catch,
                                RESULT = path_target_catch,
                                METHOD = 0,
                                OLD = 0,
                                NEW = 1,
                                SOPERATOR = 4,
                                OTHEROPT = 1,
                                OTHERS = -99999))

#Convertreclassified cells to polygons
path_target_catch_shp <- paste0(data_dir, "eu_dem/processed/basins/target_catch.shp")
rsaga.geoprocessor("shapes_grid", 6, env = env, 
                   param = list(GRID = path_target_catch ,
                                POLYGONS = path_target_catch_shp,
                                CLASS_ALL = 0,
                                CLASS_ID = 1))




target_catch <-  rgdal::readOGR(path_target_catch_shp, encoding = "UTF8")

plot(eu_dem_500_fill)
points(gauges)
plot(target_catch, add = T)








target_catch <- raster(paste0(data_dir, "eu_dem/processed/basins/target_catch.sdat"))







plot(target_raster)


target_raster <- rasterize(gauges, eu_dem_500_fill, -1)

path_target_raster <- paste0(data_dir, "eu_dem/processed/target_raster.tif")
path_target_raster_saga <- paste0(data_dir, "eu_dem/processed/target_raster.sgrd")
writeRaster(target_raster, path_target_raster, overwrite = T) #wirte as .tif
rsaga.import.gdal(in.grid = path_target_raster, out.grid = path_target_raster_saga) #transform to .sgrd

rsaga.geoprocessor(lib = 'ta_hydrology', 4,
                   param = list(TARGET = path_target_raster_saga,
                                ELEVATION = path_eu_dem_500_fill_burn,
                                AREA = paste0(data_dir, "eu_dem/processed/target_catch.sgrd"),
                                METHOD = 0))

target_catch <- raster(paste0(data_dir, "eu_dem/processed/target_catch.sdat"))

plot(target_catch)

summary(target_catch@data@values)





#opt_rivers----

rhin_gdb <- paste0(data_dir, "eu_hydro/Rhine/Rhine.gdb")
danu_gdb <- paste0(data_dir, "eu_hydro/Danube/Danube.gdb")

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
rhin_gdb_list <- ogrListLayers(rhin_gdb)
danu_gdb_list <- ogrListLayers(danu_gdb)

# Read the feature class
rhin_river_net_l <- readOGR(dsn = rhin_gdb, layer = "River_Net_l")
danu_river_net_l <- readOGR(dsn = danu_gdb, layer = "River_Net_l")

rhin_river_net_l_sel <- rhin_river_net_l[rhin_river_net_l@data$STRAHLER > 4, ]
danu_river_net_l_sel <- danu_river_net_l[danu_river_net_l@data$STRAHLER > 4, ]

plot(danu_river_net_l_sel)
plot(rhin_river_net_l_sel)

#snap point to line
gauges <- snapPointsToLines(gauges, rhin_river_net_l_sel)

