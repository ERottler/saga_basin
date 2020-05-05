###

#Watershed delineation using R
#Erwin Rottler, University of Potsdam

###

pacman::p_load(RSAGA, rgdal, raster, sp, gdalUtils, maptools, leaflet, readr, htmlwidgets, rgeos)

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

#Crop DEM 500
eu_dem_100 <- raster(path_eu_dem_100)
eu_dem_500 <- raster(path_eu_dem_500)
my_ext <- extent(eu_dem_500)
my_ext_buf <- my_ext + c(+800000, -300000, +400000, -700000) #xmin, xmax, ymin, ymax

my_box <- as(my_ext_buf, 'SpatialPolygons')
eu_dem_500_crop <- raster::crop(eu_dem_500, extent(my_box))
eu_dem_500_mask <- mask(eu_dem_500_crop, my_box)
plot(eu_dem_500_mask)

#Crop DEM 100
eu_dem_100 <- raster(path_eu_dem_100)
my_ext <- extent(eu_dem_100)
my_ext_buf <- my_ext + c(+800000, -300000, +400000, -700000) #xmin, xmax, ymin, ymax

my_box <- as(my_ext_buf, 'SpatialPolygons')
eu_dem_100_crop <- raster::crop(eu_dem_100, extent(my_box))
eu_dem_100_mask <- mask(eu_dem_100_crop, my_box)
plot(eu_dem_100_mask)

#Convert NA to zero for fill DEM with 'sim_qm_of_esp'
eu_dem_100_mask@data@values[is.na(eu_dem_100_mask@data@values)] <- 0
eu_dem_500_mask@data@values[is.na(eu_dem_500_mask@data@values)] <- 0

#Write dem
path_eu_dem_100_mask <- paste0(data_dir, "eu_dem/processed/eu_dem_100_mask.tif")
writeRaster(eu_dem_100_mask, path_eu_dem_100_mask, overwrite = T)

path_eu_dem_500_mask <- paste0(data_dir, "eu_dem/processed/eu_dem_500_mask.tif")
writeRaster(eu_dem_500_mask, path_eu_dem_500_mask, overwrite = T)

#Convert to .sgrd (SAGA GIS format)
path_eu_dem_100_saga <- paste0(data_dir, "eu_dem/processed/eu_dem_100.sgrd")
rsaga.import.gdal(in.grid = path_eu_dem_100, out.grid = path_eu_dem_100_saga)

path_eu_dem_500_saga <- paste0(data_dir, "eu_dem/processed/eu_dem_500.sgrd")
rsaga.import.gdal(in.grid = path_eu_dem_500, out.grid = path_eu_dem_500_saga)

path_eu_dem_100_mask_saga <- paste0(data_dir, "eu_dem/processed/eu_dem_mask_100.sgrd")
rsaga.import.gdal(in.grid = path_eu_dem_100_mask, out.grid = path_eu_dem_100_mask_saga)

path_eu_dem_500_mask_saga <- paste0(data_dir, "eu_dem/processed/eu_dem_mask_500.sgrd")
rsaga.import.gdal(in.grid = path_eu_dem_500_mask, out.grid = path_eu_dem_500_mask_saga)

#Fill sinks
path_eu_dem_100_mask_saga_fill <- paste0(data_dir, "eu_dem/processed/eu_dem_100_mask_fill.sgrd")
rsaga.fill.sinks(in.dem = path_eu_dem_100_mask_saga,
                 out.dem = path_eu_dem_100_mask_saga_fill,
                 method = "wang.liu.2006")

path_eu_dem_500_saga_fill <- paste0(data_dir, "eu_dem/processed/eu_dem_500_fill.sgrd")
rsaga.fill.sinks(in.dem = path_eu_dem_500_saga,
                 out.dem = path_eu_dem_500_saga_fill,
                 method = "wang.liu.2006")

path_eu_dem_500_mask_saga_fill_simqm <- paste0(data_dir, "eu_dem/processed/eu_dem_500_mask_fill_simqm.sgrd")
rsaga.geoprocessor(lib = "sim_qm_of_esp", 1, env = env,
                   param = list(DEM = path_eu_dem_500_mask_saga,
                                FILLED = path_eu_dem_500_mask_saga_fill_simqm,
                                DZFILL = 0.1))

path_eu_dem_100_mask_saga_fill_simqm <- paste0(data_dir, "eu_dem/processed/eu_dem_100_mask_fill_simqm.sgrd")
rsaga.geoprocessor(lib = "sim_qm_of_esp", 1, env = env,
                   param = list(DEM = path_eu_dem_100_mask_saga,
                                FILLED = path_eu_dem_100_mask_saga_fill_simqm,
                                DZFILL = 0.1))

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
names_sel <- c("ROCKENAU-SKA")

grdc_sel <- grdc_meta[grdc_meta$name %in% names_sel, ]             
          
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gauges_84 <- SpatialPoints(data.frame(lon = grdc_sel$longitude, lat = grdc_sel$latitude), proj4string =  crswgs84)
gauges_84 <- SpatialPointsDataFrame(data.frame(lon = grdc_sel$longitude, lat = grdc_sel$latitude), 
                       data = data.frame(data =rep(-1, length(grdc_sel$longitude))), proj4string =  crswgs84)
gauges    <- spTransform(gauges_84, dem_E30N20@crs)
# st_write(gauges, paste0(data_dir, "eu_dem/processed/basins/gauges.shp"), delete_dsn = T)

eu_dem_500_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_500_fill.sdat"))
plot(eu_dem_500_fill)
points(gauges)

#buffer around gauge (optional)
gauges_buf <- buffer(gauges, width = 2000)

#delineation----

eu_dem_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_100_mask_fill_simqm.sdat"))
# eu_dem_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_100_mask_fill.sdat"))
# eu_dem_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_500_mask_fill_simqm.sdat"))

#Create target raster: gauges get values -1
target_raster <- rasterize(gauges_buf, eu_dem_fill, -1)

#Write as .tif
path_target_raster <- paste0(data_dir, "eu_dem/processed/target_raster.tif")
writeRaster(target_raster, path_target_raster, overwrite = T)

#Transform to .sgrd
path_target_raster_saga <- paste0(data_dir, "eu_dem/processed/target_raster.sgrd")
rsaga.import.gdal(in.grid = path_target_raster, out.grid = path_target_raster_saga) #transform to .sgrd

#Calculate drainage area
Sys.time()
path_eu_dem_100_mask_saga_fill_simqm <- paste0(data_dir, "eu_dem/processed/eu_dem_100_mask_fill_simqm.sgrd")
path_target_catch <- paste0(data_dir, "eu_dem/processed/basins/target_catch.sgrd")
rsaga.geoprocessor(lib = "ta_hydrology", module = 4, env = env, 
                   param = list(TARGET =  path_target_raster_saga,
                                ELEVATION = path_eu_dem_100_mask_saga_fill_simqm,
                                AREA = path_target_catch,
                                METHOD = 2,
                                CONVERGE = 0.001))
Sys.time()

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

#Convert reclassified cells to polygons
name_basin <- "rockenau"
path_target_catch_shp <- paste0(data_dir, "eu_dem/processed/basins/", name_basin, "_catch.shp")
rsaga.geoprocessor("shapes_grid", 6, env = env, 
                   param = list(GRID = path_target_catch ,
                                POLYGONS = path_target_catch_shp,
                                CLASS_ALL = 0,
                                CLASS_ID = 1))

target_catch <-  rgdal::readOGR(path_target_catch_shp, encoding = "UTF8")
plot(target_catch)

#river_network----

eu_dem_fill <- raster(paste0(data_dir, "eu_dem/processed/eu_dem_100_mask_fill_simqm.sdat"))

basin_cut_raw <- rgdal::readOGR(dsn = "/media/rottler/data2/basin_data/eu_dem/processed/basins/lobith_catch.shp")

my_box <- as(my_ext_buf, 'SpatialPolygons')
dem_cro <- raster::crop(eu_dem_fill, basin_cut_raw)
dem_sub <- mask(dem_cro, basin_cut_raw)

plot(dem_sub)

#Write as .tif
path_dem_flow_accu <- paste0(data_dir, "eu_dem/processed/dem_flow_accu.tif")
writeRaster(dem_sub, path_dem_flow_accu, overwrite = T)

#Convert to .sgrd (SAGA GIS format)
path_dem_flow_accu_saga <- paste0(data_dir, "eu_dem/processed/dem_flow_accu.sgrd")
rsaga.import.gdal(in.grid = path_dem_flow_accu, out.grid = path_dem_flow_accu_saga)

#Flow network
path_river_network <- paste0(data_dir, "eu_dem/processed/basins/river_network.sgrd")
rsaga.geoprocessor(lib = "ta_channels", module = 5, env = env, 
                   param = list(SEGMENTS =  path_river_network,
                                DEM = path_dem_flow_accu_saga,
                                THRESHOLD = 9))

river_network_raw <- rgdal::readOGR(dsn = "/media/rottler/data2/basin_data/eu_dem/processed/basins/river_network.shp")

plot(river_network_raw)



#old----



#corp DEM sub-basin area
my_ext <- extent(target_catch)
my_ext_buf <- my_ext + c(-20000, +30000, -30000, +7000)

my_box <- as(my_ext_buf, 'SpatialPolygons')
dem_cro <- raster::crop(eu_dem_fill, extent(my_box))
dem_sub <- mask(dem_cro, my_box)

#BAFU basins
basins <-  rgdal::readOGR(dsn = "/media/rottler/data2/basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp")
basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(eu_dem_fill, asText = T))

plot(dem_sub, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)
points(gauges)
plot(target_catch, add = T, border = "red3")
plot(basin_base, add = T, border = "blue3")




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

