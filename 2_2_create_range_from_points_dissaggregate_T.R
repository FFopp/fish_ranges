
### =========================================================================
### Set directories and source functions
### =========================================================================

source("../raw_fish_ranges_scripts/range_fun.R")
output_dir <- file.path('../fish_ranges_output_range')
library(rgdal)
library(tools)

if(interactive()){
  ######   IF NOT CLUSTER   #########
  fromhere <- 1
  tohere <- 1
  tempfiles <- 'temp01'
  input_dir <- '../raw_fish_ranges_occ_caribbean/point_selected_dissaggregated'
  ###################################
  
  
} else {
  #########       CLUSTER   #########
  args <- commandArgs(trailingOnly = TRUE)
  fromhere <- as.numeric(args[1])
  tohere <- as.numeric(args[2])
  tempfiles <- as.character(args[3])
  input_dir <- as.character(args[4])
  ###################################
}



proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

bioregions_coastal <- readOGR(dsn='../raw_fish_ranges_shapefiles/bioregion_caribbean', layer='shape_caribbean_pacific')


### =========================================================================
### Loop over species and create maps
### =========================================================================

sps <- list.files(input_dir) 
for (id in fromhere:tohere){
  
  
  ### =========================================================================
  ### Prepare data for species id
  ### =========================================================================
  
  # Get name
  spi <- sps[id]
  spi_name <- file_path_sans_ext(spi)
  
  print(paste0("working with        ", spi_name, "         [index:", id, "]"))
  
  #spina=gsub(" ","_",spi_name)
  #sdm_file <- paste0("../fish_ranges_output_sdm/ComitteeVote_",spina,".tif")
  
 # if(file.exists(sdm_file)){
 #   print(paste0("working with        ", spi_name, "         [index:", id, "]"))
  range.fun(species_name=spi_name, occ_coord = read.table(file.path(input_dir, spi), header=T), proj=crs(bioregions_coastal), Bioreg=bioregions_coastal, Bioreg_name='Id',
                          final_resolution=0.2, degrees_outlier = 10, clustered_points_outlier = 1, Climatic_layer=NA, return_raster = T,
                          buffer_width_point=2, buffer_increment_point_line=1, buffer_width_polygon=3, dest_output=output_dir,
                          cover_threshold=0.1, dir_temp=file.path(output_dir, tempfiles), method='PCM', overwrite=T, write_raster = T, write_plot=T, desaggregate_points=0.00083)
 # } else{
 #   
 #   print(paste0("no SDM for        ", spi_name, "         [index:", id, "]"))
 # }
  
   
}

unlink(file.path(output_dir, tempfiles), recursive=T)
