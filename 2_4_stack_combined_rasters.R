
# stack all rasters together and calculate species richness


library(raster)
input_dir <- file.path('../fish_ranges_output_combined')
folder_names <- paste(1:6, 'models', sep='_')


library(doParallel)
cl <- makeCluster(detectCores())/2
registerDoParallel(cl) #register 8 nodes

foreach(i = folder_names, .packages='raster') %dopar% {
#for(i in folder_names){
  raster_files <- list.files(file.path(input_dir, i), full.names = T, pattern='.tif')
  richness_stack <- stack(raster_files)
  writeRaster(richness_stack, filename = file.path(input_dir, paste0('species_richness_stack', i, '.tif')), prj=T)
  richness_raster <- sum(richness_stack, na.rm=T)
  #richness_raster <- calc(richness_stack, sum, na.rm=T) #sum() is almost twice as fast
  writeRaster(richness_raster, filename = file.path(input_dir, paste0('species_richness', i, '.tif')), prj=T)
}

