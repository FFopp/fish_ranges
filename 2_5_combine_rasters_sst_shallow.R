
# stack all rasters together and calculate species richness

library(raster)
input_dir <- file.path('../fish_ranges_output_combined_sst')
#folder_names <- paste(1:6, 'models', sep='_')

raster_files <- list.files(file.path(input_dir, 'species_rasters_sst'), full.names = T, pattern='.tif')
vecm <- 1:length(raster_files)
vecms <- split(vecm, ceiling(seq_along(vecm)/100))

library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl) #register 8 nodes

foreach(i = names(vecms), .packages='raster') %dopar% {
  #for(i in names(vecms)){
  richness_stack <- stack(raster_files[vecms[[i]]])
  writeRaster(richness_stack, filename = file.path(input_dir,  'richness_rasters_sst', paste0('species_richness_stack', i, '.tif')), prj=T)
  richness_raster <- sum(richness_stack, na.rm=T)
  #richness_raster <- calc(richness_stack, sum, na.rm=T) #sum() is almost twice as fast
  writeRaster(richness_raster, filename = file.path(input_dir, 'richness_rasters_sst', paste0('species_richness', i, '.tif')), prj=T)
  #cat(i, 'done \n')
}



richness_total_stack <- stack(file.path(input_dir, 'richness_rasters_sst', paste0('species_richness', names(vecms), '.tif')))
richness_total <- sum(richness_total_stack, na.rm=T)
writeRaster(richness_total, file.path(input_dir, 'richness_rasters_sst', 'species_richness_combined_sst.tif'), prj=T)

