# stack all rasters together and calculate species richness

lib_vect <- c('raster', 'tools', 'doParallel')
sapply(lib_vect, library, character.only=T)

cl <- makeCluster(8)
registerDoParallel(cl) #register 8 nodes

method <- c('glm', 'gam', 'gbm', 'rf')
levels <- c('complex', 'interm', 'simple')

for(j in c('sst', 'sbt')){
  for(k in method){
    for(l in levels){

      input_dir <- file.path(paste0('../fish_ranges_output_combined_', j), k, l, paste0('species_rasters_', j))
      output_dir <- file.path(paste0('../fish_ranges_output_combined_', j), k, l, paste0('richness_rasters_', j))
      if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive=T)
      }

      raster_files <- list.files(input_dir, full.names = T, pattern='.tif')
      vecm <- 1:length(raster_files)
      vecms <- split(vecm, ceiling(seq_along(vecm)/100))

      foreach(i = names(vecms), .packages='raster') %dopar% {
        #for(i in names(vecms)){
        richness_stack <- stack(raster_files[vecms[[i]]])
        writeRaster(richness_stack, filename = file.path(output_dir, paste0('species_richness_stack_', i, '.tif')), overwrite=T)
        richness_raster <- sum(richness_stack, na.rm=T)
        #richness_raster <- calc(richness_stack, sum, na.rm=T) #sum() is almost twice as fast
        writeRaster(richness_raster, filename = file.path(output_dir, paste0('species_richness_', i, '.tif')), overwrite=T)
        #cat(i, 'done \n')
      }
      
      
      
      richness_total_stack <- stack(file.path(output_dir, paste0('species_richness_', names(vecms), '.tif')))
      richness_total <- sum(richness_total_stack, na.rm=T)
      writeRaster(richness_total, file.path(output_dir, paste0('species_richness_combined_', j, '.tif')), overwrite=T)
      unlink(file.path(output_dir, paste0('species_richness_', names(vecms), '.tif')), recursive=T)
      unlink(file.path(output_dir, paste0('species_richness_stack_', names(vecms), '.tif')), recursive=T)

    } #l-loop
    cat('method', k, j, 'done \n')
  } #k-loop
} # j loop