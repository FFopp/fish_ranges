####################################################
###################   STEP 1 #######################
#         Put rasters in seperate folders          #
####################################################
####################################################
method <- c('glm', 'gam', 'gbm', 'rf')
levels <- c('complex', 'interm', 'simple')

library(raster)
input_sdm <- '../fish_ranges_output_sdm_sst'

for(i in method){
  species_modelled <- list.dirs(input_sdm, full.names = F)
  species_modelled <- species_modelled[!is.element(species_modelled, "")]
  
  for(spi in species_modelled){
    selected_models <- list.files(file.path(input_sdm, spi), full.names = T, pattern=paste0('_', i, '-'))
    
    for(j in levels){
      output_dir <- file.path(paste0('../fish_ranges_output_sdm_sst_', i), j)
      if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive = T)
      }
      mod_compexity <- paste0('-', j, '.tif')
      model <- selected_models[grepl(mod_compexity, selected_models, fixed=T)]
      file.copy(model, paste0(file.path(output_dir, spi), '.tif'))
    }
  
    
    #selected_models_stack <- stack(selected_models)
    #names(selected_models_stack) <- c('complex', 'intermediate', 'simple')
    #writeRaster(selected_models_stack, paste0(file.path(output_dir, spi), '.tif'))
  }
  
  rm(output_dir)
}

rm(input_sdm)


####################################################
###################   STEP 2 #######################
#         combine range and sdm rasters            #
####################################################
####################################################

lib_vect <- c('doParallel', 'raster', 'tools')
sapply(lib_vect, library, character.only=T)

cl <- makeCluster(12)
registerDoParallel(cl) #register 8 nodes

source("../raw_fish_ranges_scripts/combine_fun.R")
input_range <- '../fish_ranges_output_range/raster_range'

method <- c('glm', 'gam', 'gbm', 'rf')
levels <- c('complex', 'interm', 'simple')

for(i in method){
  for(j in levels){
    sps <- list.files(input_range)
    input_sdm <- file.path(paste0('../fish_ranges_output_sdm_sst_', i), j)
    output_dir <- file.path('../fish_ranges_output_combined_sst', i, j)
    if(!dir.exists(output_dir)){
      dir.create(output_dir, recursive = T)
    }
    
    foreach (spi = sps, .packages=c('tools', 'raster')) %dopar% {
      spi_name <- file_path_sans_ext(spi)
      print(paste0("working with        ", spi_name))
      
      #set name of sdm raster
      sdm_name <-  gsub(' ', '_', spi)
      
      if(file.exists(file.path(input_sdm, sdm_name))){
        sdm_raster <- raster(file.path(input_sdm, sdm_name))
        range_raster <- raster(file.path(input_range, spi))
        
        #j <- ifelse(max(values(sdm_stack), na.rm=T)>3, 4, ifelse(max(values(sdm_stack), na.rm=T)>2, 3, ifelse(max(values(sdm_stack), na.rm=T)>1, 2, 1)))
        #j <- ifelse(max(values(sdm_stack), na.rm=T)>1, 2, 1)
        #sdm_raster <- sdm_stack>=j
        combine.range.sdm.fun(species_name = spi_name, rast_sdm = sdm_raster, rast_range = range_raster, return_raster = F, write_plot = F, overwrite=T,
                              final_resolution = res(sdm_raster)[1], dest_output = file.path(output_dir, 'species_rasters'))
        
      } else {
        cat('No SDM available for', spi_name, '\n')
      }
      
    } #end sps-loop (parallelised)
  } #end j-level loop
} #end i-method loop

  
