#build composition matrix for SST and SBT
setwd('/Users/Fabian/Documents/ETH/Landscape Ecology/Paper - Fish Ranges Caribbean 2020/data_2020_04_23/raw_fish_ranges_occ_caribbean')
lib_vect <- c('raster', 'tools', 'doParallel')
sapply(lib_vect, library, character.only=T)

bathymetry <- raster('../raw_fish_ranges_rasters/depth/gb_depth_caribbean_pacific.asc')
species_traits_caribbean <- readRDS('raw/fish_traits_camille_2019/Traits_Caraibean_Fish.rds')
species_traits_pacific <- readRDS('raw/fish_traits_camille_2019/Traits_pacific_fish.rds')

#see if species that are in both databases have the same traits
sp_car_pac <- species_traits_caribbean[intersect(rownames(species_traits_caribbean), rownames(species_traits_pacific)),]
sp_pac_car <- species_traits_pacific[intersect(rownames(species_traits_caribbean), rownames(species_traits_pacific)),]
sp_car_pac_comb <- as.data.frame(cbind(rownames(sp_car_pac), rownames(sp_pac_car), sp_car_pac[,1:2], sp_pac_car[, 1:2]))

species_traits_pacific <- species_traits_pacific[!is.element(rownames(species_traits_pacific), rownames(sp_car_pac_comb)),]
species_traits_combined <- rbind(species_traits_caribbean, species_traits_pacific)

#folders <- file.path(rep(c('glm', 'gam', 'gbm', 'rf'), each=3), c('complex', 'interm', 'simple'))
method <- c('glm', 'gam', 'gbm', 'rf')
levels <- c('complex', 'interm', 'simple')

cl <- makeCluster(12)
registerDoParallel(cl)

depth_threshold <- 200
for(k in method){
  for(l in levels){
    for(i in c('sst', 'sbt')){
      input_dir <- file.path(paste0('../fish_ranges_output_combined_', i),k, l, 'species_rasters')
      output_dir <- file.path(paste0('../fish_ranges_output_combined_', i), k, l, paste0('species_rasters_', i))
      if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive = T)
      }
      
      species_names <- file_path_sans_ext(list.files(input_dir, pattern='.tif'))
      
      foreach(j = species_names, .packages=c('raster', 'tools')) %dopar% {
        
        sp_name_camille <- gsub(' ', '_', j)
        if(is.element(sp_name_camille, rownames(species_traits_combined))){ #take data from trait table
          
          min_depth <- as.numeric(species_traits_combined[rownames(species_traits_combined)==sp_name_camille, 1][[1]])
          max_depth <- as.numeric(species_traits_combined[rownames(species_traits_combined)==sp_name_camille, 2][[1]])
          cat('min =', min_depth, 'max = ', max_depth, 'for species', j, '\n')
          sp_depth_trait <- mean(c(min_depth, max_depth))
          
        } else{ #take ranges from observation data
          
          cat('No traits available for species', j, '--> using occurrence data \n')
          sp_occ <- read.table(file.path('point_selected', paste0(j, '.txt')))
          
          if(is.null(nrow(sp_occ))){
            cat('No occurrence data found! \n')
            stop()
          } 
          
          sp_depth_trait <- abs(median(extract(bathymetry, sp_occ[, c('x', 'y')]), na.rm=T))
          
        }
        
        if(i=='sst' & sp_depth_trait<depth_threshold){
          file.copy(file.path(input_dir, paste0(j, '.tif')), file.path(output_dir, paste0(j, '.tif')))
        }
        if (i=='sbt' & sp_depth_trait>=depth_threshold){
          file.copy(file.path(input_dir, paste0(j, '.tif')), file.path(output_dir, paste0(j, '.tif')))
        }
      } #end j=species (parallelised) 
    } #end i=precictor
  } #end l=levels
} #end k=mehtod

rm(list=ls())
