#############################################################
## Description: combine sdm and range rasters
##
## Date: 2019-01-30
## Author: Fabian Fopp (fabian.fopp@usys.ethz.ch)
#############################################################


# Describe function:
# This function combines the rasters produced by the range mapping function and sdm function and creates a new
# raster with the intersection of the two. In case the two rasters have a different resolution/projection, the range raster
# will be resampled to match the sdm raster.

#describe parameters
# species_name: character string of the species name. E.g. "Anemone nemorosa"
# rast_range: raster produced by the range mapping function
# rast_sdm: raster produced by the species distribution model function
# final_resolution: determines the final resolution of the species range raster.
# write_plot: should a plot of the range be saved (T/F)?
# write_raster: should the range raster be saved (T/F)?
# dest_output: path to where all rasters, plots, log files should be saved (if write_plot/write_raster=T)
# return_raster: should the raster be returned by the funciton (T/F)
# overwrite: if the plot/raster for this species already exists, should it be overwritten (T/F)?



combine.range.sdm.fun <- function(species_name, rast_sdm, rast_range, final_resolution, write_plot=F, dest_output=paste('output', Sys.Date(), sep='_'), 
                                  write_raster=T, return_raster=T, overwrite=F){
  
  library(raster)
  
  if(class(rast_range)!='RasterLayer' | class(rast_sdm)!='RasterLayer'){
    warning('Range or sdm raster does not exist. No output file created')
    return(NULL)
  } else {
    
    if(final_resolution<res(rast_range)[1] | final_resolution<res(rast_sdm)[1]){
     # warning('Resolution of range raster (', round(res(rast_range)[1],3), ') or sdm raster (', round(res(rast_sdm)[1],3), ') is coarser than desidered final resolution (', final_resolution, ').
     #           ', round(order(final_resolution<res(rast_range)[1], final_resolution<res(rast_sdm)[1])[2],3), ' will be used as final resolution')
      final_resolution <-  res(rast_sdm)[1] #is rast_range really completely irrelevant?
    }

    rast_range <- resample(rast_range, rast_sdm)
    rast_range <- rast_range>=0.5
    rast_combined <- rast_range + rast_sdm
    #rast_combined <- rast_combined==2 #for binary rasters
    rast_combined[rast_combined<=1] <- NA #for probability rasters
    rast_combined <- rast_combined-1 #for probability rasters
    
    if(res(rast_combined)[1]<final_resolution){
      rast_combined <- aggregate(rast_combined, fact=final_resolution/res(rast_combined)[1])
      rast_combined <- rast_combined>=0.5
    }
    
    
    #plot combined map
    if(write_plot==T){
      if(!dir.exists(file.path(dest_output, 'plot_combined'))){
        dir.create(file.path(dest_output, 'plot_combined'), recursive=T)
      }
      if(file.exists(file.path(dest_output, 'plot_combined', paste0(species_name, ".png"))) & overwrite==F){
        stop('Plot already exists. Use overwrite=T to overwrite it.')
      }
      png(file.path(dest_output, 'plot_combined', paste0(species_name, ".png")), height=2000, width=3000, pointsize = 24)
      plot(rast_combined, col=c(rgb(0,0,0,0), rgb(0,1,0,1)), legend=F, main=paste0(species_name, ' - range and sdm combined'), cex.main=3, axes=F, box=F) 
      maps::map('world', add=T)
      dev.off()
    }
    
    #save species range raster
    if(write_raster==T){
      if(!dir.exists(file.path(dest_output))){ #, 'raster_combined'
        dir.create(file.path(dest_output), recursive=T)#, 'raster_combined'
      }
      if(file.exists(file.path(dest_output, paste0(species_name, ".tif"))) & overwrite==F){ #, 'raster_combined'
        stop('Raster already exists. Use overwrite=T to overwrite it.')
      }
      writeRaster(rast_combined,filename=file.path(dest_output, paste0(species_name, ".tif")), overwrite=TRUE, prj=T) #, 'raster_combined'
    }
    
    if(return_raster==T){
      return(rast_combined)
    }
    
    
    rm(rast_combined, rast_range, rast_sdm)
    cat("########### End of computation for species: ",species_name," #######################", "\n") 
  
    gc()
  
  
  }
}
