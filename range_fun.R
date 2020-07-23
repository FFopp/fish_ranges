##################################################
## Description: create range maps based on bioregions and mark outliers
## 
## Date: 2018-06-23 15:05
## Author: Fabian Fopp (fabian.fopp@usys.ethz.ch) adapted 
## from Oskar Hagen (oskar@hagen.bio) adapted 
## from Camille Albouy (albouycamille@gmail.com)
##################################################

range.fun <- function (species_name, occ_coord, proj, Climatic_layer, Bioreg, Bioreg_name = 'ECO_NAME', final_resolution, degrees_outlier=3, clustered_points_outlier=1,
                       buffer_width_point=0.5, buffer_increment_point_line=0.5, buffer_width_polygon=0.1, cut_off=0.05, method="PCM", cover_threshold=0.3,
                       write_plot=F, write_raster=T, return_raster=T, overwrite=F, dir_temp=file.path(paste('output', Sys.Date(), sep='_'), 'temp'), desaggregate_points=F, dest_output=paste('output', Sys.Date(), sep='_')){
  
  #describe function
  # This function estimates species ranges based on occurrence data, bioregions and a climatic layer. It first deletes outliers from the 
  # observation dataset. it then creates a polygon (convex hull) with a user specified buffer around all the observations of one bioregion. 
  # If there there is only one observation in a bioregion, a buffer around this point will be created. If all points in a bioregion are on 
  # a line, the function will also create a buffer around these points but the buffer size increases with the number of points in the line.
  #PCM-method
  # In the conversion from the spatial polygons to the final range raster, every raster cell that is covered by a certain threshold (default=0.3)
  # will be marked as present. If the cell is covered by less than this threshold (but higher than 0), the cell will be marked as absence
  # if the climatic value of this cell is below the 0.05 or above the 0.95 cliamtic quantile of all observations (quantile can be changed). 
  # Otherwise the cell will be marked as present.
  #ACM-method
  # In the conversion from the spatial polygons to the final range raster, the cell will be marked as absence if the climatic value of 
  # this cell is below the 0.05 or above the 0.95 cliamtic quantile of all observations (quantile can be changed). 
  # Otherwise the cell will be marked as present.
  
  
  #describe parameters
  # species_name: character string of the species name. E.g. "Anemone nemorosa"
  # occ_coord: a dataframe containing with two columns containing the coordinates of all observations of a given species. 
  # proj: Spatial projection in which the coordinates of the occurrence data (input) are stored. The output raster of the species range will have the same projection.
  # Climatic_layer: climate raster (e.g. temperature) used to improve the distribution range (by rejecting cells with unsuitable climate).
  # Bioreg: shapefile containg different bioregions (convex hulls will be classified on a bioreg basis)
  # final_resolution: determines the final resolution of the species range raster.
  # degrees_otulier: distance threshold (degrees) for outlier classification. If the nearest minimal distance to the next point is larger than this threshold, it will be considered as an outlier. 
  # clustered_points_outlier: maximum number of points which are closer to each other than the degrees_outlier, but should still be considered as outliers.
  # buffer_width_point: buffer (in degrees) which will be applied around single observations.
  # buffer_increment_point_line: how much should the buffer be increased for each point on a line.
  # buffer_width_polygon: buffer (in degrees) which will be applied around distribution polygons (for each bioregion)
  # cut_off: quantile of temperatures (on each site) which are not considered to be suitable for the species. e.g: cut_off=0.05: lowest and highest 5% of temperature distribution
  # method: PCM (percentage cells method) method uses temperature filter only at the edge of range raster, ACM (all cells method) method for all cells.
  # cover threshold: only if method=FF. Specifies the threshold of proportion of covered area by the range polygon above which the corresponding raster will be classified as "present" without considering the temperature.
  # write_plot: should a plot of the range be saved (T/F)?
  # write_raster: should the range raster be saved (T/F)?
  # dest_output: path to where all rasters, plots, log files should be saved (if write_plot/write_raster=T)
  # return_raster: should the raster be returned by the funciton (T/F)
  # overwrite: if the plot/raster for this species already exists, should it be overwritten (T/F)?
  # dir_temp: where should the temporary text file for the convex hull be saved? (text file will be deleted again)
  # Bioreg_name: how is the slot containing the bioregion names called?
  # desaggreagte_points: should close points be desaggregated? Speeds up clustering
  
  lib_vect <- c("raster","rgdal","sp","maptools","rgbif","shape","geometry","rgeos", "FNN", 'maps', 'grDevices', 'ClusterR', 'mclust', 'ecospat', 'matrixStats')
  sapply(lib_vect,require,character.only=TRUE)
  
  if(!dir.exists(dir_temp)){
    dir.create(dir_temp, recursive=T)
  }
  
  if(class(Climatic_layer)=='RasterLayer'){
    if(final_resolution<res(Climatic_layer)[1]){
      warning('Resolution of climatic layer (', round(res(Climatic_layer)[1],3), ') is coarser than desidered final resolution (', final_resolution, ').
              ', round(res(Climatic_layer)[1],3), ' will be used as final resolution')
      final_resolution <- res(Climatic_layer)[1]
    }
  }
  
  
  # Grid making
  grd <-  SpatialGrid(GridTopology(cellcentre.offset=c( -179.5, -89.5), cellsize=c(final_resolution, final_resolution),cells.dim=c(360,180)/final_resolution), proj4string=proj)
  
  
  if(class(Climatic_layer)=='RasterLayer'){
    Climatic_layer <- projectRaster(Climatic_layer, raster(grd)) 
  }
  
  #remove duplicates
  occ_coord <- unique(occ_coord) 
  occ_coord <- na.omit(occ_coord)
  if(desaggregate_points!=F){
    occ_coord <- ecospat.occ.desaggregation(occ_coord, min.dist=desaggregate_points)
  }
  
  
  #### START INTERNAL FUNCTION ####
  conv_function <- function (x=coord_2,proj=proj){
    
    if (nrow(x)<3){ #check number of observations points in each bioregion, if <3 create point buffer
      sp_coord <- SpatialPoints(x,proj4string=proj)
      return(gBuffer(sp_coord,width=buffer_width_point+(nrow(x)-1)*buffer_increment_point_line)@polygons[[1]]) 
      
    } else {
      #test if points are on line, if yes create point buffer around points
      is_line <- 0
      for(i in 2:(nrow(x)-1)){
        dxc <- x[i,1]-x[i-1,1]
        dyc <- x[i,2]-x[i-1,2]
        dx1 <- x[i+1,1]-x[i-1,1]
        dy1 <- x[i+1,2]-x[i-1,2]  
        is_line[i-1] <- dxc*dy1-dyc*dx1
      }
      
      if(all(abs(is_line)==0)){ 
        cat('Bioreg=', g, nrow(x), 'points laying on one line. Using buffer width of ', buffer_width_point+(nrow(x)-1)*buffer_increment_point_line, '\n')
        sp_coord <- SpatialPoints(x,proj4string=proj)
        return(gBuffer(sp_coord,width=buffer_width_point+(nrow(x)-1)*buffer_increment_point_line)@polygons[[1]])
      }
      else { #if they are not on a line, create a convex hull
        
        arg_file <- paste0('QJ Fx TO ', file.path(dir_temp, 'vert.txt'))
        vert0<-convhulln(x, arg_file)
        vert1<-scan(file.path(dir_temp,'vert.txt'),quiet=T);file.remove(file.path(dir_temp,'vert.txt'))
        vert2<-(vert1+1)[-1]
        FE_vert<-row.names(x)[vert2]
        coord_conv <- x[FE_vert,]  
        
        
        coord_conv <-rbind(coord_conv,coord_conv[1,])
        P1 <- Polygons(srl=list(Polygon(coord_conv,hole=FALSE)),ID="PolygA")
        P1 <- SpatialPolygons(Srl=list(P1),proj4string=proj)
        return(gBuffer(P1,width=buffer_width_polygon)@polygons[[1]])
      }
      
      #end of ifelse  
    }# end of ifeslse
  } # end of conv_function
  
  
  plot.occ <- function(){
    if(!dir.exists(file.path(dest_output, 'log', 'plots'))){
      dir.create(file.path(dest_output, 'log', 'plots'), recursive=T)
    }
    png(file.path(dest_output, 'log', 'plots', paste0(species_name, '.png')), width=3000, height=2000)
    plot(Bioreg)
    points(occ_coord, col='red', pch=4, cex=1.5, lwd=1.2)
    points(occ_coord, col='red', pch=0, cex=1.5, lwd=1.2)
    dev.off()
  }
  
  #### END INTERNAL FUNCTION ####
  
  if (nrow(occ_coord)<=clustered_points_outlier+1){
    plot.occ()
    warning('Too few occurrences.')
    if(!dir.exists(file.path(dest_output, 'log'))){
      dir.create(file.path(dest_output, 'log'), recursive=T)
    }
    writeLines(c('Too few occurrences', species_name), con=file.path(dest_output, 'log', paste0(species_name, '.txt')), sep=' ')
    return(NULL)
  } else{ 
    
    cat("########### Start of computation for species: ",species_name," #######################", "\n") 
    
    #create distance matrix...
    mat_dist <- as.matrix(knn.dist(occ_coord, k=clustered_points_outlier))
    
    
    #mark outliers
    cond <- apply(mat_dist, 1, function(x) x[clustered_points_outlier])>degrees_outlier
    rm(mat_dist) 
    
    print(paste0(sum(cond), " outlier's from " ,nrow(occ_coord), " | proportion from total points: ", round((sum(cond)/nrow(occ_coord))*100,0), "%"))
    
    occ_coord_mod <- occ_coord[!cond,]

    
    if(nrow(occ_coord_mod)==0){
      warning('Too few occurrences within outlier threshold.')
      if(!dir.exists(file.path(dest_output, 'log'))){
        dir.create(file.path(dest_output, 'log'), recursive=T)
      }
      writeLines(c('Too few occurrences within outlier threshold for', species_name), con=file.path(dest_output, 'log', paste0(species_name, '.txt')), sep=' ')
      plot.occ()
      cat("########### End of computation for species: ",species_name," #######################", "\n") 
      return(NULL)
      
    } else{

      #correct rownames 
      rownames(occ_coord_mod) <- 1:dim(occ_coord_mod)[1]
      
      occ_points <- SpatialPoints(coords=occ_coord_mod,proj4string=proj)
      
      points_poly_dist <- suppressWarnings(gDistance(occ_points, Bioreg, byid=TRUE))  #this gives the distances of all points to the closest bioregion (0==is inside a bioregion)
      

      #option one: create buffer around polygons:
      points_poly_buffer <- colMins(points_poly_dist)
      points_poly_buffer <- max(points_poly_buffer)
      
      
      
      cat("### Projection adjustement for bioregion shapefile...", "\n") 
      unique <- unique(Bioreg@data[[Bioreg_name]])
      gc()
      
      cat("### Interscetion between occurences and bioregions ...", "\n") 
      
      library(ClusterR) #load package again (stupid cluster node 16)
      SP_dist <- list()
      a <- 0 #g= 481, 489, 529, 534, 548, 556, 558, 559
      for(g in 1:length(unique)) {
        #cat(g, '\n')
        tmp <- as(gSimplify(Bioreg[Bioreg@data[[Bioreg_name]] == unique[g],],tol=0.001,topologyPreserve=TRUE),"SpatialPolygons")
        a <- data.frame(occ_coord_mod[names(na.omit(over(occ_points, tmp, fn=NULL))),, drop=F])
        if (nrow(a)==0) {
          SP_dist[[g]] <- NA
        } else {
          #cat("g=",g,'\n')
          
          if(nrow(a)<3){
            k <- 1
            cluster_k <- kmeans(a,k)
            cluster_k$clusters <- cluster_k$cluster 
          } else {
            m_clust <- Mclust(a, verbose=F) #to determine number of clusters
            k <- m_clust$G #k=number of clusters
            
            while(k>nrow(a)-2){k <- k-1} #reduce k if necessary so that KMeans_rcpp() will run
            if(k==0){k <- 1}
            
            cluster_k <- KMeans_rcpp(a, k, num_init = 20, initializer = 'random')
            
            while(length(unique(cluster_k$clusters))<k){
              k <- k-1
              cluster_k <- KMeans_rcpp(a, k, num_init = 20, initializer = 'random')
            }
            
          }
         
          polygons_list <- list() #empty list to store all polygons later
        
          for(i in 1:k){
            #cat('cluster =', i, '\n')
            a_temp <- a[cluster_k$clusters==i,] #kmeans (with number of clusters from mcluster)
            polygons_list[[i]] <- suppressWarnings(gIntersection(gBuffer(SpatialPolygons(Srl=list(conv_function(a_temp,proj=proj))), width=0),tmp)) #zero buffer to avoid error
            polygons_list[[i]]$ID <- i
          }  
          
          
          SP_dist[[g]] <- Reduce(rbind, polygons_list)
          #SP_dist[[g]] <- suppressWarnings(gIntersection(SP_dist[[g]], tmp)) #to avoid that buffer extends to "new" polygon. should not be a problem anyway and causes script to crash
          
          if(class(SP_dist[[g]])=='SpatialCollections'){
            SP_dist[[g]] <- SP_dist[[g]]@polyobj #only keep SpatialPolygon
          }
          
        } # end of if
        
      } # end of SP_dist
      
      
      L <- SP_dist[!is.na(SP_dist)]
      
      if(length(L)==0){
        plot.occ()
        warning('No occurrences within Bioregions. Empty raster produced.')
        if(!dir.exists(paste(file.path(dest_output, 'log'), 'plots', sep='/'))){
          dir.create(paste(file.path(dest_output, 'log'), 'plots', sep='/'), recursive=T)
        }
        writeLines(c('No occurrences within Bioregions. Empty raster produced for', species_name), con=file.path(dest_output, 'log', paste0(species_name, '.txt')), sep=' ')
        
      } else{
        
       shp_species <- Reduce(rbind, L)
        
        
        if (method=='PCM'){ 
          cat("### Using Percentage Cells Method ###", "\n")
          range_raster <- rasterize(shp_species, raster(grd), getCover=T)
          range_dataframe <- as.data.frame(range_raster, xy=T)
          
          if(class(Climatic_layer)=='RasterLayer'){
            
            
            clim_all <- extract(Climatic_layer, occ_coord_mod[, 1:2]) #store all climate at observed points
            moderate_climate <- clim_all[clim_all<quantile(clim_all,probs=1-cut_off,na.rm=TRUE) & 
                                           clim_all>quantile(clim_all,probs=cut_off,na.rm=TRUE)] #store moderate temperature
            moderate_climate <- moderate_climate[!is.na(moderate_climate)]
            range_dataframe$clim <- extract(Climatic_layer, range_dataframe[, c('x', 'y')]) #extract climate of all raster cells
            
            range_dataframe$occ <- ifelse(range_dataframe$layer>=(cover_threshold), 1, 
                                          ifelse(range_dataframe$clim>=min(moderate_climate) & range_dataframe$clim<=max(moderate_climate) & range_dataframe$layer>0, 1,
                                                 ifelse(is.na(range_dataframe$clim), NA, 0)))
            
            
          } else { 
            range_dataframe$occ <- ifelse(range_dataframe$layer>=cover_threshold, 1, 0)
          }
          rast_range <- rasterFromXYZ(range_dataframe[, c('x', 'y', 'occ')], res=final_resolution)
        } 
        else if (method=='ACM') { 
          cat("### Using All Cells Method ###", "\n")
          pts_d <- do.call(rbind,lapply(shp_species@polygons,function(z){z@Polygons[[1]]@coords}))
          
          if(class(Climatic_layer)=='RasterLayer'){
            a <- extract(Climatic_layer,occ_points)
            
            test <- Climatic_layer>quantile(a,probs=cut_off,na.rm=TRUE) #why crs(Climatic_layer)=crs(a)!=crs(test)???
            test_2 <-  Climatic_layer<quantile(a,probs=1-cut_off,na.rm=TRUE ) 
            yop <- test_2+test
          }
          
          cat("### Raster making ###", "\n")
          grd <- raster(grd); grd[] <-0
          rast_cel <- unique(cellFromXY(grd,pts_d)) #TODO CHEck if this was compleatly useless. much faster!!
          data <- do.call(rbind,cellFromPolygon(grd,shp_species,weights=TRUE))[,1]
          a<-grd; a[] <- 0
          a[][unique(c(data,rast_cel))]<-1
          
          if(class(Climatic_layer)=='RasterLayer'){
            yop_2 <- crop(yop,a) 
            yop_2 <- yop_2>1
             
            rast_range <- ((yop_2 +a)==2)
          } else{
            rast_range <- a
          }
          
          
        } else{
          stop("Please choose a valid method!")
        }
        
        crs(rast_range) <- proj
        
        #plot species range
        if(write_plot==T){
          if(!dir.exists(file.path(dest_output, 'plot_range'))){
            dir.create(file.path(dest_output, 'plot_range'), recursive=T)
          }
          if(file.exists(file.path(dest_output, 'plot_range', paste0(species_name, '.png'))) & overwrite==F){
            stop('Plot already exists. Use overwrite=T to overwrite it.')
          }
          png(file.path(dest_output, 'plot_range', paste0(species_name, '.png')), height=2000, width=3000)
          plot(Bioreg, col=gray.colors(20, start = 0.5, end = 0.9, alpha=NULL))
          plot(rast_range, col=c(rgb(0,0,0,0), rgb(0,1,0,1)), legend=F, add=T)#, maxpixels=10000*6000) 
          #map('world', add=T)
          points(occ_coord_mod, pch=3, col='red', lwd=1.5, cex=1.5)
          dev.off()
        }
        
        if(write_raster==T){
          if(!dir.exists(file.path(dest_output, 'raster_range'))){
            dir.create(file.path(dest_output, 'raster_range'), recursive=T)
          }
          if(file.exists(file.path(dest_output, 'raster_range', paste0(species_name, ".tif"))) & overwrite==F){ ### add here condition that nrow(occ_coord_mod)>0
            stop('Raster already exists. Use overwrite=T to overwrite it.')
          }
          writeRaster(rast_range,filename=file.path(dest_output, 'raster_range', paste0(species_name, ".tif")), overwrite=TRUE, prj=T)
        }
        
        if(return_raster==T){
          cat("########### End of computation for species: ",species_name," #######################", "\n") 
          return(rast_range)
        }
        rm(rast_range)
      } #end else L!=list()
    }# ende else 'enough occurrences'
    
    
    gc()
  }
  }
