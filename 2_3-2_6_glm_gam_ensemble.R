#make a gam_glm ensemble

####################################################
###################   STEP 0 #######################
#                 COMMON SETUP                     #
####################################################
####################################################
lib_vect <- c('doParallel', 'raster', 'tools', 'rgdal', 'plyr')
sapply(lib_vect, library, character.only=T)

cl <- makeCluster(12)
registerDoParallel(cl) #register 8 nodes

method <- 'glm_gam'
levels <- 'simple'



for(k in c('sbt', 'sst')){
  for(j in levels){
    input_sdm <- file.path(paste0('../fish_ranges_output_sdm_', k))
    sps <- list.dirs(input_sdm)
    sps <- sps[-1]
    output_glm <- file.path(paste0('../fish_ranges_output_sdm_', k, '_glm'), j)
    output_gam <- file.path(paste0('../fish_ranges_output_sdm_', k, '_gam'), j)
    if(!dir.exists(output_glm)){
      dir.create(output_glm, recursive = T)
    }
    if(!dir.exists(output_gam)){
      dir.create(output_gam, recursive = T)
    }
    for(i in sps){
      input_raster_gam <- list.files(i, pattern='gam-simple', full.names = T)
      input_raster_glm <- list.files(i, pattern='glm-simple', full.names = T)
    file.copy(from=input_raster_gam, to=output_gam)
    file.copy(from=input_raster_glm, to=output_glm)
    }
  }
}






####################################################
###################   STEP 1 #######################
#        2_3 - combine range and sdm rasters       #
####################################################
####################################################


source("../raw_fish_ranges_scripts/combine_fun.R")
input_range <- '../fish_ranges_output_range/raster_range'


for(k in c('sbt', 'sst')){
  for(j in levels){
    sps <- list.files(input_range)
    input_sdm_glm <- file.path(paste0('../fish_ranges_output_sdm_', k, '_glm'), j)
    input_sdm_gam <- file.path(paste0('../fish_ranges_output_sdm_', k, '_gam'), j)
    output_dir <- file.path(paste0('../fish_ranges_output_combined_', k), method, j)
    if(!dir.exists(output_dir)){
      dir.create(output_dir, recursive = T)
    }
    
    foreach (spi = sps, .packages=c('tools', 'raster')) %dopar% {
      spi_name <- file_path_sans_ext(spi)
      print(paste0("working with        ", spi_name))
      
      #set name of sdm raster
      sdm_name <-  gsub(' ', '_', spi_name)

      sdm_name_glm <- paste0('BinaryPredictions_', sdm_name, '_PA_geographic_w1_x_0.2env_strata_Mod_glm-simple.tif')
      sdm_name_gam <-  paste0('BinaryPredictions_', sdm_name, '_PA_geographic_w1_x_0.2env_strata_Mod_gam-simple.tif')
      
      if(file.exists(file.path(input_sdm_glm, sdm_name_glm)) & file.exists(file.path(input_sdm_gam, sdm_name_gam))){
        sdm_raster_glm <- raster(file.path(input_sdm_glm, sdm_name_glm))
        sdm_raster_gam <- raster(file.path(input_sdm_gam, sdm_name_gam))
        sdm_raster <- sdm_raster_glm + sdm_raster_gam
        sdm_raster <- sdm_raster==2
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
}




####################################################
###################   STEP 2 #######################
#          2_4 - sort sst and sbt rasters          #
####################################################
####################################################

bathymetry <- raster('../raw_fish_ranges_rasters/depth/gb_depth_caribbean_pacific.asc')
species_traits_caribbean <- readRDS('../raw_fish_ranges_occ_caribbean/raw/fish_traits_camille_2019/Traits_Caraibean_Fish.rds')
species_traits_pacific <- readRDS('../raw_fish_ranges_occ_caribbean/raw/fish_traits_camille_2019/Traits_pacific_fish.rds')

#see if species that are in both databases have the same traits
sp_car_pac <- species_traits_caribbean[intersect(rownames(species_traits_caribbean), rownames(species_traits_pacific)),]
sp_pac_car <- species_traits_pacific[intersect(rownames(species_traits_caribbean), rownames(species_traits_pacific)),]
sp_car_pac_comb <- as.data.frame(cbind(rownames(sp_car_pac), rownames(sp_pac_car), sp_car_pac[,1:2], sp_pac_car[, 1:2]))

species_traits_pacific <- species_traits_pacific[!is.element(rownames(species_traits_pacific), rownames(sp_car_pac_comb)),]
species_traits_combined <- rbind(species_traits_caribbean, species_traits_pacific)


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
      
      #foreach(j = species_names, .packages=c('raster', 'tools')) %dopar% {
      for(j in species_names){  
        sp_name_camille <- gsub(' ', '_', j)
        if(is.element(sp_name_camille, rownames(species_traits_combined))){ #take data from trait table
          
          min_depth <- as.numeric(species_traits_combined[rownames(species_traits_combined)==sp_name_camille, 1][[1]])
          max_depth <- as.numeric(species_traits_combined[rownames(species_traits_combined)==sp_name_camille, 2][[1]])
          cat('min =', min_depth, 'max = ', max_depth, 'for species', j, '\n')
          sp_depth_trait <- mean(c(min_depth, max_depth))
          
        } else{ #take ranges from observation data
          
          cat('No traits available for species', j, '--> using occurrence data \n')
          sp_occ <- read.table(file.path('../raw_fish_ranges_occ_caribbean/point_selected', paste0(j, '.txt')))
          
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



####################################################
###################   STEP 3 #######################
#            2_5 - combine rasters              #
####################################################
####################################################


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




####################################################
###################   STEP 4 #######################
#           2_6 - sum richness and plot            #
####################################################
####################################################

coast <- readOGR("../raw_fish_ranges_shapefiles/coastline/GSHHS_l_L1_disolve.shp")
source("mapping_functions.r")

k <- method
l <- levels

richness_sst <- raster(file.path('../fish_ranges_output_combined_sst', k, l, 'richness_rasters_sst', 'species_richness_combined_sst.tif'))
richness_sbt <- raster(file.path('../fish_ranges_output_combined_sbt', k, l, 'richness_rasters_sbt', 'species_richness_combined_sbt.tif'))
richness_total <- sum(richness_sbt, richness_sst, na.rm=T)
writeRaster(richness_total, file.path('../fish_ranges_output_combined',  paste0('species_richness_combined_', k, '_', l, '.tif')))


#################          PLOT          #################
richness_total <- raster(file.path('../fish_ranges_output_combined',  paste0('species_richness_combined_', k, '_', l, '.tif')))

richness_total[richness_total<10] <- NA
species_richness <- as.data.frame(richness_total, xy=T)
names(species_richness)[3] <- 'layer'

#### Mapping
png(filename=file.path('../fish_ranges_output_combined',  paste0('species_richness_combined_', k, '_', l, '.png')), width=20,height=18,units="cm",res=600)
#layout(matrix(c(1,2,3,4), ncol=2), heights=c(7,2))
layout(matrix(c(1,1,1,2), ncol=1, byrow=F))
layout.show(2)
#detach(package:Hmisc, unload=T)
buff <- 3
x_lim <- c(-117.9937+buff, -53.17071-buff)
y_lim <- c(-7.978883+buff, 37.00712-buff)


### Fig a
par(mar=c(1.5,3.5,3.2,3), xpd=F)
col <- rev(c("#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4"))
brks <- seq(0, 800, by=1)
# brks_legend <- c(0, 100, 200, 275, 350, 475)
# #brks_legend <- c(0, 100, 200, 460, 470, 475)
# col.tot <- colorRampPalette(col)
# col_tot <- col.tot(length(brks))
# #col_values <- c(50, 150, 240.....)#automate this or make gradient
# col_legend <- c(col_tot[is.element(brks, brks_legend[-length(brks_legend)])], col[length(col)])

#brks_legend <- c(0, 475, by=100)
#brks <- seq(0, round_any(max(species_richness$layer, na.rm=T), 50, f=ceiling), length.out=6)
Map_the_world_2(Data_PA=species_richness$layer,coord_X=species_richness$x, coord_Y=species_richness$y,
                col=col,breaks=brks,include.lowest=TRUE,xlim=x_lim,
                ylim=y_lim,legend=FALSE,coastline=TRUE,names_fig="",cex_point=0.39,
                cex_names_fig=1.5,x_names=1.8,y_names=39,font_txt=2)
par(xpd=NA)
text(-85, -12, labels=paste('Species richness', k, l), cex=1.2)


#legend
par(mar=c(2,3,1,3))
#legend_image <- as.raster(matrix(colorRampPalette(col)(length(brks)-1), ncol=1))
#legend_df <- as.data.frame(matrix(cbind(seq(0, 475, length.out = (length(brks)-1)), rep(0.3, (length(brks)-1)), colorRampPalette(col)(length(brks)-1)), ncol=3))
#legend_df[,1] <- as.numeric(as.character(legend_df[,1]))
#legend_df[,2] <- as.numeric(as.character(legend_df[,2]))
#legend_image <- rasterFromXYZ(legend_df)

legend_df <- as.data.frame(cbind(seq(0, length(brks)-2, length.out=(length(brks)-1)), rep(0.25, (length(brks)-1)), seq(0, 466, length.out = (length(brks)-1))))
legend_image <- rasterFromXYZ(legend_df, res=0.25)
#plot(legend_image, col=colorRampPalette(col)(length(brks)-1), ext=c(0,475, 0.249, 0.251))
#plot(legend_image, col=NA, box=F, axes=F, legend=F)
plot.new()
plot(legend_image, legend.only=T, col=colorRampPalette(col)(length(brks)-1), horizontal=T, smallplot=c(0.2, 0.8, 0.4, 0.5), 
     axis.args=list(labels = seq(0, 475, length.out = 6)))

#smallplot=c(-110, -100, -115, -100)
#   #  axis.args=list(at=seq(0, 475, length.out = 6), labels=seq(0, 475, length.out = 6)), cex.axis=0.6)
#
#plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
#text(x = seq(0,1,length.out=6), y=0.5, labels = seq(0, 475, length.out = 6))

#rasterImage(legend_image, xleft=0.05, xright=10, ybottom=0, ytop=0.3, angle=-90)

legend("bottom", legend=c("Land","Pixels with less","than 10 species"), col=c("gray25","gray80","white"),
       pch=15, pt.cex=1.5, cex=1.2, xjust=0.5, bty="n", ncol=3)


### shared legend
#par(mar=c(2,3,1,1))
#plot(rnorm(100),type="n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
#legend("center",legend=c(levels(cut(species_richness$layer, breaks=brks_legend,include.lowest=TRUE,dig.lab=4)),
#                         "Land","Pixels with less","than 10 species"),col=c(col_legend,"gray25","gray80","white"),
#       pch=15,pt.cex=1.5,cex=1.2,xjust=0.5,text.col="black",bty="n",ncol=4) #title="Predicted species richness"

dev.off()

#brks_legend[-length(brks_legend)]