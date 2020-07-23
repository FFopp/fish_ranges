
library(raster)
lib_vect <- c('raster', 'plyr', 'rgdal')
sapply(lib_vect, library, character.only=T)

coast <- readOGR("../raw_fish_ranges_shapefiles/coastline/GSHHS_l_L1_disolve.shp")
source("mapping_functions.r")

method <- c('glm', 'gam', 'gbm', 'rf')
levels <- c('complex', 'interm', 'simple')

for(k in method){
  for(l in levels){

    richness_sst <- raster(file.path('../fish_ranges_output_combined_sst', k, l, 'richness_rasters_sst', 'species_richness_combined_sst.tif'))
    richness_sbt <- raster(file.path('../fish_ranges_output_combined_sbt', k, l, 'richness_rasters_sbt', 'species_richness_combined_sbt.tif'))
    richness_total <- sum(richness_sbt, richness_sst, na.rm=T)
    writeRaster(richness_total, file.path('../fish_ranges_output_combined',  paste0('species_richness_combined_', k, '_', l, '.tif')))
    
    #richness_sst <- raster('fish_ranges_output_combined_sst/richness_rasters_sst/species_richness_combined_sst.tif')
    #richness_sbt <- raster('fish_ranges_output_combined_sbt/richness_rasters_sbt/species_richness_combined_sbt.tif')
    #richness_total <- sum(richness_sbt, richness_sst, na.rm=T)
    #writeRaster(richness_total, 'fish_ranges_output_combined/species_richness_combined.tif')
    

    
    #################          PLOT          #################
    #richness_total <- raster(file.path('../fish_ranges_output_combined',  paste0('species_richness_combined_', k, '_', l, '.tif')))

    
    richness_total[richness_total<10] <- NA
    species_richness <- as.data.frame(richness_total, xy=T)
    names(species_richness)[3] <- 'layer'
    
    #### Mapping
    png(filename=file.path('../fish_ranges_output_combined',  paste0('species_richness_combined_', k, '_', l, '.png')), width=20,height=18,units="cm",res=600)
    #layout(matrix(c(1,2,3,4), ncol=2), heights=c(7,2))
    layout(matrix(c(1,1,1,2), ncol=1))
    layout.show(2)
    #detach(package:Hmisc, unload=T)
    buff <- 3
    x_lim <- c(-117.9937+buff, -53.17071-buff)
    y_lim <- c(-7.978883+buff, 37.00712-buff)
    
    
    ### Fig a
    par(mar=c(1.5,3.5,3.2,1), xpd=F)
    col <- rev(c("#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4"))
    #brks<- c(0,100,200,350,450,550)
    brks <- seq(0, round_any(max(species_richness$layer, na.rm=T), 50, f=ceiling), length.out=6)
    Map_the_world_2(Data_PA=species_richness$layer,coord_X=species_richness$x, coord_Y=species_richness$y,
                    col=col,breaks=brks,include.lowest=TRUE,xlim=x_lim,
                    ylim=y_lim,legend=FALSE,coastline=TRUE,names_fig="",cex_point=0.39,
                    cex_names_fig=1.5,x_names=1.8,y_names=39,font_txt=2)
    par(xpd=NA)
    text(-85, -12, labels=paste('Species richness', k, l), cex=1.2)
    
    
    ### shared legend
    par(mar=c(2,3,1,1))
    plot(rnorm(100),type="n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
    legend("center",legend=c(levels(cut(species_richness$layer, breaks=brks,include.lowest=TRUE,dig.lab=4)),
                             "Land","Pixels with less","than 10 species"),col=c(colorRampPalette(col)(length(brks)-1),"gray25","gray80","white"),
           pch=15,pt.cex=1.5,cex=1.2,xjust=0.5,text.col="black",bty="n",ncol=4) #title="Predicted species richness"
    
    dev.off()
  }
}
