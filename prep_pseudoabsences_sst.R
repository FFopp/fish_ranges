##################################################
## Description: points to range raster 
## 
## Date: 2018-03-26 13:57:16
## Authors: Fabian Fopp (fabian.fopp@usys.ethz.ch)
##          Philipp Brun (philipp.brun@wsl.ch)
##################################################

### =========================================================================
### Set directories and source functions
### =========================================================================

#setwd("/Users/Fabian/Documents/ETH/Landscape Ecology/Treemap/treemap_occurrence_mapping/treemap_occurrence_mapping-Occurrence_mapping_pb/scripts")
#setwd("C:/Users/brunp/Desktop/fish_ranges_philipp Kopie/raw_fish_ranges_scripts/")

output_dir <- file.path('../fish_ranges_output_sdm_sst')

cmon.files=list.files("common/",full.names = T)
sapply(cmon.files,source)

lib_vect <- c("raster","rgdal","maptools", "tools","dismo","cluster","class",
              "gam","gbm","randomForest","ROCR","parallel")

sapply(lib_vect,require,character.only=TRUE)
env.stk= stack('../raw_fish_ranges_rasters/depth/gb_depth_caribbean_pacific.asc',
               '../raw_fish_ranges_rasters/sstmean/sstmean_caribbean_pacific.asc')

# Take the log of NPP

pred_sdm <- c('Bathymetry', 'Seasurface_Temp') 
names(env.stk)=pred_sdm
env.stk <-  mask(env.stk, calc(env.stk,fun = sum))


### =========================================================================
### Run the function
### =========================================================================

load('../raw_fish_ranges_occ_caribbean/pseudoabsence_data/all_obs.RData')  #all_obs_sp

# Create pseudoabsences
pseu.abs=wsl.samplePseuAbs(type="geographic",
                           n=8000,
                           env.stack=env.stk,
                           add.strat=0.2,
                           template_dir="../raw_fish_ranges_occ_caribbean/pseudoabsence_data/",
                           target.group_dir="../raw_fish_ranges_occ_caribbean/point_selected_filtered/",
                           env.strat_path="../raw_fish_ranges_occ_caribbean/pseudoabsence_data/",
                           pres=all_obs_sp,
                           geores_fact=1,
                           force_spat_thin="absences",
                           limdist=5)
### =========================================================================
### Save object
### =========================================================================

save(pseu.abs,file="../raw_fish_ranges_occ_caribbean/pseudoabsence_data/pseu.abs.sst.RData")
