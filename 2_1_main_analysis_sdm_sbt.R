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
#setwd("~/Documents/ETH/Landscape Ecology/Paper - Fish Ranges Caribbean 2020/data_2020_04_14/raw_fish_ranges_scripts")
output_dir <- file.path('../fish_ranges_output_sdm_sbt')

cmon.files=list.files("common/",full.names = T)
sapply(cmon.files,source)

lib_vect <- c("raster","rgdal","maptools", "tools","dismo","cluster","class",
              "gam","gbm","randomForest","ROCR","parallel")

sapply(lib_vect,require,character.only=TRUE)

### =========================================================================
### Definitions depending on cluster/no cluster
### =========================================================================

if(interactive()){
  ######   IF NOT CLUSTER   #########
  input_dir <- '../raw_fish_ranges_occ_caribbean/point_selected_filtered_sel'
  fromhere <- 1
  tohere <- 250
  tempfiles <- 'temp01'
  ###################################
  


} else {
  #########       CLUSTER   #########
  args <- commandArgs(trailingOnly = TRUE)
  fromhere <- as.numeric(args[1])
  tohere <- as.numeric(args[2])
  tempfiles <- as.character(args[3])
  input_dir <- as.character(args[4])
  ###################################

  
}

env.stk= stack('../raw_fish_ranges_rasters/sbt_mean/bedtemp_caribbean_pacific.asc',
               '../raw_fish_ranges_rasters/depth/gb_depth_caribbean_pacific.asc')

# Take the log of NPP

pred_sdm <- c('Seabed_Temp', 'Bathymetry') 
names(env.stk)=pred_sdm
env.stk <-  mask(env.stk, calc(env.stk,fun = sum))

# Load preped psuedoabsences
load("../raw_fish_ranges_occ_caribbean/pseudoabsence_data/pseu.abs.sbt.RData")

### =========================================================================
### Definitions
### =========================================================================

proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

### =========================================================================
### Loop over species and create maps
### =========================================================================

sps <- list.files(input_dir) 
for (id in fromhere:tohere){
  

  ### =========================================================================
  ### Prepare data for species id
  ### =========================================================================
  
  # Get name
  env.stk_i=env.stk
  spi <- sps[id]
  spi_name <- file_path_sans_ext(spi)
  
  print(paste0("working with        ", spi_name, "         [index:", id, "]"))
  
  # Load observations
  obs=read.table(file.path(input_dir, spi), header=T)
  spp=SpatialPoints(obs,proj4string = proj)
  
  # Check whether presences are associated with any NA
  xtr=extract(env.stk_i,spp)
  sna=apply(xtr,1,function(x){any(is.na(x))})
  spp=spp[which(!sna),]
  
  # Thin presences
  if(nrow(spp@coords)<3000 & nrow(spp@coords)>=30){
    spp=thin_them(spp,5)
  } else if (nrow(spp@coords)>=3000){
    spp=upsample_strategic(spp,5,nrow(spp@coords),warnig=FALSE)
  }

  # Create pseudoabsences
  pseu.abs_i=copy_PseuAbs(env.stack=env.stk_i,
                          pres=spp,
                          taxon=spi_name,
                          pseu.abs)
  
  # plot(pseu.abs_i)
  
  ### =========================================================================
  ### Decide on predictor number depending on available presence observations
  ### =========================================================================
  # Note that for the analysis I assume that the presences are already cleaned,
  # i.e., points falling in the ocean are removed, and very strong density 
  # biases are removed, e.g., with the rasterize function. If this is not done,
  # there may be problems. So better to run a preparation step beforehand.
  
  ssize=floor(length(which(pseu.abs_i@pa==1))/10)
  
  # If there are not at least 20 presence observations, do nothing
  if(ssize<1){
    next
    # If there are less than 90 observations reduce predictor set so that
    # at least 10 observations are available per predictor
  }
  
  ### =========================================================================
  ### Prepare data for modelling
  ### =========================================================================
  
  # Make PCA for block crossvalidation
  # prdi=prcomp(pseu.abs@env_vars,center=T,scale.=T)
  # bldt=prdi$x[,1:2]
  
  # library("NMOF")
  
  # define strata for blocks
  #falls weniger als 10 obs --> untere Zeile
  #blks=make_blocks(df=bldt,nstrat=3,nclusters=12,pres=pseu.abs@pa)
  wip=which(pseu.abs_i@pa==1)
  blkp=sample(rep(1:3,each=ceiling(length(wip)/3)),ceiling(length(wip)/3)*3)[1:length(wip)]
  wia=which(pseu.abs_i@pa==0)
  blka=sample(1:3,size = length(wia),replace = T)
  
  blks=rep(NA,length(pseu.abs_i@pa))
  blks[wip]=blkp
  blks[wia]=blka

  # detach("package:NMOF")
  
  ### =========================================================================
  ### Define SDMs
  ### =========================================================================
  
  ### Define formulas
  # GLMs
  form.glm.s=as.formula(paste("Presence~",paste(names(env.stk_i),collapse="+"))) # simple
  form.glm.i=as.formula(paste("Presence~",paste(paste0("poly(",names(env.stk_i),",2)"),collapse="+"))) # intermediate
  
  # GLM complex
  if(nlayers(env.stk_i)>1){
    cmbs<-combn(names(env.stk_i),2)
    pst<-apply(cmbs,2,paste,collapse=":")
    int.part=paste(pst,collapse="+")
    form.glm.c=as.formula(paste(paste("Presence~",paste(paste0("poly(",names(env.stk_i),",4)"),collapse="+")),int.part,sep="+")) # complex
  } else {
    form.glm.c=as.formula(paste(paste("Presence~",paste(paste0("poly(",names(env.stk_i),",4)"),collapse="+")))) # complex
  }

  # GAMs
  form.gam.s=as.formula(paste("Presence~",paste(paste0("s(",names(env.stk_i),",df=1.5)"),collapse="+"))) # simple
  form.gam.i=as.formula(paste("Presence~",paste(paste0("s(",names(env.stk_i),",df=3)"),collapse="+"))) # intermediate
  form.gam.c=as.formula(paste("Presence~",paste(paste0("s(",names(env.stk_i),",df=10)"),collapse="+"))) # complex
  
  # Trees
  form.trees=Presence ~ .
  
  ### define model settings
  modinp=list(multi("glm",list(formula=form.glm.s,family="binomial"),"glm-simple",step=FALSE,weight=FALSE),
              multi("glm",list(formula=form.glm.i,family="binomial"),"glm-interm",step=FALSE,weight=FALSE),
              multi("glm",list(formula=form.glm.c,family="binomial"),"glm-complex",step=FALSE,weight=FALSE),
              multi("gam",list(formula=form.gam.s,family="binomial"),"gam-simple",step=FALSE,weight=FALSE),
              multi("gam",list(formula=form.gam.i,family="binomial"),"gam-interm",step=FALSE,weight=FALSE),
              multi("gam",list(formula=form.gam.c,family="binomial"),"gam-complex",step=FALSE,weight=FALSE),
              multi("gbm",list(formula=form.trees,
                               distribution = "bernoulli",
                               interaction.depth = 5,
                               shrinkage=.005,
                               n.trees = 100),"gbm-simple",weight=FALSE),
              multi("gbm",list(formula=form.trees,
                               distribution = "bernoulli",
                               interaction.depth = 5,
                               shrinkage=.005,
                               n.trees = 300),"gbm-interm",weight=FALSE),
              multi("gbm",list(formula=form.trees,
                               distribution = "bernoulli",
                               interaction.depth = 5,
                               shrinkage=.005,
                               n.trees = 10000),"gbm-complex",weight=FALSE),
              multi("randomForest",list(formula=form.trees,ntree=500,nodesize=40),"rf-simple"),
              multi("randomForest",list(formula=form.trees,ntree=500,nodesize=20),"rf-interm"),
              multi("randomForest",list(formula=form.trees,ntree=500,nodesize=1),"rf-complex"))
  
  ### =========================================================================
  ### Run models
  ### =========================================================================
  
  modis=wsl.flex(x=pseu.abs_i,
                 replicatetype="block-cv",
                 reps=3,
                 strata=blks,
                 project="fish_diversity",
                 mod_args=modinp)
  
  ### =========================================================================
  ### evaluate
  ### =========================================================================


  
  #smev=summary(evals)
  #ord=sort(smev["tss",],decreasing=T)
  #top6=which(colnames(smev)%in%names(ord)[1:6])
  
  #modinp_top=modinp[top6]
  modinp_top=modinp
  
  ### =========================================================================
  ### fit models for prediction
  ### =========================================================================
  
  prmod=wsl.flex(x=pseu.abs_i,
                 replicatetype="none",
                 reps=1,
                 project="tree_map",
                 mod_args=modinp_top)
  
  ### =========================================================================
  ### Predict
  ### =========================================================================
  

  
  # Null prediction to assess change
  prenull=wsl.predict(prmod,
                      predat=env.stk_i,
                      #thres=thrs,
                      #thres=thrs[top6],
                      write=TRUE,
                      output_dir=file.path(output_dir, 'prob'))
  
  # Null prediction to assess change
  ### Get thresholds
  evals<-wsl.evaluate(modis,crit="maxTSS",pseuAbsCor=T)
  thrs=get_thres(evals)
 prenull=wsl.predict(prmod,
                     predat=env.stk_i,
                     thres=thrs,
                     #thres=thrs[top6],
                     write=TRUE,
                     output_dir=file.path(output_dir, 'tss'))
 
 # Null prediction to assess change
 ### Get thresholds
 evals<-wsl.evaluate(modis,crit="pp=op",pseuAbsCor=T)
 thrs=get_thres(evals)
 prenull=wsl.predict(prmod,
                     predat=env.stk_i,
                     thres=thrs,
                     #thres=thrs[top6],
                     write=TRUE,
                     output_dir=file.path(output_dir, 'op'))

  ### =========================================================================
  ### combine
  ### =========================================================================

  # Load saved binary predictions
  spina=gsub(" ","_",spi_name)
  #if(!dir.exists(paste0(output_dir,"/",spina))){
  #  dir.create(paste0(output_dir,"/",spina), recursive = T)
  #}
  flbin=list.files(paste0(output_dir,"/",spina),full.names = T,pattern="BinaryPredictions")
  stk=stack(flbin)
  
  # Combine and save
  flcom=paste0(output_dir,"/ComitteeVote_",spina,".tif")
  calc(stk,sum,filename=flcom,datatype="INT1U",overwrite=TRUE)
  

  print(paste0("COMPLETED: [ ",round(id/length(sps),2)*100, "% ]"))
}

unlink(file.path(paste('output', Sys.Date(), sep='_'), tempfiles), recursive=T)




