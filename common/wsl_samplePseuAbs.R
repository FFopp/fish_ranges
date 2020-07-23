### =========================================================================
### define wsl.samplePseuAbs
### =========================================================================
#' Sample pseudoabsences using various strategies
#'
#' Flexlible function to sample pseudoabsences with various strategies and
#' store the results in a 'wsl.pseudoabsences' object  that can be
#' passed on in the 'wsl.biodiv' pipeline.
#'
#' @param n number of pseudoabsence points desired. Default is 10000.
#' @param env.stack RasterStack with environmental layers for sampling and
#' extraction
#' @param type Desired sampling strategy. Options are 'geographic', 'density',
#' 'random', 'target.group', 'geo.strat', 'env.strat' and 'env.semi.strat'
#' (see details). Default is 'geographic'
#' @param add.strat Fraction between 0 and 1; should strategy be complemented
#' by a fraction of environmental strata.
#' @param pres SpatialPoints object. Location of presence points. Necessary for
#' 'geographic' sampling strategy and the best adding point if the downstream
#' functions will be used.
#' @param taxon Character; name of taxon of interest to keep track of meta
#' information.
#' @param geodist_fact Factor to adjust spatial autocorrelation lengths: for
#' 'geographic' pseudoabsence point patterns, values below 1 increase
#' autocorrelation length; values above 1 decrease it; for 'density' sampling
#' it is the other way around.
#' @param geores_fact Aggregation factor 'geographic' template. Larger values
#' save calculation time, but derease resolution of sampling points.
#' @param template_dir Directory where template raster for 'geographic'
#' sampling should be saved in/loaded from. If NA, nothing will be saved; if
#' provided, template will be saved in/loaded from directory depending on
#' whether a file already exists.
#' @param geo_nrep number of replicates of geographic models. More will create
#' a smoother pattern but increase computation time.
#' @param target.group_dir Directory where xy files of traget group taxa are
#' stored. Must be supplied if sampling stragegy is 'target.group', must contain
#' a column names 'x' and 'y' with coordinates in the same projection as other
#' spatial data
#' @param env.strat_path Directory where sample of environmental strata for
#' 'env.strat' sampling should be saved in/loaded from. If NA, nothing will be
#' saved; if provided, environmental strata will be saved in/loaded from
#' directory depending on whether a file already exists.
#' @param rAll should all data be read into memory for computation of environmental
#' strata? this is faster but you may run into memory issues for large rasters.
#' @param force_spat_thin Should minimum distance be enforced between points? Options
#' are 'no', 'presences', 'absences', and 'both'. By default thinning is defiend for
#' pseudoabsences from 'geographic', 'density', 'random', and 'geo.strat' methods
#' with minimum distance according to the resoultion of the template raster.
#' 'presences' takes the minimum distance criterion from the template raster over
#' to the 'presence' points; 'absences' takes the criterion over to 'env.strat',
#' 'env.semi.strat', and 'target.group'; 'both' does it for both.
#' @param limdist The minimum distance accepted for spatial thinning. Units
#' should be km if the spatial data is projected, otherwise the units of the
#' coordinate reference system used. If no value is supplied, the maximum
#' distance between two cells of the template raster will be taken.
#' @details 'geographic' samples pseudoabsences with a sampling probability.
#' inversely proportional to the geographic distance to presence observations.
#' density' samples pseudoabsences proportional to the density of presence
#' observations. 'random' samples pseudoabsences randomly with a sampling
#' probability proportional the area of the cells. 'target.group' samples
#' pseudoabsences from the presences of the taxa of the target group, attempting
#' to correct for sampling bias. It depends on a directory with taxa defined by the
#' user as target group. 'geo.strat' samples pseudoabsences geographically
#' stratified either on a plane, or on a sphere depending on the projection of
#' the supplied env.stack. 'env.strat' samples pseudoabsences environmentally
#' stratified. Points are sampled form all realized combinations of environmental
#' conditions occuring in the environmental stack that have a minimal occurrence
#' frequency. By default environmental strata are calculated based on all
#' raster layers supplied. If a directory is supplied as 'env.strat_path', a
#' large sample of stratified points will be saved to speed up computations for
#' follow-up species. If environmental strata based on different predictors than
#' supplied are preferred 'env.strat_path'can be an .RData file from a
#' previous sampling of strata form different environmental predictors. 'env.semi.strat'
#' is similar to 'env.strat' but samples environmental strata proportional to the
#' logarithm of the area they cover.
#'
#' @return an object of class 'wsl.pseudoabsences' that can be plotted and
#' passed on to wsl.flex
#' @author Philipp
#' @export
#' @examples
#'### =========================================================================
#'### Data preparation
#'### =========================================================================
#'
#'# Predictors
#'bio=getData('worldclim',var='bio',lon=16, lat=48,res=.5)
#'bio=bio[[c(1,4,12)]]
#'
#'# install.packages("rgbif")
#'library(rgbif)
#'# extract species
#'spn='Boletus aestivalis'
#'xt=as.vector(extent(bio))
#'baest <- occ_search(scientificName=spn,
#'                   hasCoordinate=TRUE,
#'                   decimalLongitude=paste0(xt[1],",",xt[3]),
#'                   decimalLatitude=paste0(xt[2],",",xt[4]))
#'
#'pbaest=baest$data[,c('decimalLongitude','decimalLatitude')]
#'baest_spp=SpatialPoints(pbaest,proj4string = crs(bio))
#'
#'# extract target group
#'targr <- occ_search(familyKey = 8789,
#'                    hasCoordinate=TRUE,
#'                    limit = 10000,
#'                    decimalLongitude=paste0(xt[1],",",xt[3]),
#'                    decimalLatitude=paste0(xt[2],",",xt[4]))
#'
#'ptargr=as.matrix(targr$data[,c('decimalLongitude','decimalLatitude')])
#'colnames(ptargr)=c("x","y")
#'
#'# create temporary directory for target.group info
#'tdir=paste0(tempdir(),"/trgr")
#'dir.create(tdir)
#'write.table(ptargr,file=paste0(tdir,"/targetxy.txt"),row.names = F)
#'
#'# create temporary directory for template raster and env strata
#'strdir=paste0(tempdir(),"/str")
#'dir.create(strdir)
#'
#'# Note that for these should not be temporary files for a real analysis.
#'
#'### =========================================================================
#'### Sample pseudoabsences
#'### =========================================================================
#'
#'# Geograhpic method with 20% env strata
#'pseu.abs1=wsl.samplePseuAbs(type="geographic",
#'                           n=5000,
#'                           env.stack=bio,
#'                           pres=baest_spp,
#'                           add.strat=0.2,
#'                           template_dir=strdir,
#'                           env.strat_path=strdir,
#'                           geodist_fact=1,
#'                           geores_fact=3,
#'                           geo_nrep=7,
#'                           taxon=spn)
#'
#'plot(pseu.abs1)
#'
#'# Only geographic with longer autocorrelation length
#'pseu.abs2=wsl.samplePseuAbs(type="geographic",
#'                            n=5000,
#'                            env.stack=bio,
#'                            pres=baest_spp,
#'                            add.strat=0,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir,
#'                            geodist_fact=.5,
#'                            geores_fact=3,
#'                            geo_nrep=7,
#'                            taxon=spn)
#'
#'plot(pseu.abs2)
#'
#'# Random and thin presences
#'pseu.abs3=wsl.samplePseuAbs(type="random",
#'                            n=5000,
#'                            env.stack=bio,
#'                            template_dir=strdir,
#'                            pres=baest_spp,
#'                            geores_fact=3,
#'                            add.strat=0,
#'                            taxon=spn,
#'                            force_spat_thin="presences")
#'
#'plot(pseu.abs3)
#'
#'# Geo.start
#'pseu.abs4=wsl.samplePseuAbs(type="geo.strat",
#'                            n=5000,
#'                            env.stack=bio,
#'                            template_dir=strdir,
#'                            pres=baest_spp,
#'                            geores_fact=3,
#'                            add.strat=0,
#'                            taxon=spn)
#'
#'plot(pseu.abs4)
#'
#'# Target group with 20% env strat
#'pseu.abs5=wsl.samplePseuAbs(type="target.group",
#'                          n=5000,
#'                          env.stack=bio,
#'                          template_dir=strdir,
#'                          target.group_dir=tdir,
#'                          env.strat_path=strdir,
#'                          geores_fact=3,
#'                          pres=baest_spp,
#'                          add.strat=0.2,
#'                          taxon=spn,
#'                          force_spat_thin="both")
#'
#'plot(pseu.abs5)
#'
#'# Environmental semi-stratified
#'pseu_abs6=wsl.samplePseuAbs(n = 5000,
#'                            env.stack=bio,
#'                            type = "env.semi.strat",
#'                            add.strat = 0,
#'                            pres = baest_spp,
#'                            taxon = spn,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir)
#'
#'plot(pseu_abs6)
#'
#'# Environmental semi-stratified with min dist
#'pseu_abs7=wsl.samplePseuAbs(n = 5000,
#'                            env.stack=bio,
#'                            type = "env.semi.strat",
#'                            add.strat = 0,
#'                            geores_fact=3,
#'                            pres = baest_spp,
#'                            taxon = spn,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir)
#'
#'plot(pseu_abs7)
#'
#'# Density dependent
#'pseu_abs8=wsl.samplePseuAbs(n = 5000,
#'                            env.stack=bio,
#'                            type = "density",
#'                            add.strat = 0,
#'                            pres = baest_spp,
#'                            taxon = spn,
#'                            geores_fact=3,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir)
#'
#'plot(pseu_abs8)
#'
#'### =========================================================================
#'### Fit SDMs
#'### =========================================================================
#'
#' # Define model settings
#'vrs=names(bio)
#'form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
#'form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
#'form.tree=as.formula(Presence ~ .)
#'
#'modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=FALSE),
#'            multi("gbm",list(formula=form.tree,
#'                 distribution = "bernoulli",
#'                 interaction.depth = 1,
#'                 shrinkage=.01,
#'                 n.trees = 3500),"gbm-simple"),
#'            multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE),
#'            multi("randomForest",list(formula=form.tree,ntree=500,maxnodes=NULL),"waud1"))
#'
#'# Fit models using wsl.flex function
#'library(gam)
#'library(gbm)
#'library(randomForest)
#'modi5=wsl.flex(x=pseu.abs5,
#'               replicatetype="block-cv",
#'               reps=3,
#'               strata=sample(1:3,nrow(pseu.abs5@env_vars),replace=TRUE),
#'               project="multitest",
#'               mod_args=modinp)
#'
#'### =========================================================================
#'### Evaluate
#'### =========================================================================
#'
#'# Evaluate the models with wsl.evaluate function
#'eval5<-wsl.evaluate(modi5,crit="maxTSS",pseuAbsCor = TRUE)
#'summary(eval5)
#'
#'### =========================================================================
#'### Predict
#'### =========================================================================
#'
#'modi_pred=wsl.flex(x=pseu.abs5,
#'                   replicatetype="none",
#'                   reps=1,
#'                   project="test_pred",
#'                   mod_args=modinp[c(1:2)])
#'
#'# Get thresholds
#'thr.5=get_thres(eval5)[1:2]
#'
#'### Make some predictions
#'pred5=wsl.predict(modi_pred,predat=bio,thres = thr.5)
#'
#'par(mfrow=c(1,2))
#'plot(pred5@predictions[[1]]$`glm-simple`,main="GLM")
#'plot(pred5@predictions[[1]]$`gbm-simple`,main="GBM")
#'
wsl.samplePseuAbs<-function(n=10000,
                            env.stack,
                            type="geographic",
                            add.strat=0,
                            pres=numeric(),
                            taxon=character(),
                            geodist_fact=1,
                            geores_fact=20,
                            template_dir=NA,
                            geo_nrep=7,
                            target.group_dir=NA,
                            env.strat_path=NA,
                            rAll=TRUE,
                            force_spat_thin="no",
                            limdist=NA){

  ### ------------------------
  ### Check input and prepare
  ### ------------------------

  if(add.strat<0 | add.strat>1){
    stop("add.strat represents the fraction of pseudoabsences that are
         sampled environmentally stratified and should be between 0 and 1!")
  }

  possibtype=c("geographic","random","target.group","geo.strat","density",
               "env.strat","env.semi.strat")
  if(length(type)!=1 | !(type%in%possibtype)){
    stop("Invalid specification of pseudoabsence sampling type!")
  }

  possibthin=c("no","presences","absences","both")
  if(length(type)!=1 | !(force_spat_thin%in%possibthin)){
    stop("Invalid specification of spatial thinning method!")
  }

  if(grepl("env",type)){
    add.strat=1
  }

  ### ------------------------
  ### generate wsl.pseudoabsences object and add meta info
  ### ------------------------

  out<-preva.meta(type="pseudoabsence")

  tpnam=type
  if(tpnam=="geographic"){
    tpnam=paste0(tpnam,"_w",geodist_fact)
  }
  if(add.strat>0 && !grepl("env",type)){
    tpnam=paste0(tpnam,"_x_",add.strat,"env_strata")
  }

  out@meta$type=tpnam
  out@meta$taxon=taxon
  out@meta$force_spat_thin=force_spat_thin
  out@meta$template_file=NA

  call=match.call()
  out@call<-call

  ### ------------------------
  ### Prepare template raster
  ### ------------------------

  if(type%in%c("geographic","random","geo.strat","density") | force_spat_thin%in%c("presences","both")){

    # no template directory is supplied, just
    # calculate from scratch
    if(is.na(template_dir)){

      rst=aggregate(env.stack[[1]],
                    fact=geores_fact,
                    fun=function(x,na.rm){
                      ifelse(all(is.na(x)),NA,1)
                    },na.rm=T)
      # if template directory is supplied load if file exists
      # otherwise writeRaster
    } else {
      ptrn=paste0("template",geores_fact,".tif")
      tmfl=list.files(template_dir,pattern=ptrn,full.names=TRUE)

      if(length(tmfl)>0){
        rst=raster(tmfl)
      } else {
        tmfl=paste0(template_dir,"/",ptrn)
        rst=aggregate(env.stack[[1]],
                      fact=geores_fact,
                      fun=function(x,na.rm){
                        ifelse(all(is.na(x)),NA,1)
                      },na.rm=T,filename=tmfl,overwrite=TRUE)
      }
    }
    crs(rst)=crs(env.stack)

    # Calculate (latitudinal) distance between cells for potential
    # downstream analyses if no minim is provided
    if(is.na(limdist)){
      proje=grepl("longlat",crs(rst))
      dpp=SpatialPoints(coordinates(rst)[c(1,1+dim(rst)[2]),],
                        proj4string = crs(rst))
      limdist=spDists(dpp,longlat = proje)[1,2]
    }

    
    # Write path to template file into meta information
    out@meta$template_file=paste0(template_dir,"/template",geores_fact,".tif")
  }

  ### ------------------------
  ### Extract and refine presences
  ### ------------------------

  if(length(pres)>0){

    xt_pres=extract(env.stack,pres)
    sna=apply(xt_pres,1,function(x){
      any(is.na(x))
    })

    if(length(which(sna))>0){
      pres=pres[-which(sna)]
      xt_pres=xt_pres[-which(sna),]
      cat(paste0(length(which(sna))," non-matching presences removed..\n"))
    }

    if(force_spat_thin%in%c("presences","both")){

      if(nrow(pres@coords)<3000){
        tpp=thin_them(pres,limdist)
      } else {
        tpp=upsample_strategic(pres,limdist,nrow(pres@coords),warnig=FALSE)
      }



      cat(paste0(length(pres)-length(tpp)," presences removed to obtain min distance of ",
                 round(limdist,digits=2),"..\n"))
      pres=tpp
    }
  }

  ### ------------------------
  ### Do the geographic sampling
  ### ------------------------

  if(type=="geographic"){

    # Sample a regular grid of abence points with
    # n_presences x geodist_fact points
    abs=sampleRegular(rst,round(length(pres)*geodist_fact),sp=T)

    # create geo absences from geo_nrep times jittering regular samples
    geo.pts=list()
    for(i in 1:geo_nrep){

      pt.abs=abs@coords[,c("x","y")]
      pt.abs=apply(pt.abs,2,jitter,factor=3)

      model.idw <- geoIDW(p=as.data.frame(pres@coords),
                          a=as.data.frame(pt.abs))

      prd <- predict(rst, model.idw,mask=TRUE)

      # Sample cell centers proportional to interpolated presence
      # probability
      nonaval=which(!is.na(values(prd)))

      prb=round(values(prd)[nonaval]*10^4)

      smp=sample(1:length(prb),size=round(n*(1-add.strat+.1)/geo_nrep),prob=prb,replace=F)

      geo.pts[[i]]=coordinates(prd)[nonaval[smp],]
    }

    # combine
    df.pseu=SpatialPoints(do.call("rbind",geo.pts),
                          proj4string = crs(rst))

    # extract and subsample to match desired number
    sp_abs=as(df.pseu,"SpatialPointsDataFrame")
    sp_abs@data=as.data.frame(extract(env.stack,df.pseu))

    ### ------------------------
    ### Do the random sampling
    ### ------------------------
  } else if (type=="random"){

    nona=which(!is.na(values(rst)))

    rnd.pts=sample(nona,
                   size=n*(1-add.strat)*1.5,
                   prob=values(suppressWarnings(raster::area(rst)))[nona])
    crds=coordinates(rst)[rnd.pts,]
    sp_abs=SpatialPointsDataFrame(coords=crds,
                                  data=as.data.frame(extract(env.stack,crds)),
                                  proj4string =rst@crs)

    ### ------------------------
    ### Do the target.group sampling
    ### ------------------------
  } else if (type=="target.group"){

    if(is.na("target.group_dir")){
      stop("target.group_dir has to be supplied for sampling type target.group!")
    }

    fls=list.files(target.group_dir,full.names=T)
    ltarpt=lapply(fls,read.table,header=TRUE)
    dftarpt=do.call("rbind",ltarpt)
    sptarpt=SpatialPoints(dftarpt,
                          proj4string = crs(env.stk))

    if(force_spat_thin%in%c("absences","both")){

      if(nrow(dftarpt)>7*n*(1-add.strat)*1.1){

        crds=upsample_thin(sptarpt,limdist,n*(1-add.strat)*1.1)

      } else {

        crds=upsample_strategic(sptarpt,limdist,n*(1-add.strat)*1.1)
      }
    } else {
      if(nrow(dftarpt)>n*(1-add.strat)*1.1){
        crds=sptarpt[sample(1:nrow(dftarpt),n*(1-add.strat)*1.1,replace=FALSE),]

      } else {
        crds=sptarpt[sample(1:nrow(dftarpt),n*(1-add.strat)*1.1,replace=TRUE),]
        cat(paste("Less target.group points available than requested. Sampling with replacemnent..."))

      }

    }

    sp_abs=SpatialPointsDataFrame(coords=crds,
                                  data=as.data.frame(extract(env.stack,crds)),
                                  proj4string =rst@crs)

    ### ------------------------
    ### Do the geo.strat sampling
    ### ------------------------
  } else if (type=="geo.strat"){

    if(grepl("longlat",rst@crs)){

      # sample regularly on the surface of a sphere
      fglb=(extent(rst)@xmax-extent(rst)@xmin)*(extent(rst)@ymax-extent(rst)@ymin)/(360*180)
      N=round(1.1*n*(1-add.strat)*length(values(rst))/length(which(!is.na(values(rst))))/fglb)
      r=1
      Nc=0
      a=4*pi*r^2/N
      d=sqrt(a)
      Mx=round(pi/d)
      dx=pi/Mx
      dy=a/dx
      pts=matrix(NA,ncol=3,nrow=N*5)

      for(m in 0:(Mx-1)){
        xx=pi*(m+0.5)/Mx
        My=round(2*pi*sin(xx)/dy)
        for(nn in 0:(My-1)){
          yy=2*pi*nn/My
          pts[Nc+1,]=c(r*sin(xx)*cos(yy),r*sin(xx)*sin(yy),r*cos(xx))
          Nc=Nc+1
        }
      }
      pts=na.omit(pts)

      # Transform to longitude latitude
      lat=atan(sqrt(pts[,2]^2+pts[,1]^2)/pts[,3])
      lat=ifelse(lat<0,lat+max(lat,na.rm = TRUE),lat-max(lat,na.rm=TRUE))
      lon=atan(pts[,2]/pts[,1])
      regpt=na.omit(unique(cbind(lon,lat)))
      regpt[,1]=regpt[,1]/pi*360
      regpt[,2]=regpt[,2]/pi*180
      sppt=SpatialPoints(regpt,proj4string = crs(rst))

      # crude extraction
      xtr1=extract(rst,sppt)
      sppt=sppt[which(!is.na(xtr1)),]

    } else {
      Ntarg=round(1.1*n*ncell(rst)/length(which(!is.na(values(rst)))))
      sppt=sampleRegular(rst,Ntarg,sp=T)
      sppt=sppt[-which(is.na(sppt@data[,1])),]

    }

    # Preprare output
    sp_abs=SpatialPointsDataFrame(coords=sppt@coords,
                                  data=as.data.frame(extract(env.stack,sppt)),
                                  proj4string =rst@crs)

  } else if (type=="density"){

    ### ------------------------
    ### Prepare point pattern object
    ### ------------------------

    xrng=extent(rst)[1:2]
    yrng=extent(rst)[3:4]

    # Define Point Pattern object to calculate
    owi=owin(xrange=xrng,yrange=yrng)
    myppp=ppp(x=pres@coords[,1],y=pres@coords[,2],window = owi)

    ### ------------------------
    ### Generate 'im' object with density info
    ### ------------------------

    lo=dim(rst)[1:2]

    x=seq(xrng[1],xrng[2],length.out = lo[2])
    y=seq(yrng[1],yrng[2],length.out = lo[1])

    dens=density(myppp,xy=list(x=x,y=y),adjust=geodist_fact/10)

    ### ------------------------
    ### Draw locations proportional to point density
    ### ------------------------

    drst=raster(dens)
    drst=resample(drst,rst)
    drst=mask(drst,rst)

    vls=values(drst)

    # Replace NA's with zero probability
    if(any(is.na(vls)) || any(vls<0)){
      vls[which(is.na(vls) | vls<0)]=0
    }

    # Sample from density distributions
    pts=sample(1:length(vls),n*(1-add.strat)*1.1,prob=vls)

    sppt=SpatialPoints(coordinates(drst)[pts,])

    # Preprare output
    sp_abs=SpatialPointsDataFrame(coords=sppt@coords,
                                  data=as.data.frame(extract(env.stack,sppt)),
                                  proj4string =rst@crs)

  }

  if(exists("sp_abs")){
    # Subsample and remove NAs
    nosna=apply(sp_abs@data,1,function(x){
      all(!is.na(x))
    })

    sp_abs=sp_abs[which(nosna),]
    if(nrow(sp_abs)>round(n*(1-add.strat))){
      sp_abs=sp_abs[sample(1:nrow(sp_abs),round(n*(1-add.strat))),]
    }

  }

  ### ------------------------
  ### Do the env.strat sampling
  ### ------------------------

  if (grepl("env",type) | add.strat>0){

    # define env strat sampling strategy
    if(add.strat<1){
      tyyp="env.strat"
    } else {
      tyyp=type
    }

    # if no strata directory is supplied, just
    # calculate from scratch
    if(is.na(env.strat_path)){

      strpts = create_envstrat(env.stk=env.stack,
                               rAll=rAll,
                               save_it=FALSE,
                               strat_dir=NA,
                               type=tyyp)

      # if strat directory is supplied load if file exists
      # otherwise write file
      # if .RData file is supplied read directly
    } else if(grepl(".RData",env.strat_path)){
      load(tmfl)
    } else{

      tmfl=list.files(env.strat_path,pattern=tyyp,full.names=TRUE)
      lyn=paste(names(env.stack),collapse="_")
      tmfl=grep(lyn,tmfl,value=TRUE)

      if(length(tmfl)>0){

        load(tmfl)
        strpts=strpts[sample(1:nrow(strpts)),]

      } else {

        strpts = create_envstrat(env.stk=env.stack,
                                 rAll=rAll,
                                 save_it=TRUE,
                                 strat_dir=env.strat_path,
                                 type=tyyp)

      }
    }

    if(force_spat_thin%in%c("absences","both")){

      strfy=upsample_thin(strpts,limdist,n*add.strat*1.1)

    } else {
      if(type=="env.semi.strat"){
        strfy=strpts[sample(1:nrow(strpts),n*add.strat*1.1),]
      } else {
        strfy=stratify(strpts,tyyp,n*add.strat*1.1)
      }
    }

    sp_abse=SpatialPointsDataFrame(coords=coordinates(strfy),
                                   data=as.data.frame(extract(env.stack,strfy)),
                                   proj4string =strpts@proj4string)

    nosna=apply(sp_abse@data,1,function(x){
      all(!is.na(x))
    })

    sp_abse=sp_abse[which(nosna),]
    sp_abse=sp_abse[sample(1:nrow(sp_abse),round(n*add.strat)),]

    if(grepl("env",type)){
      sp_abs=sp_abse
    } else {
      sp_abs=rbind(sp_abse,sp_abs)
    }

  }

  ### ------------------------
  ### Prepare output
  ### ------------------------
  
  if(length(pres)>0){
    ptout=rbind(pres,as(sp_abs,"SpatialPoints"))
    
    out@pa=c(rep(1,length(pres)),
             rep(0,nrow(sp_abs)))
    
    out@env_vars=as.data.frame(rbind(xt_pres,sp_abs@data))
    out@xy=rbind(coordinates(pres),coordinates(sp_abs))
  } else {
    ptout=as(sp_abs,"SpatialPoints")
    out@pa=rep(0,nrow(sp_abs))
    out@env_vars=sp_abs@data
    out@xy=coordinates(sp_abs)
  }

  # return
  return(out)

}
