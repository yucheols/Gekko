##### ENM for Gekko swinhonis & G. japonicus == SDMtune-based workflow
set.seed(1234)

# load libraries
library(dplyr)
library(SDMtune)
library(ENMeval)
library(ggplot2)
library(raster)
library(terra)
library(ENMTools)
library(ecospat)

#####  PART 1 ::: env data  #####
# can download & prep WorldClim 2.1 data using the code below
#ClimDatDownloadR::WorldClim.HistClim.download(save.location = './envs', parameter = c('bio'), bio.var = c(1:19),
#                                              resolution = c('30s'), version.var = c('2.1'), clipping = T,
#                                              clip.extent = c(101.32, 146.26, 20.23, 45.37),
#                                              stacking.data = F, keep.raw.zip = T)

# ....or alternatively you can just load & crop envs 
envs <- raster::stack(list.files(path = 'envs', pattern = '.tif$', full.names = T))
envs <- dropLayer(envs, 'slope')
plot(envs[[1]])

names(envs) = c('altitude', 'bio1', 'bio4', 'bio12', 'dist_built', 'dist_forest')
print(envs)

#####  PART 2 ::: occurrences  #####
# collect occs == this can be done through megaSDM, rgbif, or by reading .csv files
#megaSDM::OccurrenceCollection(spplist = c('Gekko swinhonis', 'Gekko japonicus'),
#                              output = './occs',
#                              trainingarea = c(101.32, 146.26, 20.23, 45.37))

# read occurrence points
swin <- read.csv('swinhonis_10km.csv') %>% dplyr::select(1,3,2)   # G.swinhonis
japo <- read.csv('japonicus_10km.csv') %>% dplyr::select(1,3,2)   # G.japonicus

colnames(swin) = c('species','long','lat')
colnames(japo) = c('species','long','lat')

head(swin)
head(japo)

# spatially thin occurrence points == match occurrence dataset to predictor resolution 
#swin <- SDMtune::thinData(coords = swin[, c(2,3)], env = terra::rast(envs.s[[1]]), x = 'long', y = 'lat', progress = T)
#japo <- SDMtune::thinData(coords = japo[, c(2,3)], env = terra::rast(envs.j[[1]]), x = 'long', y = 'lat', progress = T)

# check
points(swin[, c(2,3)], col = 'red')
points(japo[, c(2,3)], col = 'blue')


#####  PART 3 ::: background  #####
# swinhonis
bg.s <- dismo::randomPoints(mask = envs[[1]], n = 10000, p = rbind(swin[, c(2,3)], japo[, c(2,3)]), excludep = T) %>% as.data.frame()
colnames(bg.s) = colnames(swin[, c(2,3)])
head(bg.s)

# japonicus
bg.j <- dismo::randomPoints(mask = envs[[1]], n = 10000, p = rbind(swin[, c(2,3)], japo[, c(2,3)]), excludep = T) %>% as.data.frame()
colnames(bg.j) = colnames(bg.s) 
head(bg.j)


#####  PART 4 ::: make SWD objects  #####
# swinhonis
s.swd <- prepareSWD(species = 'swinhonis', env = terra::rast(envs), p = swin[, c(2,3)], a = bg.s, verbose = T)

# japonicus
j.swd <- prepareSWD(species = 'japonicus', env = terra::rast(envs), p = japo[, c(2,3)], a = bg.j, verbose = T)


#####  PART 5 ::: make blocks  #####
# swinhonis
s.block <- get.block(occs = swin[, c(2,3)], bg = bg.s, orientation = 'lat_lon')

# japonicus
j.block <- get.block(occs = japo[, c(2,3)], bg = bg.j, orientation = 'lat_lon')


#####  PART 6 ::: make default CV model  #####
# swinhonis
s.cv.mod <- SDMtune::train(method = 'Maxent', data = s.swd, folds = s.block, progress = T, iter = 5000)   

# japonicus
j.cv.mod <- SDMtune::train(method = 'Maxent', data = j.swd, folds = j.block, progress = T, iter = 5000)


#####  PART 7 ::: model tuning  #####
# create a looper == test with spatial block
SDMtuner <- function(sp.names, models, hypers, metric, save, interactive, progress) {
  out.models <- list()
  out.metrics <- list()
  
  # init here
  for (i in 1:length(models)) {
    run <- gridSearch(model = models[[i]], hypers = hypers, metric = metric, 
                      save_models = save, interactive = interactive, progress = progress)
    
    out.models[[i]] <- run@models
    out.metrics[[i]] <- run@results
  }
  
  ### save
  # species 1
  out.models.sp1 <<- out.models[[1]]
  out.metrics.sp1 <<- out.metrics[[1]]
  
  # species 2
  out.models.sp2 <<- out.models[[2]]
  out.metrics.sp2 <<- out.metrics[[2]]
  
  # done
  return(paste('== the ENM tuning is done for', c(sp.names[[1]], sp.names[[2]]), '=='))
}

# define parameter combinations
hypers <- list(fc = c('l', 'q', 'h', 'p', 'lq', 'lp', 'qh', 'qp', 'hp', 'lqh', 'lqp', 'lqhp', 'lqhpt'),
               reg = seq(0.5, 5, by = 0.5))

# run function here
SDMtuner(sp.names = c('Gekko swinhonis', 'Gekko japonicus'), hypers = hypers, 
         models = list(s.cv.mod, j.cv.mod), metric = 'auc', save = T, interactive = F, progress = T)


#####  PART 8 ::: model selection  #####
print(out.metrics.sp1)
print(out.metrics.sp2)

#### calculate TSS
# make function
tssCalc <- function(models, test) {
  tssbin <- list()
  
  for (i in 1:length(models)) {
    tss <- SDMtune::tss(model = models[[i]], test = test)
    tssbin[[i]] <- tss
  }
  tss.out <<- as.data.frame(cbind(tssbin))
}

#### run
# swinhonis
swin.tss <- tssCalc(models = out.models.sp1, test = T)
print(swin.tss)

# japonicus
japo.tss <- tssCalc(models = out.models.sp2, test = T)
print(japo.tss)


#### bind metrics
# swinhonis 
swin.metric <- cbind(out.metrics.sp1, swin.tss)
swin.metric$tssbin <- as.numeric(swin.metric$tssbin)
head(swin.metric)

# japonicus
japo.metric <- cbind(out.metrics.sp2, japo.tss)
japo.metric$tssbin <- as.numeric(japo.metric$tssbin)
head(japo.metric)


#### min diff_AUC & max test_AUC
# species 1 ::: Gekko swinhonis == LQHP 5 == 129th model
(opt.sp1 <- swin.metric %>% dplyr::filter(diff_AUC == min(diff_AUC)) %>%
    dplyr::filter(test_AUC == max(test_AUC)))
  

# species 2 ::: Gekko japonicus == H 5 == 120th model
(opt.sp2 <- japo.metric %>% dplyr::filter(diff_AUC == min(diff_AUC)) %>%
    dplyr::filter(test_AUC == max(test_AUC))) 
    


#####  PART 9 ::: variable importance & thresholds #####
### Variable Importance
# swinhonis
(s.varimp <- varImp(model = out.models.sp1[[129]], permut = 10, progress = T))
write.csv(s.varimp, 'varimp/swinhonis.csv')

# japonicus
(j.varimp <- varImp(model = out.models.sp2[[120]], permut = 10, progress = T))
write.csv(j.varimp, 'varimp/japonicus.csv')

### Percent Contribution
# swinhonis
(s.percon <- maxentVarImp(model = out.models.sp1[[129]]) %>% dplyr::select(1,2))
write.csv(s.percon, 'varimp/swinhonis_percent_contribution.csv')

# japonicus
(j.percon <- maxentVarImp(model = out.models.sp2[[120]]) %>% dplyr::select(1,2)) 
write.csv(j.percon, 'varimp/japonicus_percent_contribution.csv')


### Thresh
# swinhonis
s.thresh <- SDMtune::thresholds(model = combineCV(out.models.sp1[[129]]), type = 'cloglog')
write.csv(s.thresh, 'thresh/swinhonis_thresholds.csv')

# japonicus
j.thresh <- SDMtune::thresholds(model = combineCV(out.models.sp2[[120]]), type = 'cloglog')
write.csv(j.thresh, 'thresh/japonicus_thresholds.csv')


#####  PART 10 ::: response curves  #####
# build a function to pull response curve data from the SDMtune function [ plotResponse ]
respDataPull <- function(model, var, type, only_presence, marginal, species_name) {
  
  plotdata.list <- list()
  
  for (i in 1:length(var)) {
    plotdata <- plotResponse(model = model, var = var[[i]], type = type, only_presence = only_presence, marginal = marginal)
    plotdata <- ggplot2::ggplot_build(plotdata)$data
    plotdata <- plotdata[[1]]
    
    plotdata <- plotdata[, c(1:4)]
    plotdata$species <- species_name
    plotdata$var <- var[[i]]
    
    plotdata.list[[i]] <- plotdata
  }
  plotdata.df <<- dplyr::bind_rows(plotdata.list) 
}

### G.swinhonis
# pull data
swin.resp.data <- respDataPull(model = out.models.sp1[[129]], var = names(envs), type = 'cloglog', 
                               only_presence = T, marginal = T, species_name = 'swinhonis')

print(swin.resp.data)

### G.japonicus
# pull data
japo.resp.data <- respDataPull(model = out.models.sp2[[120]], var = names(envs), type = 'cloglog',
                               only_presence = T, marginal = T, species_name = 'japonicus')

print(japo.resp.data)

### customize response curve plot for publication
# set font
windowsFonts(a = windowsFont(family = 'Times New Roman'))

###### variable value extracted from actual occurrence data to compute minmax range per variable 
### swinhonis
var_swin <- read.csv('swinhonis_var_minmax.csv') %>% dplyr::select('bio1','bio4','bio12','dist_forest','dist_built','altitude')
head(var_swin)

# min
var_swin_min <- data.frame(var = c('bio1','bio4','bio12','dist_forest','dist_built','altitude'),
                           val = c(min(var_swin$bio1), min(var_swin$bio4*100), min(var_swin$bio12), 
                                   min(var_swin$dist_forest*0.008333),min(var_swin$dist_built*0.008333), 
                                   min(var_swin$altitude)))

head(var_swin_min)

# max
var_swin_max <- data.frame(var = c('bio1','bio4','bio12','dist_forest','dist_built','altitude'),
                           val = c(max(var_swin$bio1), max(var_swin$bio4*100), max(var_swin$bio12), 
                                   max(var_swin$dist_forest*0.008333),max(var_swin$dist_built*0.008333), 
                                   max(var_swin$altitude)))


head(var_swin_max)


### japonicus
var_japo <- read.csv('japonicus_var_minmax.csv') %>% dplyr::select('bio1','bio4','bio12','dist_forest','dist_built','altitude')
head(var_japo)

# min
var_japo_min <- data.frame(var =  c('bio1','bio4','bio12','dist_forest','dist_built','altitude'),
                           val = c(min(var_japo$bio1), min(var_japo$bio4*100), min(var_japo$bio12), 
                                   min(var_japo$dist_forest*0.008333), min(var_japo$dist_built*0.008333), 
                                   min(var_japo$altitude)))

head(var_japo_min)

# max
var_japo_max <- data.frame(var =  c('bio1','bio4','bio12','dist_forest','dist_built','altitude'),
                           val = c(max(var_japo$bio1), max(var_japo$bio4*100), max(var_japo$bio12), 
                                   max(var_japo$dist_forest*0.008333), max(var_japo$dist_built*0.008333), 
                                   max(var_japo$altitude)))

head(var_japo_max)


######  plot
# set theme
theme <-  theme(text = element_text(family = 'a'),
                axis.title = element_text(size = 16, face = 'bold'),
                axis.title.x = element_text(margin = margin(t = 20)),
                axis.title.y = element_text(margin = margin(r = 20)),
                axis.text = element_text(size = 12),
                strip.text = element_text(size = 14))

## G.swinhonis
swin.resp.data %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(linewidth = 1.2, color = 'red') +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = 'grey', alpha = 0.4) +
  facet_wrap(~ var, scales = 'free') +
  geom_vline(data = filter(var_swin_min, var == 'bio1'), aes(xintercept = var_swin_min[1,2]), 
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_max, var == 'bio1'), aes(xintercept = var_swin_max[1,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_min, var == 'bio4'), aes(xintercept = var_swin_min[2,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_max, var == 'bio4'), aes(xintercept = var_swin_max[2,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_min, var == 'bio12'), aes(xintercept = var_swin_min[3,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_max, var == 'bio12'), aes(xintercept = var_swin_max[3,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_min, var == 'dist_forest'), aes(xintercept = var_swin_min[4,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_max, var == 'dist_forest'), aes(xintercept = var_swin_max[4,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_min, var == 'dist_built'), aes(xintercept = var_swin_min[5,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_max, var == 'dist_built'), aes(xintercept = var_swin_max[5,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_min, var == 'altitude'), aes(xintercept = var_swin_min[6,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_swin_max, var == 'altitude'), aes(xintercept = var_swin_max[6,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  xlab('Value') + ylab('Suitability') +
  theme_bw() +
  theme

## G.japonicus
japo.resp.data %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(linewidth = 1.2, color = 'red') +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = 'grey', alpha = 0.4) +
  facet_wrap(~ var, scales = 'free') +
  geom_vline(data = filter(var_japo_min, var == 'bio1'), aes(xintercept = var_japo_min[1,2]), 
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_max, var == 'bio1'), aes(xintercept = var_japo_max[1,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_min, var == 'bio4'), aes(xintercept = var_japo_min[2,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_max, var == 'bio4'), aes(xintercept = var_japo_max[2,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_min, var == 'bio12'), aes(xintercept = var_japo_min[3,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_max, var == 'bio12'), aes(xintercept = var_japo_max[3,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_min, var == 'dist_forest'), aes(xintercept = var_japo_min[4,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_max, var == 'dist_forest'), aes(xintercept = var_japo_max[4,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_min, var == 'dist_built'), aes(xintercept = var_japo_min[5,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_max, var == 'dist_built'), aes(xintercept = var_japo_max[5,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_min, var == 'altitude'), aes(xintercept = var_japo_min[6,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  geom_vline(data = filter(var_japo_max, var == 'altitude'), aes(xintercept = var_japo_max[6,2]),
             color = 'black', linewidth = 1.6, linetype = 3) +
  xlab('Value') + ylab('Suitability') +
  theme_bw() +
  theme


## show both sp
resp.data.both <- rbind(swin.resp.data, japo.resp.data)
colnames(resp.data.both) = c('x', 'y', 'ymin', 'ymax', 'Species', 'var')
head(resp.data.both)

# recode & order species names
resp.data.both$Species = recode_factor(resp.data.both$Species, 'swinhonis' = 'G.swinhonis', 'japonicus' = 'G.japonicus')
resp.data.both$Species = factor(resp.data.both$Species, levels = c('G.swinhonis', 'G.japonicus'))

# order variables
resp.data.both$var = factor(resp.data.both$var, levels = c('bio1', 'bio4', 'bio12', 'altitude', 'dist_built', 'dist_forest'))

# plot
resp.data.both %>%
  ggplot(aes(x = x, y = y, group = Species, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c('indianred1', 'cornflowerblue')) +
  facet_wrap(~ var, scales = 'free') +
  xlab('Value') + ylab('Suitability') +
  theme_bw() +
  theme(text = element_text(family = 'a'),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')

# save
ggsave('resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


#####  PART 11 ::: model prediction  #####
#### current
# G.swinhonis
swin.pred <- SDMtune::predict(object = out.models.sp1[[129]], data = terra::rast(envs), 
                              type = 'cloglog', progress = T, clamp = T) %>% raster()
plot(swin.pred)

# G.japonicus
japo.pred <- SDMtune::predict(object = out.models.sp2[[120]], data = terra::rast(envs),
                              type = 'cloglog', clamp = T, progress = T) %>% raster()
plot(japo.pred)


#### ssp == HadGem3_GC31-LL_ssp245 ::: year 2050
#proj <- brick('proj/wc2.1_30s_bioc_MIROC6_ssp245_2021-2040.tif')

# layer extraction
#for (i in 1:nlayers(proj)) {
#  r <- proj[[i]]
#  layer <- paste0('proj/layers/', 'bio', i, '.tif')
#  writeRaster(r, filename = layer, overwrite = T)
#}

# import layers & rename
proj.clim <- raster::stack(list.files(path = '2050future', pattern = '.tif$', full.names = T))
print(proj.clim)
plot(proj.clim[[1]])

names(proj.clim) = c('bio1', 'bio4', 'bio12')

# add non-climatic variables & subset
proj.clim <- raster::stack(proj.clim, envs[[c('altitude', 'dist_built', 'dist_forest')]])  # add non-climate variable
print(proj.clim)


### check layer units
unitCheck <- function(ref.env, proj.env) {
  
  # results
  results <- list()
  
  # names
  name.ref.env <- sort(names(ref.env))
  name.proj.env <- sort(names(proj.env[[names(ref.env)]]))
  
  # units
  MinMax <- for (i in 1:length(name.ref.env)) {
    ref.val <- as.data.frame(ref.env) %>% na.omit()
    proj.val <- as.data.frame(proj.env) %>% na.omit()
    
    ref.min <- min(ref.val[[name.ref.env[[i]]]])
    ref.avg <- mean(ref.val[[name.ref.env[[i]]]])
    ref.max <- max(ref.val[[name.ref.env[[i]]]])
    proj.min <- min(proj.val[[name.proj.env[[i]]]])
    proj.avg <- mean(proj.val[[name.proj.env[[i]]]])
    proj.max <- max(proj.val[[name.proj.env[[i]]]])
    
    
    results[[i]] <- data.frame(ref.min, ref.avg, ref.max, proj.min, proj.avg, proj.max)
  }
  
  results <- do.call(rbind, results)
  results$var.name <- name.ref.env
  print(results)
}

# check 
unitCheck(ref.env = envs, proj.env = proj.clim)  

###  adjust and run the codes below if needed
#mult <- raster::stack(subset(proj.clim, c('var names to multiply')))
#mult <- mult*10

#proj.clim <- raster::stack(mult, proj.clim[['var name that was not multiplied']])

#print(envs)
#print(proj.clim)


##### predict to future climate
# G.swinhonis
swin.proj <- SDMtune::predict(object = out.models.sp1[[129]], data = terra::rast(proj.clim), 
                              type = 'cloglog', progress = T, clamp = T) %>% raster()
plot(swin.proj)

# G.japonicus
japo.proj <- SDMtune::predict(object = out.models.sp2[[120]], data = terra::rast(proj.clim),
                              type = 'cloglog', progress = T, clamp = T) %>% raster()
plot(japo.proj)


#####  PART 11a ::: make binary map  #####
# calculate 10% presence threshold
calc_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}

######  use 10p & MTSS for each sp.
### G.swinhonis
# 10p
swin.10p <- calc_threshold(sdm = swin.pred, occs = swin[, c(2,3)], type = 'p10', binary = F)
swin.bin.10p <- ecospat::ecospat.binary.model(Pred = swin.pred, Threshold = minValue(swin.10p))
plot(swin.bin.10p)

# MTSS
swin.bin.mtss <- ecospat::ecospat.binary.model(Pred = swin.pred, Threshold = s.thresh[3,2])
plot(swin.bin.mtss)

# MTSS == future
swin.bin.mtss.proj <- ecospat::ecospat.binary.model(Pred = swin.proj, Threshold = s.thresh[3,2])
plot(swin.bin.mtss.proj)

### G.japonicus
# 10p
japo.10p <- calc_threshold(sdm = japo.pred, occs = japo[, c(2,3)], type = 'p10', binary = F)
japo.bin.10p <- ecospat::ecospat.binary.model(Pred = japo.pred, Threshold = minValue(japo.10p))
plot(japo.bin.10p)

# MTSS
japo.bin.mtss <- ecospat::ecospat.binary.model(Pred = japo.pred, Threshold = j.thresh[3,2])
plot(japo.bin.mtss)

# MTSS == future
japo.bin.mtss.proj <- ecospat::ecospat.binary.model(Pred = japo.proj, Threshold = j.thresh[3,2])
plot(japo.bin.mtss.proj)

#####  PART 11b ::: export pred rasters for GIS vis  #####
# G.swinhonis
swin.full <- raster::stack(swin.pred, swin.bin.10p, swin.bin.mtss, swin.proj)
names(swin.full) = c('swin_current', 'swin_bin_10p', 'swin_bin_mtss', 'swin_ssp245_2050') 

for (i in 1:nlayers(swin.full)) {
  r <- swin.full[[i]]
  layer <- paste0('output/swinhonis/full/', names(swin.full)[[i]], '.tif')
  writeRaster(r, filename = layer, overwrite = T)
}

# G.japonicus
japo.full <- raster::stack(japo.pred, japo.bin.10p, japo.bin.mtss, japo.proj)
names(japo.full) = c('japo_current', 'japo_bin_10p', 'japo_bin_mtss', 'japo_ssp245_2050') 

for (i in 1:nlayers(japo.full)) {
  r <- japo.full[[i]]
  layer <- paste0('output/japonicus/full/', names(japo.full)[[i]], '.tif')
  writeRaster(r, filename = layer, overwrite = T)
}


#####  PART 12 ::: model evaluation  #####
#### calc TSS == use this when ENMeval is used to make models
# swinhonis == 0.567
#pred.val.s <- raster::extract(swin.pred, rbind(swin, bg.s)) %>% as.data.frame()
#colnames(pred.val.s) = 'pred.val.s'
#head(pred.val.s)

#sp.occ.s <- rbind(swin, bg.s)
#sp.occ.s$pa <- c(rep(1, nrow(swin)), rep(0, nrow(bg.s)))
#head(sp.occ.s)

#swin.TSS <- ecospat::ecospat.max.tss(Pred = pred.val.s, Sp.occ = sp.occ.s$pa)
#print(swin.TSS)


# japonicus == 0.704
#pred.val.j <- raster::extract(japo.pred, rbind(japo, bg.j)) %>% as.data.frame()
#colnames(pred.val.j) = 'pred.val.j'
#head(pred.val.j)

#sp.occ.j <- rbind(japo, bg.s)
#sp.occ.j$pa <- c(rep(1, nrow(japo)), rep(0, nrow(bg.j)))
#head(sp.occ.j)

#japo.TSS <- ecospat.max.tss(Pred = pred.val.j, Sp.occ = sp.occ.j$pa)
#print(japo.TSS)


#####  PART 13 ::: MESS  #####
## G.swinhonis
swin.ref <- raster::extract(envs, swin[, c(2,3)]) %>% as.data.frame()
swin.mess <- dismo::mess(x = proj.clim, v = swin.ref, full = F, filename = 'output/swinhonis/mess/swin_mess.tif', overwrite = T)
plot(swin.mess)

## G.japonicus
japo.ref <- raster::extract(envs, japo[, c(2,3)]) %>% as.data.frame()
japo.mess <- dismo::mess(x = proj.clim, v = japo.ref, full = F, filename = 'output/japonicus/mess/japo_mess.tif', overwrite = T)
plot(japo.mess)


#####  PART 14 ::: optional == niche identity & background test in ecospat workflow  #####
#### generate MCPs for each sp. to define species-specific background
## make MCP
#swin.mcp <- ConR::EOO.computing(swin[, c(3,2,1)], write_shp = T)
#japo.mcp <- ConR::EOO.computing(japo[, c(3,2,1)], write_shp = T)

## import MCP
swin.mcp <- rgdal::readOGR('MCP/swinhonis_MCP/swinhonis_EOO_poly.shp')
japo.mcp <- rgdal::readOGR('MCP/japonicus_MCP/japonicus_EOO_poly.shp')

#### raster PCA
# PCA
env.glob.pca <- raster.pca(env = envs, n = 2)

# PC percentage
factoextra::get_eigenvalue(env.glob.pca$pca.object)

# PCA rasters for downstream analyses
env.glob <- env.glob.pca$rasters
plot(env.glob)

#### species 1 == G.swinhonis
# env
sp1.env <- raster::extract(env.glob, swin[, c(2,3)])
sp1.env <- cbind(rep('G.swinhonis', nrow(swin[, c(2,3)])), swin[, c(2,3)], sp1.env)
sp1.env <- sp1.env[complete.cases(sp1.env), ]
colnames(sp1.env) = c('Species', colnames(swin[, c(2,3)]), names(env.glob))
head(sp1.env)

# bg env
sp1.bg.env <- raster::mask(env.glob, swin.mcp)
sp1.bg.env <- as.data.frame(rasterToPoints(sp1.bg.env))
sp1.bg.env <- cbind(rep('G.swinhonis.bg', length(sp1.bg.env[,1])), sp1.bg.env)
sp1.bg.env <- sp1.bg.env[complete.cases(sp1.bg.env), ]
colnames(sp1.bg.env) = colnames(sp1.env)
head(sp1.bg.env)


#### species 2 == G.japonicus
# env
sp2.env <- raster::extract(env.glob, japo[, c(2,3)])
sp2.env <- cbind(rep('G.japonicus', nrow(japo[, c(2,3)])), japo[, c(2,3)], sp2.env)
sp2.env <- sp2.env[complete.cases(sp2.env), ]
colnames(sp2.env) = colnames(sp1.env)
head(sp2.env)

# bg env
sp2.bg.env <- raster::mask(env.glob, japo.mcp)
sp2.bg.env <- as.data.frame(rasterToPoints(sp2.bg.env))
sp2.bg.env <- cbind(rep('G.japonicus.bg', length(sp2.bg.env[,1])), sp2.bg.env)
sp2.bg.env <- sp2.bg.env[complete.cases(sp2.bg.env), ]
colnames(sp2.bg.env) = colnames(sp1.env)
head(sp2.bg.env)


#### background env == "global area"
total.bg.env <- as.data.frame(rasterToPoints(env.glob))
total.bg.env <- cbind(rep('background', length(total.bg.env[,1])), total.bg.env)
total.bg.env <- total.bg.env[complete.cases(total.bg.env), ]
colnames(total.bg.env) = colnames(sp1.env)
head(total.bg.env)


#### define niches
sp1.niche <- ecospat.grid.clim.dyn(glob = total.bg.env[, c(4,5)], glob1 = sp1.bg.env[, c(4,5)], sp = sp1.env[, c(4,5)], R = 1000)
sp2.niche <- ecospat.grid.clim.dyn(glob = total.bg.env[, c(4,5)], glob1 = sp2.bg.env[, c(4,5)], sp = sp2.env[, c(4,5)], R = 1000)

#### calculate overlaps
ecospat.niche.overlap(z1 = sp1.niche, z2 = sp2.niche, cor = T)

#### plot niche overlaps
## modify niche plotting function to pass additional arguments for manual axis adjustments
ecospat.plot.niche.dyn.mod <- function (z1, z2, quant = 0, title = "", name.axis1 = "Axis 1", 
                                        name.axis2 = "Axis 2", interest = 1, col.unf = "green", col.exp = "red", 
                                        col.stab = "blue", colZ1 = "green3", colZ2 = "red3", transparency = 70, ...) 
{
  t_col <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], alpha = (100 - 
                                                                percent) * 255/100, names = name, maxColorValue = 255)
  }
  col.unf = t_col(col.unf, transparency)
  col.exp = t_col(col.exp, transparency)
  col.stab = t_col(col.stab, transparency)
  if (is.null(z1$y)) {
    R <- length(z1$x)
    x <- z1$x
    xx <- sort(rep(1:length(x), 2))
    y1 <- z1$z.uncor/max(z1$z.uncor)
    Y1 <- z1$Z/max(z1$Z)
    if (quant > 0) {
      Y1.quant <- quantile(z1$Z[which(z1$Z > 0)], probs = seq(0, 
                                                              1, quant))[2]/max(z1$Z)
    }
    else {
      Y1.quant <- 0
    }
    Y1.quant <- Y1 - Y1.quant
    Y1.quant[Y1.quant < 0] <- 0
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 
                                           2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 
                                           2)]
    y2 <- z2$z.uncor/max(z2$z.uncor)
    Y2 <- z2$Z/max(z2$Z)
    if (quant > 0) {
      Y2.quant <- quantile(z2$Z[which(z2$Z > 0)], probs = seq(0, 
                                                              1, quant))[2]/max(z2$Z)
    }
    else {
      Y2.quant <- 0
    }
    Y2.quant <- Y2 - Y2.quant
    Y2.quant[Y2.quant < 0] <- 0
    yy2 <- sort(rep(1:length(y2), 2))[-c(1:2, length(y2) * 
                                           2)]
    YY2 <- sort(rep(1:length(Y2), 2))[-c(1:2, length(Y2) * 
                                           2)]
    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence", ...)
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = col.unf, border = 0)
    polygon(x[xx], c(0, y2[yy2], 0, 0), col = col.exp, border = 0)
    polygon(x[xx], c(0, apply(cbind(y2[yy2], y1[yy1]), 1, 
                              min, na.exclude = TRUE), 0, 0), col = col.stab, border = 0)
    lines(x[xx], c(0, Y2.quant[YY2], 0, 0), col = colZ2, 
          lty = "dashed")
    lines(x[xx], c(0, Y1.quant[YY1], 0, 0), col = colZ1, 
          lty = "dashed")
    lines(x[xx], c(0, Y2[YY2], 0, 0), col = colZ2)
    lines(x[xx], c(0, Y1[YY1], 0, 0), col = colZ1)
    segments(x0 = 0, y0 = 0, x1 = max(x[xx]), y1 = 0, col = "white")
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, col = "white")
    seg.cat <- function(inter, cat, col.unf, col.exp, col.stab) {
      if (inter[3] == 0) {
        my.col <- 0
      }
      if (inter[3] == 1) {
        my.col <- col.unf
      }
      if (inter[3] == 2) {
        my.col <- col.stab
      }
      if (inter[3] == -1) {
        my.col <- col.exp
      }
      segments(x0 = inter[1], y0 = -0.01, y1 = -0.01, x1 = inter[2], 
               col = my.col, lwd = 4, lty = 2)
    }
    cat <- ecospat.niche.dyn.index(z1, z2, intersection = quant)$dyn
    cat <- cat[length(cat):1]
    inter <- cbind(z1$x[-length(z1$x)], z1$x[-1], cat[-1])
    apply(inter, 1, seg.cat, col.unf = col.unf, col.exp = col.exp, 
          col.stab = col.stab)
  }
  if (!is.null(z1$y)) {
    if (interest == 1) {
      plot(z1$z.uncor, col = gray(100:0/100), legend = F, 
           xlab = name.axis1, ylab = name.axis2, ...)
    }
    if (interest == 2) {
      plot(z2$z.uncor, col = gray(100:0/100), legend = F, 
           xlab = name.axis1, ylab = name.axis2, ...)
    }
    raster::image(2 * z1$w + z2$w, col = c("#FFFFFF", col.exp, 
                                                    col.unf, col.stab), add = TRUE, legend = FALSE)
    title(title)
    raster::contour(z1$Z, add = TRUE, levels = quantile(z1$Z[z1$Z > 
                                                               0], c(0, quant)), drawlabels = FALSE, lty = c(1, 
                                                                                                             2), col = colZ1)
    raster::contour(z2$Z, add = TRUE, levels = quantile(z2$Z[z2$Z > 
                                                               0], c(0, quant)), drawlabels = FALSE, lty = c(1, 
                                                                                                             2), col = colZ2)
    axis(1, labels = FALSE, lwd.ticks = 0)
    axis(2, lwd.ticks = 0, labels = FALSE)
    axis(3, labels = FALSE, lwd.ticks = 0)
    axis(4, lwd.ticks = 0, labels = FALSE)
  }
}
  
 
## export high res plot
# open device
jpeg('density.jpeg', width = 30, height = 30, units = 'cm', res = 600)
tiff('density.tiff', width = 30, height = 30, units = 'cm', res = 600)

# make plot
ecospat.plot.niche.dyn.mod(z1 = sp1.niche, z2 = sp2.niche, name.axis1 = 'PC1', name.axis2 = 'PC2', interest = 2, 
                           quant = 0.5, xlim = c(-3,4), ylim = c(-2, 1))

# close device
dev.off()

##### use ENMTools to run ecospat identity test and ecospat background test
## define range
sp1.range <- mask(envs, swin.mcp)
sp2.range <- mask(envs, japo.mcp)

## create enmtools.species objects
swinhonis <- enmtools.species(range = sp1.range[[1]], presence.points = swin[, c(2,3)], species.name = 'swinhonis')
japonicus <- enmtools.species(range = sp2.range[[1]], presence.points = japo[, c(2,3)], species.name = 'japonicus')

#### niche equivalency test
niche.eq <- enmtools.ecospat.id(species.1 = swinhonis, species.2 = japonicus, env = envs, 
                                nreps = 1000, R = 1000, bg.source = 'range', verbose = T)

print(niche.eq)

#### niche similarity test ::: swinhonis niche vs japonicus background 
niche.sim <- enmtools.ecospat.bg(species.1 = swinhonis, species.2 = japonicus, env = envs,
                                 nreps = 1000, test.type = 'symmetric', R = 1000, bg.source = 'range', verbose = T)

print(niche.sim)


######## plot results
head(niche.eq$test.results$sim)
head(niche.sim$test.results$sim)

#### D
# prep data
niche.eq.D <- dplyr::select(niche.eq$test.results$sim, 1)
niche.eq.D$test = 'Equivalency'
head(niche.eq.D)

niche.sim.D <- dplyr::select(niche.sim$test.results$sim, 1)
niche.sim.D$test = 'Similarity'
head(niche.sim.D)

D.data <- rbind(niche.eq.D, niche.sim.D)
head(D.data)
tail(D.data)

# plot
D.data %>%
  ggplot(aes(x = D)) +
  facet_wrap(~ test, scales = 'free') +
  geom_density(fill = "#69b3a2", color = NA, alpha = 0.4) +
  geom_vline(xintercept = niche.eq$test.results$obs$D, linetype = 2, linewidth = 1.2) +
  ylab('Density') + theme_bw() + 
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)))

## export 
ggsave('niche.D.jpg', width = 20, height = 10, dpi = 600, units = 'cm')


#### I 
# prep data
niche.eq.I <- dplyr::select(niche.eq$test.results$sim, 2)
niche.eq.I$test <- 'Equivalency'
head(niche.eq.I)

niche.sim.I <- dplyr::select(niche.sim$test.results$sim, 2)
niche.sim.I$test <- 'Similarity'
head(niche.sim.I)

I.data <- rbind(niche.eq.I, niche.sim.I)
head(I.data)
tail(I.data)

I.data %>%
  ggplot(aes(x = I)) +
  facet_wrap(~ test, scales = 'free') +
  geom_density(fill = "#69b3a2", color = NA, alpha = 0.4) +
  geom_vline(xintercept = niche.eq$test.results$obs$I, linetype = 2, linewidth = 1.2) +
  ylab('Density') + theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)))

## export 
ggsave('niche.I.jpg', width = 20, height = 10, dpi = 600, units = 'cm')


