##### ENM for Gekko swinhonis & G. japonicus == SDMtune-based workflow

##### REVISED CODE === points from R Korea removed for both species

# set seed
set.seed(1234)

# clear env
rm(list = ls(all.names = T))
gc()

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
# read occurrence points
swin <- read.csv('revision1/swinhonis_filtered.csv')    # G.swinhonis
japo <- read.csv('revision1/japonicus_filtered.csv')    # G.japonicus

head(swin)
head(japo)

# spatially thin occurrence points == match occurrence dataset to predictor resolution 
#swin <- SDMtune::thinData(coords = swin[, c(2,3)], env = terra::rast(envs[[1]]), x = 'long', y = 'lat', progress = T)
#japo <- SDMtune::thinData(coords = japo[, c(2,3)], env = terra::rast(envs[[1]]), x = 'long', y = 'lat', progress = T)

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
### reorganize occs & bg
# swinhonis
s.coords <- cbind(s.swd@coords, s.swd@pa)
colnames(s.coords) = c('long', 'lat', 'pa')

s.occs <- s.coords %>% dplyr::filter(pa == 1)
s.bg <- s.coords %>% dplyr::filter(pa == 0)

# japonicus
j.coords <- cbind(j.swd@coords, j.swd@pa)
colnames(j.coords) = c('long', 'lat', 'pa')  
  
j.occs <- j.coords %>% dplyr::filter(pa == 1)  
j.bg <- j.coords %>% dplyr::filter(pa == 0)

# swinhonis
s.block <- get.block(occs = s.occs, bg = s.bg, orientation = 'lat_lon')

# japonicus
j.block <- get.block(occs = j.occs, bg = j.bg, orientation = 'lat_lon')


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
  out.models.sp1 <- out.models[[1]]
  out.metrics.sp1 <- out.metrics[[1]]
  
  # species 2
  out.models.sp2 <- out.models[[2]]
  out.metrics.sp2 <- out.metrics[[2]]
  
  # combine outputs
  output <- list(models.sp1 = out.models.sp1, metrics.sp1 = out.metrics.sp1, models.sp2 = out.models.sp2, metrics.sp2 = out.metrics.sp2)
  
  # done
  paste('== the ENM tuning is done for', c(sp.names[[1]], sp.names[[2]]), '==')
  return(output)
}

# define parameter combinations
hypers <- list(fc = c('l', 'q', 'h', 'p', 'lq', 'lp', 'qh', 'qp', 'hp', 'lqh', 'lqp', 'lqhp', 'lqhpt'),
               reg = seq(0.5, 5, by = 0.5))

# run function here
tune_models <- SDMtuner(sp.names = c('Gekko swinhonis', 'Gekko japonicus'), hypers = hypers, 
                        models = list(s.cv.mod, j.cv.mod), metric = 'auc', save = T, interactive = F, progress = T)

print(tune_models$metrics.sp1)
print(tune_models$metrics.sp2)

# save tuned models as an .rds file for later use
#saveRDS(tune_models, 'revision1/rds/model_tuning.rds')


#####  PART 8 ::: model selection  #####
#### calculate TSS
# make function
tssCalc <- function(models, test) {
  tssbin <- list()
  
  for (i in 1:length(models)) {
    tss <- SDMtune::tss(model = models[[i]], test = test)
    tssbin[[i]] <- tss
  }
  tss.out <- as.data.frame(cbind(tssbin))
  return(tss.out)
}

#### run
# swinhonis
swin.tss <- tssCalc(models = tune_models$models.sp1, test = T)
print(swin.tss)

# japonicus
japo.tss <- tssCalc(models = tune_models$models.sp2, test = T)
print(japo.tss)


#### bind metrics
# swinhonis 
swin.metric <- cbind(tune_models$metrics.sp1, swin.tss)
swin.metric$tssbin <- as.numeric(swin.metric$tssbin)
head(swin.metric)

# japonicus
japo.metric <- cbind(tune_models$metrics.sp2, japo.tss)
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
(s.varimp <- varImp(model = tune_models$models.sp1[[129]], permut = 10, progress = T))
write.csv(s.varimp, 'revision1/varimp/swinhonis.csv')

# japonicus
(j.varimp <- varImp(model = tune_models$models.sp2[[120]], permut = 10, progress = T))
write.csv(j.varimp, 'revision1/varimp/japonicus.csv')

### Percent Contribution
# swinhonis
(s.percon <- maxentVarImp(model = tune_models$models.sp1[[129]]) %>% dplyr::select(1,2))
write.csv(s.percon, 'revision1/varimp/swinhonis_percent_contribution.csv')

# japonicus
(j.percon <- maxentVarImp(model = tune_models$models.sp2[[120]]) %>% dplyr::select(1,2)) 
write.csv(j.percon, 'revision1/varimp/japonicus_percent_contribution.csv')


### Thresh
# swinhonis
s.thresh <- SDMtune::thresholds(model = combineCV(tune_models$models.sp1[[129]]), type = 'cloglog')
write.csv(s.thresh, 'revision1/thresh/swinhonis_thresholds.csv')

# japonicus
j.thresh <- SDMtune::thresholds(model = combineCV(tune_models$models.sp2[[120]]), type = 'cloglog')
write.csv(j.thresh, 'revision1/thresh/japonicus_thresholds.csv')


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
  plotdata.df <- dplyr::bind_rows(plotdata.list) 
  return(plotdata.df)
}

### G.swinhonis
# pull data
swin.resp.data <- respDataPull(model = tune_models$models.sp1[[129]], var = names(envs), type = 'cloglog', 
                               only_presence = T, marginal = T, species_name = 'swinhonis')

print(swin.resp.data)

### G.japonicus
# pull data
japo.resp.data <- respDataPull(model = tune_models$models.sp2[[120]], var = names(envs), type = 'cloglog',
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
ggsave('revision1/plots/resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


#####  PART 11 ::: model prediction  #####
#### current
# G.swinhonis
swin.pred <- SDMtune::predict(object = tune_models$models.sp1[[129]], data = terra::rast(envs), 
                              type = 'cloglog', progress = T, clamp = T) %>% raster()
plot(swin.pred)

# G.japonicus
japo.pred <- SDMtune::predict(object = tune_models$models.sp2[[120]], data = terra::rast(envs),
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
swin.proj <- SDMtune::predict(object = tune_models$models.sp1[[129]], data = terra::rast(proj.clim), 
                              type = 'cloglog', progress = T, clamp = T) %>% raster()
plot(swin.proj)

# G.japonicus
japo.proj <- SDMtune::predict(object = tune_models$models.sp2[[120]], data = terra::rast(proj.clim),
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
swin.bin.10p <- ecospat::ecospat.binary.model(Pred = terra::rast(swin.pred), Threshold = minValue(swin.10p)) %>% raster()
plot(swin.bin.10p)

# MTSS
swin.bin.mtss <- ecospat::ecospat.binary.model(Pred = terra::rast(swin.pred), Threshold = s.thresh[3,2]) %>% raster()
plot(swin.bin.mtss)

# MTSS == future
swin.bin.mtss.proj <- ecospat::ecospat.binary.model(Pred = terra::rast(swin.proj), Threshold = s.thresh[3,2]) %>% raster()
plot(swin.bin.mtss.proj)

### G.japonicus
# 10p
japo.10p <- calc_threshold(sdm = japo.pred, occs = japo[, c(2,3)], type = 'p10', binary = F)
japo.bin.10p <- ecospat::ecospat.binary.model(Pred = terra::rast(japo.pred), Threshold = minValue(japo.10p)) %>% raster()
plot(japo.bin.10p)

# MTSS
japo.bin.mtss <- ecospat::ecospat.binary.model(Pred = terra::rast(japo.pred), Threshold = j.thresh[3,2]) %>% raster()
plot(japo.bin.mtss)

# MTSS == future
japo.bin.mtss.proj <- ecospat::ecospat.binary.model(Pred = terra::rast(japo.proj), Threshold = j.thresh[3,2]) %>% raster()
plot(japo.bin.mtss.proj)

#####  PART 11b ::: export pred rasters for GIS vis  #####
# G.swinhonis
swin.full <- raster::stack(swin.pred, swin.bin.10p, swin.bin.mtss, swin.proj)
names(swin.full) = c('swin_current', 'swin_bin_10p', 'swin_bin_mtss', 'swin_ssp245_2050') 

for (i in 1:nlayers(swin.full)) {
  r <- swin.full[[i]]
  layer <- paste0('revision1/output/swinhonis/', names(swin.full)[[i]], '.tif')
  writeRaster(r, filename = layer, overwrite = T)
}

# G.japonicus
japo.full <- raster::stack(japo.pred, japo.bin.10p, japo.bin.mtss, japo.proj)
names(japo.full) = c('japo_current', 'japo_bin_10p', 'japo_bin_mtss', 'japo_ssp245_2050') 

for (i in 1:nlayers(japo.full)) {
  r <- japo.full[[i]]
  layer <- paste0('revision1/output/japonicus/', names(japo.full)[[i]], '.tif')
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
swin.mess <- dismo::mess(x = proj.clim, v = swin.ref, full = F, filename = 'revision1/output/swinhonis/swin_mess.tif', overwrite = T)
plot(swin.mess)

## G.japonicus
japo.ref <- raster::extract(envs, japo[, c(2,3)]) %>% as.data.frame()
japo.mess <- dismo::mess(x = proj.clim, v = japo.ref, full = F, filename = 'revision1/output/japonicus/japo_mess.tif', overwrite = T)
plot(japo.mess)


