### occ density grid prep
# env values of the entire study area
env.glob <- as.data.frame(na.omit(envs)) 
head(env.glob)

# env values of background area
env.bg <- raster::extract(envs, rbind(bg.s, bg.j)) %>% as.data.frame()
head(env.bg)

### swinhonis
# extract env values from occs
swin.env <- raster::extract(envs, swin[, c(2,3)]) %>% as.data.frame()
head(swin.env)

### japonicus
# extract env values from occs
japo.env <- raster::extract(envs, japo[, c(2,3)]) %>% as.data.frame()
head(japo.env)

### PCA-env
pca.env <- ade4::dudi.pca(df = rbind(env.bg, swin.env, japo.env), scannf = F, nf = 2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)

# PCA scores for the whole study area
scores.glob <- pca.env$li
head(scores.glob)

# PCA scores for bg
scores.bg <- ade4::suprow(pca.env, rbind(swin.env, japo.env))$li
head(scores.bg)

# PCA scores for swinhonis
scores.swin <- ade4::suprow(pca.env, swin.env)$li
head(scores.swin)

# PCA scores for japonicus
scores.japo <- ade4::suprow(pca.env, japo.env)$li
head(scores.japo)

### prep grid
swin.grd <- ecospat.grid.clim.dyn(glob = scores.glob, glob1 = scores.bg, sp = scores.swin, R = 1000)
japo.grd <- ecospat.grid.clim.dyn(glob = scores.glob, glob1 = scores.bg, sp = scores.japo, R = 1000)

### niche equivalency test
niche.eq <- ecospat.niche.equivalency.test(z1 = swin.grd, z2 = japo.grd, rep = 1000, intersection = NA,
                                           ncores = 2,
                                           overlap.alternative = 'higher',
                                           expansion.alternative = 'lower',
                                           stability.alternative = 'higher',
                                           unfilling.alternative = 'lower')

print(niche.eq)


### niche similarity test
niche.sim <- ecospat.niche.similarity.test(z1 = swin.grd, z2 = japo.grd, rep = 1000, intersection = NA,
                                           rand.type = 1, ncores = 2,
                                           overlap.alternative = 'higher',
                                           expansion.alternative = 'lower',
                                           stability.alternative = 'higher',
                                           unfilling.alternative = 'lower')

print(niche.sim)


### plot niche
ecospat.plot.niche.dyn(z1 = swin.grd, z2 = japo.grd, interest = 2, name.axis1 = 'PC1', name.axis2 = 'PC2')

### plot niche eq/sim statistics == Schoener's D
# format data
niche.eq.D <- niche.eq$sim %>% dplyr::select(1)
niche.sim.D <- niche.sim$sim %>% dplyr::select(1)

niche.eq.D$test <- 'Equivalency'
niche.sim.D$test <- 'Similarity'

niche.data <- rbind(niche.eq.D, niche.sim.D)
head(niche.data)
tail(niche.data)

# plot
niche.data %>%
  ggplot(aes(x = D)) +
  xlim(0, 0.8) +
  xlab("Schoener's D") + ylab('Density') +
  facet_wrap(~ test, scales = 'free') +
  geom_density(fill = "#69b3a2", color = NA, alpha = 0.4) +
  geom_vline(xintercept = niche.eq$obs$D, linetype = 2, linewidth = 1.1) +
  theme_bw() + 
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)))
