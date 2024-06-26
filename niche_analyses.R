#####  PART 14 ::: optional == niche identity & background test in ecospat workflow  #####
library(ENMTools)
library(ecospat)

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


