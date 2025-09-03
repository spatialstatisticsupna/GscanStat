rm(list=ls())
library(tmap)
library(data.table)
library(RColorBrewer)
library(INLA)
library(sp)
library(spdep)
library(bigDM)
library(gridExtra)

# load data
load("exampleGscanStat_simulatedData_NA.RData")

tm_shape(carto) + tm_polygons()
    
nAreas <- length(unique(data$muni))
nYears <- length(unique(data$year))

###############################
## Using GscanStat algorithm ##
###############################
source("./GscanStat_parallelClusteringAlgorithm.R")
# we obtain the area pair-wise distance matrix
geoDist <- as_Spatial(carto) %>%
  coordinates() %>%
  as.matrix %>%
  spDists(x = ., y = ., longlat = FALSE)


# parameter for the search window in the GscanStat algorithm:
# # maximum size for the spatial window: we search within the
# # `spatialNeighbourhoodSize` nearest spatial neighbours of each area
spatialNeighbourhoodSize <- 200
# # maximum size for the temporal window: we search within
# # +- `temporalNeighbourhoodSize` time periods
temporalNeighbourhoodSize <- 5
# # no. of repetitions in the Monte Carlo test
reps.mc <- 99

  
clusters <- getSignificantSpatioTemporalClusters_parallel(data=data, geoDist=geoDist,
                                                        W=W, maxAreas=spatialNeighbourhoodSize,
                                                        timeWindowSize=temporalNeighbourhoodSize,
                                                        reps.mc=reps.mc,
                                                        parallel=TRUE,
                                                        verbose=FALSE,
                                                        maxCPUs=detectCores()-1)
# The obtained clusters may contain overlapped partitions. We filter them out getting the most significant clusters
# with a minimum significance level of`pVal_threshold`
clusters_GscanStat <- get.nonOverlapping.clusters(clusters,pVal_threshold = 0.05)

# We obtain the cluster matrix used when fitting the models
clusterIdMatrix_GscanStat <- getClusterIdentityMatrix(data, clusters_GscanStat)

#############################
## Using SaTScan algorithm ##
#############################
# to use SaTScan is needed to install not only the `rsatscan` R-library
# but also SaTScan software (https://www.satscan.org/download.html)
source("./satScan_auxFunctions.R")
# Year format for SaTScan must be YYYY: as we have values 1:8, we transform
# them to 2001:2008
data$year <- data$year+2000
satScan_output <- run_SatScan(data, carto)
clusterIdMatrix_satScan <- get_clusteringIDMatrix_satScan(satScan_output, pVal.threshold=0.05)


##############################################################
## Fitting models with global model (the whole risk surface)##
##############################################################
      
# No-cluster Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction 
model_noCluster <- STCAR_INLA(carto=carto, data=data,
                              ID.area="muni", ID.year="year", O="obs", E="exp",
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="global", inla.mode="compact")

# GscanStat Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction 
# load("ID_Matrix_GscanStat.RData")
model_GscanStat <- STCAR_INLA(carto=carto, data=data, ID.area="muni", 
                              ID.year="year", O="obs", E="exp", X=clusterIdMatrix_GscanStat,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="global", inla.mode="compact")


# SaTScan Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction 
# load("ID_Matrix_satScan.RData")
model_satScan <- STCAR_INLA(carto=carto, data=data, ID.area="muni", 
                              ID.year="year", O="obs", E="exp", X=clusterIdMatrix_satScan,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="global", inla.mode="compact")

# Obtaining model fitting scores
scores.global <- data.frame(deviance=c(model_noCluster$dic$mean.deviance,
                                model_GscanStat$dic$mean.deviance,
                                model_satScan$dic$mean.deviance),
                     pD=c(model_noCluster$dic$p.eff,
                          model_GscanStat$dic$p.eff,
                          model_satScan$dic$p.eff),
                     DIC=c(model_noCluster$dic$dic,
                           model_GscanStat$dic$dic,
                           model_satScan$dic$dic),
                     WAIC=c(model_noCluster$waic$waic,
                            model_GscanStat$waic$waic,
                            model_satScan$waic$waic),
                     LS=c(-sum(log(model_noCluster$cpo$cpo)),
                          -sum(log(model_GscanStat$cpo$cpo)),
                          -sum(log(model_satScan$cpo$cpo)))
                     )
row.names(scores.global) <- c("noCluster","GscanStat","SaTScan")
print(scores.global)

############################################
## space-time risk Maps for global models ##
############################################
time.periods <- paste("Year",1:nYears,sep = ".")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.2,1.3,Inf)
## noCluster model
model_noCluster.risk <- matrix(model_noCluster$summary.fitted.values$`0.5quant`,nAreas,nYears,byrow=F)
colnames(model_noCluster.risk) <- time.periods
carto.aux <- cbind(carto,model_noCluster.risk)
# list of maps for each time period
maps <- lapply(time.periods, function(period) {
  tm_shape(carto.aux) +
    tm_polygons(col=period, palette=paleta,
                style="fixed", breaks=values, interval.closure="left",
                lwd=0, border.col=rgb(0, 0, 0, 0)) +
    tm_layout(legend.show=FALSE,
              main.title=period,
              main.title.position="center")
})
# Leyend in a different empty map
leyend <- tm_shape(carto.aux) +
  tm_polygons(col=time.periods[1], palette=paleta,
              style="fixed", breaks=values, interval.closure="left",
              border.col=NA, title="Risk - noCluter model") + 
  tm_layout(legend.only=TRUE,
            legend.position=c("center","center"),
            legend.text.size=1.2,           
            legend.title.size=1.4)          
# Combine the maps in a 3x3 grid
Map.noCluster <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.noCluster)

## GscanStat model
model_GscanStat.risk <- matrix(model_GscanStat$summary.fitted.values$`0.5quant`,nAreas,nYears,byrow=F)
colnames(model_GscanStat.risk) <- time.periods
carto.aux <- cbind(carto,model_GscanStat.risk)
# list of maps for each time period
maps <- lapply(time.periods, function(period) {
  tm_shape(carto.aux) +
    tm_polygons(col=period, palette=paleta,
                style="fixed", breaks=values, interval.closure="left",
                lwd=0, border.col=rgb(0, 0, 0, 0)) +
    tm_layout(legend.show=FALSE,
              main.title=period,
              main.title.position="center")
})
# Leyend in a different empty map
leyend <- tm_shape(carto.aux) +
  tm_polygons(col=time.periods[1], palette=paleta,
              style="fixed", breaks=values, interval.closure="left",
              border.col=NA, title="Risk - GscanStat model") + 
  tm_layout(legend.only=TRUE,
            legend.position=c("center","center"),
            legend.text.size=1.2,           
            legend.title.size=1.4)          
# Combine the maps in a 3x3 grid with the leyend
Map.noCluster <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.noCluster)

## satScam model
model_satScan.risk <- matrix(model_satScan$summary.fitted.values$`0.5quant`,nAreas,nYears,byrow=F)
colnames(model_satScan.risk) <- time.periods
carto.aux <- cbind(carto,model_satScan.risk)
# list of maps for each time period
maps <- lapply(time.periods, function(period) {
  tm_shape(carto.aux) +
    tm_polygons(col=period, palette=paleta,
                style="fixed", breaks=values, interval.closure="left",
                lwd=0, border.col=rgb(0, 0, 0, 0)) +
    tm_layout(legend.show=FALSE,
              main.title=period,
              main.title.position="center")
})
# Leyend in a different empty map
leyend <- tm_shape(carto.aux) +
  tm_polygons(col=time.periods[1], palette=paleta,
              style="fixed", breaks=values, interval.closure="left",
              border.col=NA, title="Risk - SaTScan model") + 
  tm_layout(legend.only=TRUE,
            legend.position=c("center","center"),
            legend.text.size=1.2,           
            legend.title.size=1.4)          
# Combine the maps in a 3x3 grid with the leyend
Map.noCluster <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.noCluster)


################################################################
## Fitting models with a divide-and-conquer approach with K=1 ##
################################################################
carto.divideAndConquer <- random_partition(carto=carto, rows=2, columns=2, max.size=NULL)

tm_shape(carto.divideAndConquer) + tm_polygons(fill="ID.group",
                                               fill.scale=tm_scale(values="brewer.blues"),
                                               fill.legend=tm_legend(show=FALSE))

# No-cluster Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction and no cluster effect
model_noCluster <- STCAR_INLA(carto=carto.divideAndConquer, data=data,
                              ID.area="muni", ID.year="year", ID.group = "ID.group",
                              O="obs", E="exp",
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="partition", k=1, inla.mode="compact")

# GscanStat Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction and cluster effect given by GscanStat partition
# load("ID_Matrix_GscanStat.RData")
model_GscanStat <- STCAR_INLA(carto=carto.divideAndConquer, data=data,
                              ID.area="muni", ID.year="year", ID.group = "ID.group",
                              O="obs", E="exp", X=clusterIdMatrix_GscanStat,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="partition", k=1, inla.mode="compact")


# SaTScan Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction and cluster effect given by SaTScan partition
# load("ID_Matrix_satScan.RData")
model_satScan <- STCAR_INLA(carto=carto.divideAndConquer, data=data,
                            ID.area="muni", ID.year="year", ID.group = "ID.group",
                            O="obs", E="exp", X=clusterIdMatrix_satScan,
                            spatial="BYM2", temporal="rw1", interaction="TypeIV",
                            model="partition", k=1, inla.mode="compact")

# Obtaining model fitting scores
# Note that the global LS score cannot be calculated when using a divide-and-conquer approach
scores.partition <- data.frame(deviance=c(model_noCluster$dic$mean.deviance,
                                          model_GscanStat$dic$mean.deviance,
                                          model_satScan$dic$mean.deviance),
                               pD=c(model_noCluster$dic$p.eff,
                                    model_GscanStat$dic$p.eff,
                                    model_satScan$dic$p.eff),
                               DIC=c(model_noCluster$dic$dic,
                                     model_GscanStat$dic$dic,
                                     model_satScan$dic$dic),
                               WAIC=c(model_noCluster$waic$waic,
                                      model_GscanStat$waic$waic,
                                      model_satScan$waic$waic)
                               )
row.names(scores.partition) <- c("noCluster","GscanStat","SaTScan")
print(scores.partition)

###############################################
## space-time risk Maps for partition models ##
###############################################
time.periods <- paste("Year",1:nYears,sep = ".")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.2,1.3,Inf)
## noCluster model
model_noCluster.risk <- matrix(model_noCluster$summary.fitted.values$`0.5quant`,nAreas,nYears,byrow=F)
colnames(model_noCluster.risk) <- time.periods
carto.aux <- cbind(carto,model_noCluster.risk)
# list of maps for each time period
maps <- lapply(time.periods, function(period) {
  tm_shape(carto.aux) +
    tm_polygons(col=period, palette=paleta,
                style="fixed", breaks=values, interval.closure="left",
                lwd=0, border.col=rgb(0, 0, 0, 0)) +
    tm_layout(legend.show=FALSE,
              main.title=period,
              main.title.position="center")
})
# Leyend in a different empty map
leyend <- tm_shape(carto.aux) +
  tm_polygons(col=time.periods[1], palette=paleta,
              style="fixed", breaks=values, interval.closure="left",
              border.col=NA, title="Risk - noCluter model") + 
  tm_layout(legend.only=TRUE,
            legend.position=c("center","center"),
            legend.text.size=1.2,           
            legend.title.size=1.4)          
# Combine the maps in a 3x3 grid
Map.noCluster <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.noCluster)

## GscanStat model
model_GscanStat.risk <- matrix(model_GscanStat$summary.fitted.values$`0.5quant`,nAreas,nYears,byrow=F)
colnames(model_GscanStat.risk) <- time.periods
carto.aux <- cbind(carto,model_GscanStat.risk)
# list of maps for each time period
maps <- lapply(time.periods, function(period) {
  tm_shape(carto.aux) +
    tm_polygons(col=period, palette=paleta,
                style="fixed", breaks=values, interval.closure="left",
                lwd=0, border.col=rgb(0, 0, 0, 0)) +
    tm_layout(legend.show=FALSE,
              main.title=period,
              main.title.position="center")
})
# Leyend in a different empty map
leyend <- tm_shape(carto.aux) +
  tm_polygons(col=time.periods[1], palette=paleta,
              style="fixed", breaks=values, interval.closure="left",
              border.col=NA, title="Risk - GscanStat model") + 
  tm_layout(legend.only=TRUE,
            legend.position=c("center","center"),
            legend.text.size=1.2,           
            legend.title.size=1.4)          
# Combine the maps in a 3x3 grid with the leyend
Map.noCluster <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.noCluster)

## satScam model
model_satScan.risk <- matrix(model_satScan$summary.fitted.values$`0.5quant`,nAreas,nYears,byrow=F)
colnames(model_satScan.risk) <- time.periods
carto.aux <- cbind(carto,model_satScan.risk)
# list of maps for each time period
maps <- lapply(time.periods, function(period) {
  tm_shape(carto.aux) +
    tm_polygons(col=period, palette=paleta,
                style="fixed", breaks=values, interval.closure="left",
                lwd=0, border.col=rgb(0, 0, 0, 0)) +
    tm_layout(legend.show=FALSE,
              main.title=period,
              main.title.position="center")
})
# Leyend in a different empty map
leyend <- tm_shape(carto.aux) +
  tm_polygons(col=time.periods[1], palette=paleta,
              style="fixed", breaks=values, interval.closure="left",
              border.col=NA, title="Risk - SaTScan model") + 
  tm_layout(legend.only=TRUE,
            legend.position=c("center","center"),
            legend.text.size=1.2,           
            legend.title.size=1.4)          
# Combine the maps in a 3x3 grid with the leyend
Map.noCluster <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.noCluster)
