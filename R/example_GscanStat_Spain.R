rm(list=ls())
library(sf)
library(tmap)
library(tidyverse)
library(devtools)
library(bigDM)

# load data
load("exampleGscanStat_simulatedData_SP_AggregatedAreas.RData")

tm_shape(carto) + tm_polygons()

nAreas <- length(unique(data$muni))
nYears <- length(unique(data$year))
source("GscanStat_parallelClusteringAlgorithm.R")

###############################
## Using GscanStat algorithm ##
###############################
source("./GscanStat_parallelClusteringAlgorithm.R")



# parameter for the search window in the GscanStat algorithm:
# # maximum size for the spatial window: 2000
# # maximum size for the temporal window: 5
# # no. of repetitions in the Monte Carlo test: 99

GScanStat_output <- run_GscanStat(data=data, W=W, spatialWindow.maxSize=2000,
                                  temporalWindow.maxSize=5, reps.mc=99,
                                  pVal_threshold = 0.05)
clusters_GscanStat <- get_clusterIDs_GscanStat(data,GScanStat_output)
 




################################################################
## Fitting models with a divide-and-conquer approach with K=1 ##
################################################################
# The map is divided into 15 regions, corresponding to the
# autonomous communities of Spain. This information is contained
# in the `ID.CCAA` variable of the `carto` object.

tm_shape(carto) + tm_polygons(fill="ID.CCAA",
                              fill.scale=tm_scale(values="brewer.blues"),
                              fill.legend=tm_legend(show=FALSE))

# No-cluster Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction and no cluster effect
model_noCluster <- STCAR_INLA(carto=carto, data=data,
                              ID.area="ID.new", ID.year="year", ID.group = "ID.CCAA",
                              O="obs", E="exp",
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="partition", k=1, inla.mode="compact")

# GscanStat Model
## BYM2 priors for spatial effects, random-walk order 1 for temporal effect, 
## Type IV for space-time interaction and cluster effect given by GscanStat partition
# load("ID_Matrix_GscanStat.RData")
model_GscanStat <- STCAR_INLA(carto=carto, data=data,
                              ID.area="ID.new", ID.year="year", ID.group = "ID.CCAA",
                              O="obs", E="exp", X=clusters_GscanStat,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="partition", k=1, inla.mode="compact")
# data_GscanStat <- cbind(data,clusters_GscanStat)
# clusters_GscanStat_names <- names(clusters_GscanStat)
# model_GscanStat <- STCAR_INLA(carto=carto, data=data_GscanStat,
#                               ID.area="ID.new", ID.year="year", ID.group = "ID.CCAA",
#                               O="obs", E="exp", X=clusters_GscanStat_names,
#                               spatial="BYM2", temporal="rw1", interaction="TypeIV",
#                               model="partition", k=1, inla.mode="compact")

# Obtaining model fitting scores
# Note that the global LS score cannot be calculated when using a divide-and-conquer approach
scores.partition <- data.frame(deviance=c(model_noCluster$dic$mean.deviance,
                                          model_GscanStat$dic$mean.deviance),
                               pD=c(model_noCluster$dic$p.eff,
                                    model_GscanStat$dic$p.eff),
                               DIC=c(model_noCluster$dic$dic,
                                     model_GscanStat$dic$dic),
                               WAIC=c(model_noCluster$waic$waic,
                                      model_GscanStat$waic$waic)
)
row.names(scores.partition) <- c("noCluster","GscanStat")
print(scores.partition)

###############################################
## space-time risk Maps for partition models ##
###############################################
time.periods <- paste("Year",1:nYears,sep = ".")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.2,1.3,Inf)
## noCluster model
model_noCluster.risk <- matrix(exp(model_noCluster$summary.linear.predictor$`0.5quant`),nAreas,nYears,byrow=F)
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
model_GscanStat.risk <- matrix(exp(model_GscanStat$summary.linear.predictor$`0.5quant`),nAreas,nYears,byrow=F)
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
Map.GsacanStat <- tmap_arrange(c(maps, list(leyend)), ncol=3, nrow=3)
print(Map.GsacanStat)


