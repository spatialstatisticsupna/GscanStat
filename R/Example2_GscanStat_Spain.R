rm(list=ls())
library(bigDM)
library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
library(gridExtra)
library(INLA)
library(parallel)
library(purrr)
# library(rsatscan)
library(RColorBrewer)
library(sf)
library(sp)
library(spdep)
library(tibble)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


####################################
## Load data and cartography file ##
####################################
load("../Data/DataSpain_Cancer.Rdata")

str(data)
head(carto)

tmap_mode("view")
tm_shape(carto) + 
  tm_polygons() + 
  tm_basemap("OpenStreetMap")


## Compute spatial binary adjacency matrix
W <- bigDM::connect_subgraphs(carto)$W


## No. of areas and time points
nAreas <- length(unique(data$ID))
nYears <- length(unique(data$year))


#################################
## Run the GscanStat algorithm ##
#################################
source("./GscanStat_parallelClusteringAlgorithm.R")

## Parameters for the searching window
# spatialWindow.maxSize: maximum size for the spatial window
# temporalWindow.maxSize: maximum size for the temporal window
# reps.mc: no. of repetitions in the Monte Carlo test


## CAUTION: Computationally very intensive!!!
run <- FALSE

if(run){
  GScanStat_output <- run_GscanStat(data=data, carto=carto, W=W,
                                    spatialWindow.maxSize=2000,
                                    temporalWindow.maxSize=5,
                                    reps.mc=99,
                                    pVal_threshold=0.05)
  print(GScanStat_output)
  
  clusters_GscanStat <- get_clusterIDs_GscanStat(data, GScanStat_output)
}else{
  load("./Example2_GscanStat_Spain.Rdata")
}
str(clusters_GscanStat)


## Plot of the clusters identified by the GscanStat algorithm 
cluster.ID <- Reduce("+",lapply(1:ncol(clusters_GscanStat), function(x) x*as.numeric(as.character(clusters_GscanStat[,x]))))
cluster.ID <- as.data.frame(matrix(cluster.ID,nAreas, nYears, byrow=F))
colnames(cluster.ID) <- paste("Year", 1:nYears, sep=".")

carto.aux <- cbind(carto,cluster.ID)

Map.clusters <- tm_shape(carto.aux) +
  tm_polygons(fill=paste("Year", 1:nYears, sep="."),
              fill.scale = tm_scale_categorical(labels=c("No cluster","Cluster1: low-risk","Cluster2: high-risk"),
                                                values=c("lightgray","lightblue","red")),
              fill.legend=tm_legend("", show=TRUE, position=tm_pos_out("right","center")),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="GscanStat algorithm") + 
  tm_facets(nrow=3, ncol=3)

tmap_mode("plot")
print(Map.clusters)


##########################################################
## Fit local spatio-temporal models using bigDM package ##
##########################################################
help("STCAR_INLA", package="bigDM")

## We fit 1st-order nb models dividing the map of Spain into 15 Autonomous Communities
tm_shape(carto) + 
  tm_polygons(fill="CCAA",
              fill.scale=tm_scale(values="brewer.spectral"),
              fill.legend=tm_legend(show=FALSE))

## No-cluster model
model_noCluster <- STCAR_INLA(carto=carto, data=data,
                              ID.area="ID", ID.year="year", O="obs", E="exp",
                              ID.group="CCAA", model="partition", k=1,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              inla.mode="compact", compute.fitted.values=TRUE)

## Model including the clustering structure obtained through the GscanStat algorithm
model_GscanStat <- STCAR_INLA(carto=carto, data=data,
                              ID.area="ID", ID.year="year", O="obs", E="exp",
                              X=clusters_GscanStat,
                              ID.group="CCAA", model="partition", k=1,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              inla.mode="compact", compute.fitted.values=TRUE)

model_GscanStat$summary.fixed


######################################
## Compute model selection criteria ##
######################################
compute.MSC <- function(model){
  data.frame(deviance=model$dic$mean.deviance,
             pD=model$dic$p.eff,
             DIC=model$dic$dic,
             WAIC=model$waic$waic)
}

MODELS <- list('noCluster'=model_noCluster,
               'GscanStat'=model_GscanStat)

do.call(rbind,lapply(MODELS, compute.MSC))


##########################################################
## Maps of posterior median estimates of relative risks ##
##########################################################
time.periods <- paste("Year", 1:nYears, sep=".")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.2,1.3,Inf)


## noCluster model
####################
model_noCluster.RR <- matrix(model_noCluster$summary.fitted.values$`0.5quant`,
                             nAreas, nYears, byrow=F)
colnames(model_noCluster.RR) <- time.periods

carto.noCluster <- cbind(carto, model_noCluster.RR)

Map.noCluster <- tm_shape(carto.noCluster) +
  tm_polygons(fill=time.periods,
              fill.scale=tm_scale(values=paleta, breaks=values, interval.closure="left"),
              fill.legend=tm_legend("Risk", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="noCluster Model") + 
  tm_facets(nrow=3, ncol=3)

print(Map.noCluster)


## GscanStat model
####################
model_GscanStat.RR <- matrix(model_GscanStat$summary.fitted.values$`0.5quant`,
                             nAreas, nYears, byrow=F)
colnames(model_GscanStat.RR) <- time.periods

carto.GscanStat <- cbind(carto, model_GscanStat.RR)

Map.GscanStat <- tm_shape(carto.GscanStat) +
  tm_polygons(fill=time.periods,
              fill.scale=tm_scale(values=paleta, breaks=values, interval.closure="left"),
              fill.legend=tm_legend("Risk", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="GscanStat Model") + 
  tm_facets(nrow=3, ncol=3)

print(Map.GscanStat)
