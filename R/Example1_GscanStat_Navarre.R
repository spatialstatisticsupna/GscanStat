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
library(rsatscan)
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


## Select the spatial units corresponding to the Autonomous Region of Navarra
data.NAV <- data |> filter(substr(ID,1,2)=="31")
carto.NAV <- carto |> filter(substr(ID,1,2)=="31")

tmap_mode("view")
tm_shape(carto.NAV) + 
  tm_polygons() + 
  tm_basemap("OpenStreetMap")
tmap_mode("plot")


## Compute spatial binary adjacency matrix
W <- bigDM::connect_subgraphs(carto.NAV)$W


## No. of areas and time points
nAreas <- length(unique(data.NAV$ID))
nYears <- length(unique(data.NAV$year))


#################################
## Run the GscanStat algorithm ##
#################################
source("./GscanStat_parallelClusteringAlgorithm.R")

## Parameters for the searching window
# spatialWindow.maxSize: maximum size for the spatial window
# temporalWindow.maxSize: maximum size for the temporal window
# reps.mc: no. of repetitions in the Monte Carlo test

GScanStat_output <- run_GscanStat(data=data.NAV,
                                  carto=carto.NAV,
                                  W=W, 
                                  spatialWindow.maxSize=50,
                                  temporalWindow.maxSize=5,
                                  reps.mc=99,
                                  pVal_threshold=0.05)
print(GScanStat_output)

clusters_GscanStat <- get_clusterIDs_GscanStat(data.NAV, GScanStat_output)
str(clusters_GscanStat)


## Plot of the clusters identified by the GscanStat algorithm 
cluster.ID <- Reduce("+",lapply(1:ncol(clusters_GscanStat), function(x) x*as.numeric(as.character(clusters_GscanStat[,x]))))
cluster.ID <- as.data.frame(matrix(cluster.ID,nAreas, nYears, byrow=F))
colnames(cluster.ID) <- paste("Year", 1:nYears, sep=".")

carto.aux <- cbind(carto.NAV,cluster.ID)

tm_shape(carto.aux) +
  tm_polygons(fill=paste("Year", 1:nYears, sep="."),
              fill.scale = tm_scale_categorical(labels=c("No cluster",paste("Cluster",1:ncol(clusters_GscanStat))),
                                                values=c("lightgray", brewer.pal(ncol(clusters_GscanStat), "Set1"))),
              fill.legend=tm_legend("", show=TRUE, position=tm_pos_out("right","center")),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="GscanStat algorithm") + 
  tm_facets(nrow=3, ncol=3)


###############################
## Run the SaTScan algorithm ##
###############################
source("./SaTScan_auxFunctions.R")

# NOTE: To use SaTScan, you need to install not only the 'rsatscan' R package 
#       but also the SaTScan software itself (https://www.satscan.org/download.html)

satScan_output <- run_satScan(data=data.NAV |> mutate(year=year+2000), # To ensure compatibility with SaTScan,
                              carto=carto.NAV,
                              sslocation="C:/Program Files/SaTScan/")
print(satScan_output)

clusters_satScan <- get_clusterIDs_satScan(satScan_output, carto=carto.NAV, pVal.threshold=0.05)
str(clusters_satScan)


## Plot of the clusters identified by the SaTScan algorithm 
cluster.ID <- Reduce("+",lapply(1:ncol(clusters_satScan), function(x) x*as.numeric(as.character(clusters_satScan[,x]))))
cluster.ID <- as.data.frame(matrix(cluster.ID,nAreas, nYears, byrow=F))
colnames(cluster.ID) <- paste("Year", 1:nYears, sep=".")

carto.aux <- cbind(carto.NAV,cluster.ID)

tm_shape(carto.aux) +
  tm_polygons(fill=paste("Year", 1:nYears, sep="."),
              fill.scale = tm_scale_categorical(labels=c("No cluster",paste("Cluster",1:ncol(clusters_satScan))),
                                                values=c("lightgray", brewer.pal(ncol(clusters_satScan), "Set1"))),
              fill.legend=tm_legend("", show=TRUE, position=tm_pos_out("right","center")),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="SaTScan algorithm") + 
  tm_facets(nrow=3, ncol=3)


#############################################################
## Fit (global) spatio-temporal models using bigDM package ##
#############################################################
help("bigDM")

## No-cluster model
model_noCluster <- STCAR_INLA(carto=carto.NAV, data=data.NAV,
                              ID.area="ID", ID.year="year", O="obs", E="exp",
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="global", inla.mode="compact")

## Model including the clustering structure obtained through the GscanStat algorithm
model_GscanStat <- STCAR_INLA(carto=carto.NAV, data=data.NAV,
                              ID.area="ID", ID.year="year", O="obs", E="exp",
                              X=clusters_GscanStat,
                              spatial="BYM2", temporal="rw1", interaction="TypeIV",
                              model="global", inla.mode="compact")

model_GscanStat$summary.fixed


## Model including the clustering structure obtained through the SaTScan software
model_satScan <- STCAR_INLA(carto=carto.NAV, data=data.NAV,
                            ID.area="ID", ID.year="year", O="obs", E="exp",
                            X=clusters_satScan,
                            spatial="BYM2", temporal="rw1", interaction="TypeIV",
                            model="global", inla.mode="compact")

model_satScan$summary.fixed


######################################
## Compute model selection criteria ##
######################################
compute.MSC <- function(model){
  data.frame(deviance=model$dic$mean.deviance,
             pD=model$dic$p.eff,
             DIC=model$dic$dic,
             WAIC=model$waic$waic,
             LS=-sum(log(model$cpo$cpo)))
}

MODELS <- list('noCluster'=model_noCluster,
               'GscanStat'=model_GscanStat,
               'SaTScan'=model_satScan)

do.call(rbind,lapply(MODELS, compute.MSC))


##########################################################
## Maps of posterior median estimates of relative risks ##
##########################################################
tmap_mode("plot")

time.periods <- paste("Year", 1:nYears, sep=".")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.77,0.83,0.91,1,1.1,1.2,1.3,Inf)


## noCluster model
####################
model_noCluster.RR <- matrix(model_noCluster$summary.fitted.values$`0.5quant`,
                             nAreas, nYears, byrow=F)
colnames(model_noCluster.RR) <- time.periods

carto.noCluster <- cbind(carto.NAV, model_noCluster.RR)

Map.noCluster <- tm_shape(carto.noCluster) +
  tm_polygons(fill=time.periods,
              fill.scale=tm_scale(values=paleta, breaks=values, interval.closure="left"),
              fill.legend=tm_legend("Risk", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="noCluster model") + 
  tm_facets(nrow=3, ncol=3)

print(Map.noCluster)


## GscanStat model
####################
model_GscanStat.RR <- matrix(model_GscanStat$summary.fitted.values$`0.5quant`,
                             nAreas, nYears, byrow=F)
colnames(model_GscanStat.RR) <- time.periods

carto.GscanStat <- cbind(carto.NAV, model_GscanStat.RR)

Map.GscanStat <- tm_shape(carto.GscanStat) +
  tm_polygons(fill=time.periods,
              fill.scale=tm_scale(values=paleta, breaks=values, interval.closure="left"),
              fill.legend=tm_legend("Risk", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="GscanStat model") + 
  tm_facets(nrow=3, ncol=3)

print(Map.GscanStat)


## SaTScan model
##################
model_satScan.RR <- matrix(model_satScan$summary.fitted.values$`0.5quant`,
                           nAreas, nYears, byrow=F)
colnames(model_satScan.RR) <- time.periods

carto.satScan <- cbind(carto.NAV, model_satScan.RR)

Map.satScan <- tm_shape(carto.satScan) +
  tm_polygons(fill=time.periods,
              fill.scale=tm_scale(values=paleta, breaks=values, interval.closure="left"),
              fill.legend=tm_legend("Risk", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              fill.free=FALSE, col_alpha=0) + 
  tm_title(text="SatStan model") + 
  tm_facets(nrow=3, ncol=3)

print(Map.satScan)
