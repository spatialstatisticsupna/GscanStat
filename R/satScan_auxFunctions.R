## These functions are designed to use the SaTScan software within the simulation
## experiments described in the paper. They are not intended as general-purpose 
## code for using SaTScan in other contexts.

library(purrr)
library(sf)
library(tidyverse)
library(rsatscan)

###########################################################################################################
## run_SatScan:
## Generates the files and configurations needed to run SatScan clustering on the
## data and the given simulation number. The data object must contain the following variables:
##  muni: char. municipality code
##  year: int. year number of the simulated data, sequential values from 1 to the last year
##  sim: int. simulation number (optional)
##  obs: int. observed values
##  exp: num. expected values
##  pop: num. population for the municipality and year
###########################################################################################################

run_SatScan <- function(datos, carto, i.sim=NULL){
  carto <- st_as_sf(carto)
  # We transform the projection of the object to obtain lon/lat
  carto <- st_transform(carto, crs="+proj=longlat +datum=WGS84")
  carto$long <- st_coordinates(st_centroid(carto))[,1]
  carto$lat <- st_coordinates(st_centroid(carto))[,2]

  # If there are multiple simulations in the dataset, we select the given simulation
  if(!is.null(i.sim)){
    datos <- datos %>% filter(sim==i.sim) %>% mutate(year = 2000+ year)
  }


  datos.geo <- data.frame(muni=carto$muni, lat=carto$lat, long=carto$long)
  datos.cas <- datos %>% select(muni,obs,year)
  datos.pop <- datos %>% select(muni,year,pop)

  td = tempfile(pattern = "auxDir_satScan", tmpdir = "./")
  dir.create(td, recursive = TRUE)

  write.cas(datos.cas,td,"auxSaTScanFile")
  write.geo(datos.geo,td,"auxSaTScanFile")
  write.pop(datos.pop,td,"auxSaTScanFile")

  # write options
  # If the error "Error in ss.options(reset = TRUE) : object 'ssenv' not found" occurs, it means that
  # the environment variable 'ssenv' created when loading rsatscan has been deleted.
  # You need to detach the package and load it again.
  if(!exists("ssenv")){
    detach("package:rsatscan", unload = TRUE)
    library(rsatscan)
  }


  years <- sort(unique(datos.cas$year))

  invisible(ss.options(reset=TRUE))
  ss.options(list(CaseFile="auxSaTScanFile.cas",
                  StartDate=paste0(years[1],"/01/01"),
                  EndDate=paste0(years[length(years)],"/12/31"),
                  PopulationFile="auxSaTScanFile.pop",
                  CoordinatesFile="auxSaTScanFile.geo",
                  CoordinatesType=0, #long/lat coordinates
                  AnalysisType=3, # retorsprospective space-time
                  ScanAreas=3,  # hight and low risk
                  MaxSpatialSizeInPopulationAtRisk=50, # max spatial size 50% of population at risk (default and maximum allowed value)
                  MaxTemporalSizeInterpretation=0,MaxTemporalSize=90 # temporal size % of time periods (max allowed 90%)
                  ))
  ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))

  write.ss.prm(td,"auxSatScanFile")
  

  satScan_obj <- satscan(prmlocation=td, prmfilename="auxSatScanFile",
                         sslocation="/usr/local/bin/", ssbatchfilename="satscan",
                         verbose=FALSE)

  ## delete the temporary directory
  unlink(td, recursive = TRUE)

  return(satScan_obj)
}


#######################################################################################################
## get_clusteringPartition_satScan:
## Given the object returned by SatScan, creates a vector indicating which cluster each area belongs to
## for the spatio-temporal dataset (nAreas * nYears)
#######################################################################################################

get_clusteringPartition_satScan <- function(satScanOutput, pVal.threshold=0.05, nAreas=NULL, nYears=NULL){

  if(is.null(nAreas)){
    nAreas <- nrow(satScanOutput$rr)  
  }
  if(is.null(nYears)){
    # Extract dates
    start_date <- sub("StartDate=", "", satScan_output$prm[grepl("^StartDate=", satScan_output$prm)])
    end_date   <- sub("EndDate=", "",   satScan_output$prm[grepl("^EndDate=",   satScan_output$prm)])
    # Convert to Date
    start_date <- as.Date(start_date)
    end_date   <- as.Date(end_date)
    # Compute number of years
    nYears <- as.numeric(format(end_date, "%Y")) - as.numeric(format(start_date, "%Y")) + 1
  }
  
  # nClust <- length(cluster_duration$CLUSTER)
  nClust <- length(unique(satScanOutput$gis$CLUSTER))
  
  # cluster_duration <-
  #   satScanOutput$gis %>%
  #   filter(P_VALUE < pVal.threshold) %>%
  #   mutate(START_YEAR = as.integer(substr(START_DATE, 4, 4)),
  #          END_YEAR = as.integer(substr(END_DATE, 4, 4))) %>%
  #   select(CLUSTER, START_YEAR, END_YEAR, NUMBER_LOC)
  
  
  # cluster_duration <-
  #   satScanOutput$col %>%
  #   filter(P_VALUE < pVal.threshold) %>%
  #   mutate(START_YEAR = as.integer(substr(START_DATE, 4, 4)),
  #          END_YEAR = as.integer(substr(END_DATE, 4, 4))) %>%
  #   select(CLUSTER, START_YEAR, END_YEAR, NUMBER_LOC)
  
  # getting the time frame for each cluster
  time_frame_lines <- satScan_output$main[grep("^\\s*Time frame", satScan_output$main)]
  # getting the dates
  dates_only <- sub(".*: ", "", time_frame_lines)
  # getting the start and end date
  dates_split <- strsplit(dates_only, " to ")
  
  # Convertimos en data.frame
  time_frames <- do.call(rbind, dates_split)
  years_last_digit <- data.frame(
    START_YEAR = as.numeric(substr(time_frames[,1], 4, 4)),
    END_YEAR   = as.numeric(substr(time_frames[,2], 4, 4)),
    CLUSTER    = 1:nClust) 
    
  spatialCluster <-
    satScanOutput$gis %>% mutate(LOC_ID=as.numeric(as.character(LOC_ID))) %>% select(LOC_ID, CLUSTER)
  

  # spatialCluster$CLUSTER <- as.integer(spatialCluster$CLUSTER)
  # years_last_digit$CLUSTER <- as.integer(years_last_digit$CLUSTER)  

  # include cluster duration information
  cluster_duration <- merge(spatialCluster, years_last_digit, by = "CLUSTER", all.x = TRUE)

  clustering_partition_satScan <- numeric(nAreas*nYears)
  for(clust.i in 1:nClust){
    # cat("cluster ",clust.i,"/",nClust)
    index_clust.i <- spatialCluster %>% filter(CLUSTER==clust.i) %>% select(LOC_ID) %>% pull() %>% match(as.numeric(carto$muni))
    clust <- numeric(nAreas)
    clust[index_clust.i] <- clust.i

    # Start and end of the cluster clust.i
    years <- numeric(nYears)
    years[years_last_digit[clust.i,'START_YEAR']:years_last_digit[clust.i,'END_YEAR']] <- 1
    # Using the Kronecker product: if years==0, a vector of 0s with the dimension of clust is created;
    # if years==1, the clust vector is used, and all are concatenated
    aux_clusteringPartition <- kronecker(years,clust)
    clustering_partition_satScan <- clustering_partition_satScan + aux_clusteringPartition
  }
  return(clustering_partition_satScan)
}

#######################################################################################################
## get_clusteringIDMatrix_satScan:
## Given the object returned by SatScan, creates a matrix with as many columns as clusters,
## indicating in each column whether the area belongs to the cluster or not.
## The matrix is for the spatio-temporal dataset (nAreas * nYears)
#######################################################################################################
get_clusteringIDMatrix_satScan <- function(satScanOutput, pVal.threshold=0.05, nAreas=NULL, nYears=NULL){
  
  if(is.null(nAreas)){
    nAreas <- nrow(satScanOutput$rr)  
  }
  if(is.null(nYears)){
    # Extract dates
    start_date <- sub("StartDate=", "", satScan_output$prm[grepl("^StartDate=", satScan_output$prm)])
    end_date   <- sub("EndDate=", "",   satScan_output$prm[grepl("^EndDate=",   satScan_output$prm)])
    # Convert to Date
    start_date <- as.Date(start_date)
    end_date   <- as.Date(end_date)
    # Compute number of years
    nYears <- as.numeric(format(end_date, "%Y")) - as.numeric(format(start_date, "%Y")) + 1
  }
  
  # nClust <- length(cluster_duration$CLUSTER)
  nClust <- length(unique(satScanOutput$gis$CLUSTER))
  
  # cluster_duration <-
  #   satScanOutput$gis %>%
  #   filter(P_VALUE < pVal.threshold) %>%
  #   mutate(START_YEAR = as.integer(substr(START_DATE, 4, 4)),
  #          END_YEAR = as.integer(substr(END_DATE, 4, 4))) %>%
  #   select(CLUSTER, START_YEAR, END_YEAR, NUMBER_LOC)
  
  # getting the time frame for each cluster
  time_frame_lines <- satScan_output$main[grep("^\\s*Time frame", satScan_output$main)]
  # getting the dates
  dates_only <- sub(".*: ", "", time_frame_lines)
  # getting the start and end date
  dates_split <- strsplit(dates_only, " to ")
  
  # Convertimos en data.frame
  time_frames <- do.call(rbind, dates_split)
  years_last_digit <- data.frame(
    START_YEAR = as.numeric(substr(time_frames[,1], 4, 4)),
    END_YEAR   = as.numeric(substr(time_frames[,2], 4, 4)),
    CLUSTER    = 1:nClust) 
  
  spatialCluster <-
    satScanOutput$gis %>% mutate(LOC_ID=as.numeric(as.character(LOC_ID))) %>% select(LOC_ID, CLUSTER)
  
  
  # spatialCluster$CLUSTER <- as.integer(spatialCluster$CLUSTER)
  # years_last_digit$CLUSTER <- as.integer(years_last_digit$CLUSTER)  
  
  # include cluster duration information
  cluster_duration <- merge(spatialCluster, years_last_digit, by = "CLUSTER", all.x = TRUE)
  
  clustering_IDMatrix_satScan <- NULL
  for(clust.i in 1:nClust){
    # cat("cluster ",clust.i,"/",nClust)
    index_clust.i <- spatialCluster %>% filter(CLUSTER==clust.i) %>% select(LOC_ID) %>% pull() %>% match(as.numeric(carto$muni))
    clust <- numeric(nAreas)
    clust[index_clust.i] <- clust.i
    
    # Start and end of the cluster clust.i
    years <- numeric(nYears)
    years[years_last_digit[clust.i,'START_YEAR']:years_last_digit[clust.i,'END_YEAR']] <- 1
    # Using the Kronecker product: if years==0, a vector of 0s with the dimension of clust is created;
    # if years==1, the clust vector is used, and all are concatenated
    aux_clusteringPartition <- kronecker(years,clust)
    clustering_IDMatrix_satScan <- cbind(clustering_IDMatrix_satScan, aux_clusteringPartition)
  }
  colnames(clustering_IDMatrix_satScan) <- paste0("C",1:nClust)
  return(clustering_IDMatrix_satScan)
}

