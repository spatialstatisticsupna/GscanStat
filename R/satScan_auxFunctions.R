library(purrr)
library(sf)
library(tidyverse)
library(rsatscan)

#' Runs SaTScan clustering on a spatio-temporal dataset.
#'
#' This function prepares the necessary files and configurations to execute a
#' SaTScan analysis. It takes a data frame containing spatio-temporal
#' information and a cartography object, then generates the required input
#' files for the SatScan software. IMPORTANT: This function is designed to use 
#' the SaTScan software within the simulation experiments described in the paper
#' "Improving Disease Risk Estimation in Small Areas by Accounting for 
#' Spatio-Temporal Local Discontinuities". They are not intended as 
#' general-purpose code for using SaTScan in other contexts.
#'
#' @param datos A data frame containing the spatio-temporal data. It must
#'   include the following columns:
#'   \itemize{
#'     \item \strong{muni}: A character vector with municipality codes.
#'     \item \strong{year}: An integer vector representing the year.
#'     \item \strong{sim}: An integer vector indicating the simulation number
#'       (optional).
#'     \item \strong{obs}: An integer vector of the observed case counts.
#'     \item \strong{exp}: A numeric vector of the expected case counts.
#'     \item \strong{pop}: A numeric vector of the population for each
#'       municipality and year.
#'   }
#' @param carto An object containing the cartographic information (e.g.,
#'   municipality coordinates and polygons) required by SaTScan.
#' @param i.sim An optional integer specifying the simulation number to be
#'   processed. Defaults to `NULL`.
#'
#' @return An objecto with the SaTScan output as provided by function `satscan`
#' from the `rsatscan` package.
#'
#' @details The main purpose of this function is to automate the preparation of
#'   input files, facilitating the process of running a SaTScan analysis in
#'   experiments with simulated data, as described in the paper 
#'   "Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal 
#'   Local Discontinuities". **WARNING:** This function cannot be used with real data, 
#'   since internal age and sex standardization would require providing 
#'   `age` and `sex` covariates.
#'
#' @examples
#' # Assuming 'my_data' and 'my_carto' are prepared
#' # run_SatScan(datos = my_data, carto = my_carto)
run_SatScan <- function(datos, carto, i.sim = NULL) {
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


#' Creates a data.frame object indicating cluster membership from a SatScan output.
#'
#' Given the output object from a SaTScan analysis, this function generates a
#' data.frame object. The data.frame has as many factors (columns) as significant clusters 
#' and indicates for each area (rows) whether it belongs to a specific cluster (1) or not (0).
#' The output data.frame is designed for spatio-temporal datasets, with dimensions
#' corresponding to `no. of Areas * no. of Time periods`.
#'
#' @param satScanOutput An object returned by a SaTScan analysis. This object
#'   is expected to contain information about detected clusters, including their
#'   locations, sizes, and p-values.
#' @param pVal.threshold A numeric value specifying the significance threshold for
#'   the p-value. Only clusters with a p-value less than this threshold will be 
#'   included in the output data.frame. The default is 0.05.
#'
#' @return A data.frame with a number of columns equal
#'   to the number of significant clusters (those meeting the `pVal.threshold`).
#'   For each area i at time period t, a value of 1 indicates membership in the cluster, 
#'   and 0 indicates non-membership.
#'
#' @examples
#' # Assuming 'satScanResult' is a valid object from a SatScan run
#' # cluster_matrix <- get_clusterIDs_satScan(
#' #   satScanResult,
#' #   pVal.threshold = 0.05)
get_clusterIDs_satScan <- function(satScanOutput, pVal.threshold=0.05){
  
  # compute the number of areas (spatial dimension)
    nAreas <- nrow(satScanOutput$rr)  
    # compute the number of years (temporal dimension)
    # Extract dates
    start_date <- sub("StartDate=", "", satScan_output$prm[grepl("^StartDate=", satScan_output$prm)])
    end_date   <- sub("EndDate=", "",   satScan_output$prm[grepl("^EndDate=",   satScan_output$prm)])
    # Convert to Date
    start_date <- as.Date(start_date)
    end_date   <- as.Date(end_date)
    # Compute number of years
    nYears <- as.numeric(format(end_date, "%Y")) - as.numeric(format(start_date, "%Y")) + 1
  
    nClust <- length(unique(satScanOutput$gis$CLUSTER))
  
 
  # getting the time frame for each cluster
  time_frame_lines <- satScan_output$main[grep("^\\s*Time frame", satScan_output$main)]
  # getting the dates
  dates_only <- sub(".*: ", "", time_frame_lines)
  # getting the start and end date
  dates_split <- strsplit(dates_only, " to ")
  
  # getting start and end year
  time_frames <- do.call(rbind, dates_split)
  years_last_digit <- data.frame(
    START_YEAR = as.numeric(substr(time_frames[,1], 4, 4)),
    END_YEAR   = as.numeric(substr(time_frames[,2], 4, 4)),
    CLUSTER    = 1:nClust) 
  
  spatialCluster <-
    satScanOutput$gis %>% mutate(LOC_ID=as.numeric(as.character(LOC_ID))) %>% select(LOC_ID, CLUSTER)

  # include cluster duration information
  cluster_duration <- merge(spatialCluster, years_last_digit, by = "CLUSTER", all.x = TRUE)
  
  clustering_IDMatrix_satScan <- NULL
  for(clust.i in 1:nClust){
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
# coercing to data.frame
    clusters <- as.data.frame(clustering_IDMatrix_satScan) %>%
    dplyr::mutate(across(everything(), as.factor))
  
  return(clusters)
}

