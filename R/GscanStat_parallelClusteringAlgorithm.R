library(INLA)
# library(maptools)
library(spdep)
library(sp)
library(RColorBrewer)
library(parallel)
library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
library(sf)
library(tmap)



##' Extend the spatial adajacency matrix to spatio-temporal adjacency matrix
##' @param W_ an adjacency n x n 'matrix' with the spatial neighborhood
##' @param nblocks integer, number of considered time periods (number of blocks in the spatio-temporal neighoborhood matrix)
##' @return a sparse (nblocks*k x nblocks*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
timeExtended_W <- function(W_, nblocks) {

        # Ait, Ait+1 and Ait-1 are neighbors
        n <- nrow(W_)
        blocks <- diag(nblocks)
        diag(blocks) <- 0
        diag(blocks[-1,]) <- 1
        diag(blocks[,-1]) <- 1

        result1 <- kronecker(blocks, diag(n))
        result2 <- kronecker(diag(nblocks),W_)
        return(result1+result2)
}

###################################################################################################

###################################################################################################

# Calculates the maximum spatio-temporal neighborhood window, considering the maxAreas closest
# spatial neighbors over a time period from t-timeWindowSize to t+timeWindowSize,
# where t is the time associated with each area

# i : number of the area where the window is centered
# t : time instant where the window is centered
spatioTemporalMaxKnn_Ait <- function(i,t,geoDist,maxAreas, timeWindowSize, n.areas=nrow(geoDist), n.years){

        zones_spatial <- order(geoDist[i,])[1:maxAreas]

        # vector of indices with the spatial neighborhood at time t
        spatialNeighborhood_t <- zones_spatial + (t-1)*n.areas

        # add the previous time instants within the timeWindowSize
        prev.t <- NULL
        if(t>1){
                prev.t <- seq(t-1,max(1,t-timeWindowSize))
        }
        post.t <- NULL
        if(t<n.years){
                post.t <- seq(t+1,min(n.years,t+timeWindowSize))
        }
        other.t <- c(prev.t,post.t)

        remaining.areas <- as.vector(sapply(other.t, FUN=function(t.aux,zones_spatial,n.areas){zones_spatial + (t.aux-1)*n.areas}, zones_spatial, n.areas))
        neighborWindow_spatioTemporal <- c(spatialNeighborhood_t, remaining.areas)
        return(neighborWindow_spatioTemporal)
}


##############################################
# Calculation of lambda according to the description in Tango's book, but on a logarithmic scale
# (suitable for low-risk clusters).
# It is also valid for the restricted version of the statistic, since the areas with mid-pvalue
# have already been removed previously during the search.

calculateLogLambda.spatioTemporal.poisson <- function(Zwindow,data,lowRiskZone=FALSE,n=sum(data$obs)){

        nZ <- sum(data[Zwindow,'obs'])
        eZ <- sum(data[Zwindow,'exp'])

        ratio.Z <- nZ/eZ
        ratio.outZ <- (n-nZ)/(n-eZ)
        if((!lowRiskZone && (ratio.Z<ratio.outZ) )|| (lowRiskZone && (ratio.Z>ratio.outZ))){
                return(.Machine$double.xmin)
        }
        else{
                # To avoid numerical calculation issues with logs, we add a very small value: 1e-10
                return(nZ*log(ratio.Z + 1e-10) + (n-nZ)*log(ratio.outZ + 1e-10))
        }
}

##############################################
calculateLogLambda.spatioTemporal.poisson.multipleCandidateAreas <- function(Zwindow,candidateZs,data,lowRiskZone=FALSE,n=sum(data$obs)){


        nZ <- sum(data$obs[Zwindow]) +  data$obs[candidateZs]
        eZ <- sum(data$exp[Zwindow]) + data$exp[candidateZs]

        ratio.Z <- nZ/eZ
        ratio.outZ <- (n-nZ)/(n-eZ)

        lambda <- nZ*log(ratio.Z + 1e-10) + (n-nZ)*log(ratio.outZ + 1e-10)
        condition <- ((!lowRiskZone) * (ratio.Z<ratio.outZ)) + (lowRiskZone * (ratio.Z>ratio.outZ))

        result <- condition*.Machine$double.xmin + (1-condition)*lambda
        return(result)
}

#####################################################################################################
## Given a vector k_nearest with the nearest neighbors (in space and time) to an area Ai
## (the first in the vector), obtains the subset of k_nearest starting from Ai that is connected and
## yields the best value of the lambda statistic (performing a local search, that is, adding
## one by one the best neighbor until the maximum lambda is reached)

greedySearch_bestConnectedSpatioTemporalNeighborhood <- function(k_nearest,
                                                                 adjacency_matrix,
                                                                 data,
                                                                 n = sum(data$obs)){
        ## We start the search in area Ai at time t (this is given by the first value of k_nearest,
        ## which corresponds to the index of the data where Ai is located at time t)
        Ai <- k_nearest[1]

        lowRiskCluster <- data$lowRisk[Ai]

        Z.final<- c(Ai)
        areas   <- k_nearest[-1]


        lambda <- calculateLogLambda.spatioTemporal.poisson(Z.final,data,lowRiskZone=lowRiskCluster)
        stop =FALSE

        while(!stop){
                # Regions of areas connected to some area in Z.final
                index <- rowSums(adjacency_matrix[areas,Z.final, drop=FALSE])>0

                # neighbours contains all possible neighbors of the areas in Z.final
                neighbours <- areas[index]

                if(length(neighbours) > 0){
                        ## We calculate all the lambdas for each possible neighbor
                        lambda.neighbours <- calculateLogLambda.spatioTemporal.poisson.multipleCandidateAreas(Zwindow=Z.final,
                                                                                                              candidateZs=neighbours,
                                                                                                              data=data,
                                                                                                              lowRiskZone=lowRiskCluster,
                                                                                                              n=n)

                        bestArea.i   <- which.max(lambda.neighbours)
                        bestArea     <- neighbours[bestArea.i]
                        bestLambda <- lambda.neighbours[bestArea.i]
                        # cat("Area ", Ai, "neighboursLength = ", length(neighbours), " bestLambda = ", bestLambda, "  currentLambda = ", lambda, "\n")
                        if(bestLambda > lambda ){
                                lambda <- bestLambda
                                Z.final <- c(Z.final,bestArea)
                                areas <- areas[areas!=bestArea]
                        }
                        else{
                                stop <- TRUE
                        }
                }
                else{
                        stop <- TRUE
                }
        }

        result <- tibble(areas=list(Z.final), lambda=lambda, lowRiskCluster=lowRiskCluster)
        return(result)
}

# # ###################################################################################################
# # # ###################################################################################################
getClusters.SpatioTemporalScanStat.greedy_Parallel <- function(data,
                                                               geoDist,
                                                               maxAreas,
                                                               timeWindowSize,
                                                               W_spatioTemp,
                                                               verbose=TRUE,
                                                               maxCPUs = detectCores()-1,
                                                               areas.to.process=ceiling(nrow(data)/(maxCPUs)),
                                                               computationCluster=NULL){
        n.obs <- sum(data$obs)
        n.areas <- data %>% filter(year==data$year[1]) %>% nrow(.)
        n.years <- data %>% dplyr::select(year) %>% unique(.) %>% nrow(.)

        # if cluster is not created and given as parameter, create a cluster for parallel computing
        if(is.null(computationCluster)){
                #create a cluster for parallel computing in FORK mode - only for linux systems!!
                cl <- makeCluster(maxCPUs, type = "FORK")
                registerDoParallel(cl)
        } else{
                cl <- computationCluster
        }
        if(verbose){
                cat("Running in parallel using a maximun of",maxCPUs, " CPUs, and processing", areas.to.process, "areas per CPU\n\n")
                ## create a temp folder to save logs from each process
                tmpFolder.name <- paste0("TMP_",Sys.info()[["nodename"]],"_PID",Sys.getpid())
                dir.create(tmpFolder.name)
                cat("log files for every process are saved in", tmpFolder.name, "temp folder\n")
                }


        results <- foreach(Aij = seq(1,nrow(data),areas.to.process)) %dopar% {
        # results <- foreach(Aij = seq(1,10000,areas.to.process)) %dopar% {
                result <- NULL

                if(verbose){
                        sink(paste0("./",tmpFolder.name,"/logCPU_",Sys.getpid(),".log"))
                        cat("Processing areas from", Aij,"to", min(Aij+areas.to.process-1,n.areas*n.years),"\n\n")
                }

                tiempo <- system.time({
                        for(Aij_local in Aij:min(Aij+areas.to.process-1,n.areas*n.years)){
                                if(verbose){cat("Area",Aij_local,"\n")}
                                i <- Aij_local %% (n.areas)
                                if(i==0){i<-n.areas}
                                t <- min(floor(Aij_local/n.areas)+1,n.years)

                                knn_maxNeighbourhood <- spatioTemporalMaxKnn_Ait(i,t,geoDist,maxAreas, timeWindowSize, n.areas=n.areas, n.years=n.years)
                                aux <- greedySearch_bestConnectedSpatioTemporalNeighborhood(k_nearest=knn_maxNeighbourhood, adjacency_matrix=W_spatioTemp, data=data, n=n.obs)
                                aux$index <- Aij_local
                                aux$i <- i
                                aux$t <- t
                                result <- rbind(result,aux)
                        }
                })
                result$time <- tiempo['elapsed']

                if(verbose){
                        cat("\nElapsed time:", tiempo['elapsed']/60,"min \n")
                        sink()
                }

                return(result)
        }

        # if the computation cluster has been created in this method, we stop the cluster
        if(is.null(computationCluster)){
                stopCluster(cl)
        }

        if(verbose){
                cat("All the areas have been processed looking for clusters\n\n")
                unlink(tmpFolder.name, recursive = TRUE)
        }

        zones <- do.call(rbind,results) %>% arrange(index)
        return(zones)
}

# # ###################################################################################################

get.nonOverlapping.clusters <- function(Zs, pVal_threshold = 0.05){
        Zs_filtrado <- Zs %>% arrange(pVal,desc(lambda)) %>% slice(1:max(10, Zs %>% filter(pVal < pVal_threshold) %>% nrow(.)))

        nonOverlap <- logical(nrow(Zs_filtrado))
        nonOverlap[1] <- TRUE
        aux <- Zs_filtrado$areas[[1]]
        for(i in 2:length(nonOverlap)){
                cluster_st_i <- Zs_filtrado$areas[[i]]
                is.overlap <- cluster_st_i %in% aux
                overlap <- any(is.overlap)

                if(!overlap){
                        nonOverlap[i] <- TRUE
                        aux <- c(aux, cluster_st_i)
                }
        }

        return(Zs_filtrado[nonOverlap,])
}


# # ###################################################################################################
# We use this to include a column with the cluster identification in the dataset.
# Returns the dataset with a new column.
# Due to model fitting issues, in cases where there are low-risk clusters with all areas having 0 cases,
# these clusters are removed.

getClusterIdentityMatrix <- function(data,Zs, max_pVal=0.05){

        Zs <- Zs %>% filter(pVal < max_pVal)

        n.data  <- nrow(data)
        n.clust <- nrow(Zs)

        if(n.clust==0){
                warning("There is no significant spatio-temporal cluster in the provided data")
                return(NULL)
        }

        clusterIdMatrix <- matrix(0,nrow = n.data, ncol = n.clust)

        for(iClust in 1:n.clust){
                clust <- Zs$areas[iClust][[1]]
                clusterIdMatrix[clust,iClust] <- 1
        }

        # We check if there is any cluster with all 0s in the SMR. When we have a
        # cluster with 0s, model fitting causes problems. Therefore,
        # we will not consider such clusters.
        structural0_cluster <- apply(clusterIdMatrix, MARGIN=2, FUN=function(ids, data){all(data$obs[as.logical(ids)]==0)}, data=data)

        if(any(structural0_cluster)){
          warning(paste0("Cluster(s) ", which(structural0_cluster), " contain(s) only areas with SMR=0. The cluster(s) is/are removed  !!!\n"))
          clusterIdMatrix <- clusterIdMatrix[,!structural0_cluster, drop = FALSE]
        }
        return(clusterIdMatrix)
}
# # ###################################################################################################
# Monte Carlo test to obtain the distribution of max logLambda for both lowRisk and highRisk
# in parallel computing using the maximum number of CPUs when performing the greedy search over the areas
mcTest.logLambda.poisson_Parallel <- function(data,
                                              W_spatioTemp,
                                              geoDist,
                                              maxAreas,
                                              timeWindowSize,
                                              rep=999,
                                              verbose.steps=1,
                                              maxCPUs = detectCores()-1,
                                              computationCluster=NULL){

        n.obs <- sum(data$exp)
        n <- nrow(data)

        lambdaDist <- tibble(lowRisks =numeric(), highRisks=numeric())

        tiempo <- system.time({
                for(i in 1:rep) {
                        set.seed(i)
                        data$obs <- rmultinom(1,n.obs,data$exp)
                        # The value of whether it is a lowRisk cluster or not can vary when generating the random data
                        data$lowRisk <- data$obs < data$exp
                        if(i%%verbose.steps == 0){
                                cat("MonteCarlo - it",i,": dataset created. Looking for clusters in dataset ",i,"...\n")
                        }
                        t.mcIt <- system.time({
                                Zs <- getClusters.SpatioTemporalScanStat.greedy_Parallel(data=data,
                                                                                         geoDist = geoDist,
                                                                                         maxAreas = maxAreas,
                                                                                         timeWindowSize = timeWindowSize,
                                                                                         W_spatioTemp = W_spatioTemp,
                                                                                         verbose = (i%%verbose.steps == 0),
                                                                                         maxCPUs = maxCPUs,
                                                                                         computationCluster = computationCluster)
                        })
                        if(i%%verbose.steps == 0){
                                cat(paste0("it",i," cluster detection finished (time:",t.mcIt['elapsed']/60,"min\n\n"))
                        }

                        lambdaLowRisk  <- Zs %>% filter(lowRiskCluster) %>% dplyr::select(lambda) %>% max()
                        lambdaHighRisk <- Zs %>% filter(!lowRiskCluster) %>% dplyr::select(lambda) %>% max()
                        if(i%%verbose.steps == 0){
                                cat(paste0("it",i,": lambdaHighRisk=", lambdaHighRisk, ";  lambdaLowRisk=", lambdaLowRisk,"\n\n"))
                        }
                        lambdaDist <- rbind(lambdaDist, tibble_row(lowRisks=lambdaLowRisk, highRisks=lambdaHighRisk))


                }
        })
        cat("MonteCarlo total processing time: ", tiempo['elapsed']/60," min\n")
        return(lambdaDist)
}


#####################################################################################################
# Runs Monte Carlo test iterations only between first.MC.It and last.MC.It and stores the
# max logLambda for both lowRisk and highRisk for each iteration. Later,
# all results can be combined to construct the logLambda distribution for lowRisk and highRisk.
# In parallel computing using the maximum number of CPUs when performing the greedy search over the areas
mcTest.logLambda.poisson.PartialIts_Parallel <- function(data,
                                                         W_spatioTemp=NULL,
                                                         geoDist,
                                                         maxAreas,
                                                         timeWindowSize,
                                                         first.MC.It,
                                                         last.MC.It,
                                                         restricted=NULL,
                                                         verbose.steps=1,
                                                         maxCPUs = detectCores()-1,
                                                         computationCluster=NULL,
                                                         file.name.out){


        n.areas <- data %>% filter(year==data$year[1]) %>% nrow(.)
        n.years <- data %>% dplyr::select(year) %>% unique(.) %>% nrow(.)

        if(is.null(W_spatioTemp)){
                W_spatioTemp.filename <- paste0("W_spatioTem_Map_nAreas",n.areas,"_nYears",n.years,".RData")
                if(file.exists(W_spatioTemp.filename)){
                        warning(paste0("Loading spatio-temporal neighborhood matrix from file ",W_spatioTemp.filename))
                        load(W_spatioTemp.filename)
                } else {
                        ## compute the spatio-temporal neighborhood matrix
                        if(verbose){cat("computing spatio-temporal neighborhood\n")}
                        W_spatioTemp <- timeExtended_W(W, n.years)
                        #save(W_spatioTemp,file=W_spatioTemp)
                }
        }

        ## If we use the restricted statistic, we remove from the neighbors matrix those
        ## areas that have a mid.pval >= the given value
        if(!is.null(restricted)){
                ## Using the restricted log-likelihood = disconnecting areas with high two-tail mid-pvalue
                cat("RESTRICTED STATISTIC: removing low-likely areas from neighborhoods\n")
                data <- data %>% mutate(mid.pval = 2*pmin(ppois(obs+1,exp,lower.tail = FALSE)+0.5*dpois(obs+1,exp),ppois(obs-1,exp,lower.tail = TRUE)+0.5*dpois(obs-1,exp)))
                areas.to.disconnect <- (data$mid.pval >= restricted)
                W_spatioTemp[areas.to.disconnect, ] <- 0
                W_spatioTemp[ , areas.to.disconnect] <- 0
        }

        n.obs <- sum(data$exp)
        n <- nrow(data)
        lambdaDist <- tibble(lowRisks =numeric(), highRisks=numeric(),MC.it.time.mins=numeric())

        tiempo <- system.time({
                for(i in first.MC.It:last.MC.It) {
                        set.seed(i)
                        data$obs <- rmultinom(1,n.obs,data$exp)
                        # The value of whether it is a lowRisk cluster or not can vary when generating random data
                        data$lowRisk <- data$obs < data$exp
                        if(i%%verbose.steps == 0){
                                cat("MonteCarlo - it",i,": dataset created. Looking for clusters in dataset ",i,"...\n")
                        }
                        t.mcIt <- system.time({
                                Zs <- getClusters.SpatioTemporalScanStat.greedy_Parallel(data=data,
                                                                                         geoDist = geoDist,
                                                                                         maxAreas = maxAreas,
                                                                                         timeWindowSize = timeWindowSize,
                                                                                         W_spatioTemp = W_spatioTemp,
                                                                                         verbose = (i%%verbose.steps == 0),
                                                                                         computationCluster = computationCluster)
                        })
                        if(i%%verbose.steps == 0){
                                cat(paste0("it",i," cluster detection finished (time:",t.mcIt['elapsed']/60,"min)\n\n"))
                        }

                        lambdaLowRisk  <- Zs %>% filter(lowRiskCluster) %>% dplyr::select(lambda) %>% max()
                        lambdaHighRisk <- Zs %>% filter(!lowRiskCluster) %>% dplyr::select(lambda) %>% max()
                        if(i%%verbose.steps == 0){
                                cat(paste0("it",i,": lambdaHighRisk=", lambdaHighRisk, ";  lambdaLowRisk=", lambdaLowRisk,"\n\n"))
                        }
                        lambdas <- tibble_row(lowRisks=lambdaLowRisk, highRisks=lambdaHighRisk, MC.it.time.mins = t.mcIt['elapsed']/60)
                        lambdaDist <- rbind(lambdaDist, lambdas)
                        save(lambdas,file=paste0(file.name.out,"_it",i,".RData"))


                }
        })
        cat("MonteCarlo partial processing time (it",first.MC.It,":",last.MC.It,")  : ", tiempo['elapsed']/60," min\n")
        lambdaDist$average.MC.it.time.mins <- tiempo['elapsed']/60 /nrow(lambdaDist)
        save(lambdaDist,file=paste0(file.name.out,"_it_",first.MC.It,"-",last.MC.It,".RData"))
}


# # ###################################################################################################
# Calculates the maximum likelihood windows for each studied area, returning the results in a list
# (Zs) saved within the given file name.
#
# restricted controls whether the restricted likelihood ratio (NULL) is used or not
# (giving a value between 0 and 1). If a value is given, it is used to filter the p-value
# (probability that a Poisson distribution with mean equal to the expected cases
# generates a value as observed or more extreme).
#
# computation cluster can be used when parallel=TRUE to pass an already created cluster
# and reuse it instead of creating and registering a cluster each time. This should save time.
getSpatioTemporalClustersAij_parallel <- function(data, geoDist, W,
                                                  maxAreas=15,
                                                  timeWindowSize=2,
                                                  W_spatioTemp=NULL,
                                                  restricted = NULL,
                                                  computationCluster = NULL,
                                                  verbose=FALSE,
                                                  file.out){

        ## We modify the dataset to include information on whether we have a lowRisk cluster
        ## and about the mid.pval
        data <- data %>% mutate(lowRisk = (obs < exp))


        n.areas <- data %>% filter(year==data$year[1]) %>% nrow(.)
        n.years <- data %>% dplyr::select(year) %>% unique(.) %>% nrow(.)

        if(is.null(W_spatioTemp)){
                W_spatioTemp.filename <- paste0("W_spatioTem_Map_nAreas",n.areas,"_nYears",n.years,".RData")
                if(file.exists(W_spatioTemp.filename)){
                        warning(paste0("Loading spatio-temporal neighborhood matrix from file ",W_spatioTemp.filename))
                        load(W_spatioTemp.filename)
                } else {
                        ## compute the spatio-temporal neighborhood matrix
                        if(verbose){cat("computing spatio-temporal neighborhood\n")}
                        W_spatioTemp <- timeExtended_W(W, n.years)
                        #save(W_spatioTemp,file=W_spatioTemp)
                }
        }


        ## If we use the restricted statistic, we remove from the neighbors matrix those
        ## areas that have a mid.pval >= the given value
        if(!is.null(restricted)){
                ## Using the restricted log-likelihood = disconnecting areas with high mid.pval
                if(verbose){cat("RESTRICTED STATISTIC: removing low-likely areas from neighborhoods\n")}
                data <- data %>% mutate(mid.pval = 2*pmin(ppois(obs+1,exp,lower.tail = FALSE)+0.5*dpois(obs+1,exp),ppois(obs-1,exp,lower.tail = TRUE)+0.5*dpois(obs-1,exp)))
                areas.to.disconnect <- (data$mid.pval >= restricted)
                W_spatioTemp[areas.to.disconnect, ] <- 0
                W_spatioTemp[ , areas.to.disconnect] <- 0
        }

        cat("Evaluating windows with greedy-scanStatistic \n")
        tiempo <- system.time({
                Zs <- getClusters.SpatioTemporalScanStat.greedy_Parallel(data=data,
                                                                         geoDist=geoDist,
                                                                         maxAreas=maxAreas,
                                                                         timeWindowSize = timeWindowSize,
                                                                         W_spatioTemp=W_spatioTemp,
                                                                         verbose=verbose,
                                                                         maxCPUs = maxCPUs,
                                                                         computationCluster = computationCluster)
        })

        cat("Time evaluating areas: ", tiempo['elapsed']/60, " min\n")
        Zs$global.time.mins <- tiempo['elapsed']/60
        save(Zs, file = paste(file.out,".RData"))
        return(Zs)
}

#####################################################################################################
# Reads data from files for the clusters identified for each area and the lambda values
# from the Monte Carlo test iterations, and generates the list of non-overlapping clusters and their p-value.
# The files are obtained using the method getSpatioTemporalClustersAij_parallel for the area clusters
# and mcTest.logLambda.poisson.PartialIts_Parallel for the lambda values from the MC test.
#
# clusterAij : string with the file name containing the Zs object with the ML cluster for each area
# lambdasMC : vector of strings with the file names of the Monte Carlo iterations to consider

getSignificantSpatioTemporalClusters_fromFiles <- function(clusterAij, lambdasMC){

        cat("Loading ML clusters for areas Aij from file\n")
        load(clusterAij)

        lambdaDist <- NULL
        reps.mc <- length(lambdasMC)

        lambdasMC <- paste0("SpainExp_cluster_spN2000tmpN5_MC_it", 1:99,".RData")

        cat("Obtaining lambda distribution from files\n")
        for(mcItsName.i in lambdasMC){
                load(mcItsName.i)
                lambdaDist <- rbind(lambdaDist,lambdas)
        }



        cat("Calculating p-values\n")
        pVal <- apply(Zs, MARGIN=1, FUN = function(x,lambdaDist){
                if(x$lowRiskCluster){
                        lambdaDist %>% dplyr::select(lowRisks) %>% filter(lowRisks > x$lambda) %>% nrow() /(reps.mc+1)
                }
                else{
                        lambdaDist %>% dplyr::select(highRisks) %>% filter(highRisks > x$lambda) %>% nrow() /(reps.mc+1)
                }
        }, lambdaDist=lambdaDist)


        Zs$pVal <- pVal

        # We combine the high-risk and low-risk clusters and sort them first by p-value and then by size.
        # For clusters with equal significance, priority is given to the larger ones.
        Zs <- Zs %>% mutate(clusterSize=sapply(areas, FUN=function(x){length(x)})) #%>% arrange(pVal, desc(clusterSize))

        return(Zs)
}

#####################################################################################################
# restricted controls whether the restricted likelihood ratio (NULL) is used or not
# (giving a value between 0 and 1). If a value is provided, it is used to filter the p-value
# (probability that a Poisson distribution with mean equal to the expected number of cases
# generates a value as observed or more extreme).
#
# computation cluster can be used when parallel=TRUE to pass an already created cluster
# and reuse it instead of creating and registering a cluster each time. This should save time.
getSignificantSpatioTemporalClusters_parallel <- function(data, geoDist, W,
                                                          maxAreas=15, timeWindowSize=2, reps.mc=999,
                                                          restricted = NULL,
                                                          parallel=FALSE,
                                                          computationCluster = NULL,
                                                          verbose=FALSE,
                                                          maxCPUs=detectCores()-1){

        ## We modify the dataset to include information on whether we have a lowRisk cluster
        ## and about the mid.pval
        data <- data %>% mutate(lowRisk = (obs < exp))


        n.areas <- data %>% filter(year==data$year[1]) %>% nrow(.)
        n.years <- data %>% dplyr::select(year) %>% unique(.) %>% nrow(.)

        W_spatioTemp.filename <- paste0("W_spatioTem_Map_nAreas",n.areas,"_nYears",n.years,".RData")
        if(file.exists(W_spatioTemp.filename)){
                warning(paste0("Loading spatio-temporal neighborhood matrix from file ",W_spatioTemp.filename))
                load(W_spatioTemp.filename)
        } else {
                ## compute the spatio-temporal neighborhood matrix
                if(verbose){cat("computing spatio-temporal neighborhood\n")}
                W_spatioTemp <- timeExtended_W(W, n.years)
                #save(W_spatioTemp,file=W_spatioTemp)
        }


        ## If we use the restricted statistic, we remove from the neighbors matrix those
        ## areas that have a mid.pval >= the given value
        if(!is.null(restricted)){
                ## Using the restricted log-likelihood = disconnecting areas with high mid.pval
                if(verbose){cat("RESTRICTED STATISTIC: removing low-likely areas from neighborhoods\n")}
                data <- data %>% mutate(mid.pval = 2*pmin(ppois(obs+1,exp,lower.tail = FALSE)+0.5*dpois(obs+1,exp),ppois(obs-1,exp,lower.tail = TRUE)+0.5*dpois(obs-1,exp)))
                areas.to.disconnect <- (data$mid.pval >= restricted)
                W_spatioTemp[areas.to.disconnect, ] <- 0
                W_spatioTemp[ , areas.to.disconnect] <- 0
        }

        cat("Evaluating windows with greedy-scanStatistic \n")
        tiempo <- system.time({
                Zs <- getClusters.SpatioTemporalScanStat.greedy_Parallel(data=data,
                                                                         geoDist=geoDist,
                                                                         maxAreas=maxAreas,
                                                                         timeWindowSize = timeWindowSize,
                                                                         W_spatioTemp=W_spatioTemp,
                                                                         verbose=verbose,
                                                                         maxCPUs = maxCPUs,
                                                                         computationCluster = computationCluster)
        })
        cat("Time evaluating areas: ", tiempo['elapsed']/60, " min\n")


        cat("\n\nPerforming MC test:\n")
        lambdaDist <- mcTest.logLambda.poisson_Parallel(data=data,
                                                        # zones_maxNeighbourhood=zones_maxNeighbourhood,
                                                        geoDist=geoDist,
                                                        maxAreas=maxAreas,
                                                        timeWindowSize = timeWindowSize,
                                                        W_spatioTemp=W_spatioTemp,
                                                        rep=reps.mc,
                                                        verbose.steps=1,
                                                        computationCluster=computationCluster)

        cat("MonteCarlo Test finished, calculating p-values\n")
        pVal <- apply(Zs, MARGIN=1, FUN = function(x,lambdaDist){
                if(x$lowRiskCluster){
                        lambdaDist %>% dplyr::select(lowRisks) %>% filter(lowRisks > x$lambda) %>% nrow() /(reps.mc+1)
                }
                else{
                        lambdaDist %>% dplyr::select(highRisks) %>% filter(highRisks > x$lambda) %>% nrow() /(reps.mc+1)
                }
        }, lambdaDist=lambdaDist)


        Zs$pVal <- pVal

        # We combine the high-risk and low-risk clusters and sort them by p-value and then by size.
        # Among clusters with the same significance, we give priority to the larger ones.
        Zs <- Zs %>% mutate(clusterSize=sapply(areas, FUN=function(x){length(x)})) #%>% arrange(pVal, desc(clusterSize))

        return(Zs)
}
