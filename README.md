# Greedy Scan Statistic (GscanStat) algorithm

This repository contains the R code for implementing the GscanStat clustering algorithm and reproducing the results from the case study described in the paper *“Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities”* (Santafé et al., 2025).

## Table of Contents

-   [Data](#data)
-   [R code](#r-code)
-   [Acknowledgements](#Acknowledgements)
-   [References](#References)

# Data

Overall cancer (all sites) mortality data for the male population in continental Spain (excluding the Balearic and Canary Islands) are provided. The data are aggregated into 3-year periods spanning 1999–2022 (i.e., 1999–2001, 2002–2004, …, 2020–2022).

In addition, to avoid introducing substantial variability that could hinder the performance of the GscanStat algorithm in detecting clusters, neighboring municipalities were aggregated (while respecting the administrative boundaries of the Autonomous Communities) until each resulting spatial unit contained at least 16 observed cases over the entire study period.

The [**DataSpain_Cancer.Rdata**](https://github.com/spatialstatisticsupna/GscanStat/blob/master/Data/DataSpain_Cancer.Rdata) file contains both the data and the cartographic information corresponding to the final configuration, which comprises 2.470 regions.

This .Rdata file contains the following objects:

-   `carto`: `sf` object containing the polygon geometries of the spatial units

-   `data`: `tibble` object with 19.760 rows and 5 columns

    -   `ID`: character vector with the IDs of the spatial units
    -   `year`: numeric vector representing the time period (1=1999-2001, ..., 8=2020-2022)
    -   `obs`: observed number of cancer deaths
    -   `exp`: expected number of cancer deaths (calculated using internal age-standardization)
    -   `pop`: population at risk

# R code

[Here](https://github.com/spatialstatisticsupna/GscanStat/tree/main/R) we provide R code to fit the following spatio-temporal models:

-   SaTScan model (Kulldorff, 2021)

-   GscanStat model (Santafé et al., 2025)

The [`bigDM`](https://github.com/spatialstatisticsupna/bigDM) package is used to fit local spatio-temporal models through a divide-and-conquer strategy (Orozco-Acosta et al., 2023), incorporating the significant risk clustering structure identified by the GscanStat algorithm. Version 0.5.7 of `bigDM` has been developed specifically for this purpose.

Two examples are provided:

-   [Example1_GscanStat_Navarre.R](https://github.com/spatialstatisticsupna/GscanStat/tree/main/R/Example1_GscanStat_Navarre.R) illustrates a small-scale data analysis for the Autonomous Region of Navarre (58 areas across 8 time points). First, the GscanStat and SaTScan algorithms are applied to detect significant clusters. Subsequently, spatio-temporal models with BYM2 spatial prior, RW1 temporal prior and Type IV interaction are fitted to estimate relative risks.

-   [Example2_GscanStat_Spain.R](https://github.com/spatialstatisticsupna/GscanStat/tree/main/R/Example2_GscanStat_Spain.R) demonstrates a large-scale data analysis for the entire Spanish regions (2.470 areas across 8 time points). First, the GscanStat and SaTScan algorithms are applied to detect significant clusters. Subsequently, local spatio-temporal models are fitted using the divide-and-conquer strategy implemented in the bigDM package to estimate relative risks.

    **WARNING:** This analysis is computationally intensive. As a reference, obtaining the final clustering partition with the GscanStat algorithm takes approximately 12.2 hours on an Intel(R) Xeon(R) Silver 4316 processor with 80 CPUs at 2.30 GHz and 256 GB of RAM.

Additional scripts:

-   [GscanStat_parallelClusteringAlgorithm.R](https://github.com/spatialstatisticsupna/GscanStat/tree/main/R/GscanStat_parallelClusteringAlgorithm.R) contains the functions needed to run the Greedy Scan Statistics (GscanStat) algorithm for cluster detection.

-   [SaTScan_auxFunctions.R](https://github.com/spatialstatisticsupna/GscanStat/tree/main/R/SaTScan_auxFunctions.R) provides auxiliary functions to running the SaTScan software via the `rsatscan` package.

# Acknowledgements

This work has been supported by project PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033 (Spanish Ministry of Science, Innovation and Universities, AEI).

![plot](https://github.com/spatialstatisticsupna/GscanStat/blob/main/miciu-aei.png)


# References

[Kulldorff, M (2001). Prospective time-periodic geographical disease surveillance using a scan statistic. *Journal of the Royal Statistical Society: Series A (Statistics in Society)*, 164, 61-72.](https://www.jstor.org/stable/pdf/2680534)

[Orozco-Acosta, E., Adin, A., and Ugarte, M.D. (2023). Big problems in spatio-temporal disease mapping: methods and software. *Computer Methods and Programs in Biomedicine*, 231, 107403.](https://doi.org/10.1016/j.cmpb.2023.107403)

Santafé, G., Adin, A., and Ugarte, M.D. (2025). Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities. *arXiv preprint*.
