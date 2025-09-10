# Greedy Scan Statistic (GscanStat) algorithm

This repository contains the R code and examples to fit the models described in the paper entitled "Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities" (Santafé et al., 2025)

## Table of Contents
- [R Code - GscanStat algorithm](#R-code)
- [R Code - Example](#Examples)
- [Acknowledgements - Example](#Example)
- [References](#References)




# R code
R code to fit the following spatio-temporal disease mapping models is provided:
- SaTScan model (Kulldorff, 2021)
- GscanStat model (Santafé et al., 2025)
  
The code can be found [here](https://github.com/spatialstatisticsupna/GscanStat_article/tree/main/R).

The [BigDM](https://github.com/spatialstatisticsupna/bigDM) package is used to fit spatio-temporal models. In addition, the divide-and-conquer approach (Orozco-Acosta et al., 2023) implemented in the [BigDM](https://github.com/spatialstatisticsupna/bigDM) package is applied when dealing with large risk maps.  

The file [satScan_auxFunctions.R](https://github.com/spatialstatisticsupna/GscanStat_article/blob/main/R/satScan_auxFunctions.R) provides auxiliary functions to run the SaTScan software within the simulation experiments described in the paper *"Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities"*. These functions can only be used with simulated data, since they do not implement internal age and sex standardization. They are not intended as general-purpose code for using SaTScan outside the simulated experiments presented in the paper.

# Examples
Two examples are provided [here](https://github.com/spatialstatisticsupna/GscanStat_article/tree/main/R):

* The first example corresponds to a small-scale problem where the SaTScan and GscanStat models are fitted using a simulated spatio-temporal risk map for the autonomous community of Navarre (Spain) at the municipality level (265 municipalities). The [example_GscanStat_Navarre.R](https://github.com/spatialstatisticsupna/GscanStat_article/blob/main/R/example_GscanStat_Navarre.R) script provides this code and also includes examples of fitting spatio-temporal models using the divide-and-conquer approach.
* The second example corresponds to a large-scale problem where the GscanStat model is fitted using a simulated spatio-temporal risk map for Spain. In this case, small neighboring municipalities are aggregated into supramunicipality areas, as described in the paper *"Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities"* (Santafé et al., 2025), to mitigate the high variability observed in small municipalities with very low observed/expected death counts. The final map with aggregated areas consists of 2,470 areas. The [example_GscanStat_Spain.R](https://github.com/spatialstatisticsupna/GscanStat_article/blob/main/R/example_GscanStat_Spain_.R) script provides the code to fit spatio-temporal models using the divide-and-conquer approach.  
  **WARNING:** Since this is a large-scale problem, the GscanStat algorithm may take a long time to obtain the final clustering partition. As a reference, the algorithm requires approximately 12.2 hours to compute the clustering partition on an Intel(R) Xeon(R) Silver 4316 processor with 80 CPUs at 2.30 GHz and 256 GB of RAM.


# Acknowledgements
This work has been supported by the Spanish Ministry of Science and Innovation - State Research Agency (PID2020-113125RB-I00). 2021-2025

# References
[Kulldorff, M (2001). Prospective time-periodic geographical disease surveillance using a scan statistic. Journal of the Royal Statistical Society, A164:61-72.](https://www.jstor.org/stable/pdf/2680534)

[Orozco-Acosta, E., Adin, A., and Ugarte, M.D. (2023). Big problems in spatio-temporal disease mapping: methods and software. Computer Methods and Programs in Biomedicine, 231, 107403.](https://doi.org/10.1016/j.cmpb.2023.107403)

Santafé, G., Adin, A., and Ugarte, M.D. (2025). Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities.
