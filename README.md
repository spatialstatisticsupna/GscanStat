# Greedy Scan Statistic (GscanStat) algorithm

This repository contains the R code and examples to fit the models described in the paper entitled "Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities" (Santafé et al., 2025)

## Table of Contents
- [R Code - GscanStat algorithm](#R-code)
- [R Code - Example](#Example)
- [Acknowledgements - Example](#Example)
- [References](#References)




# R code
R code to fit the following spatio-temporal disease mapping models:
- SaTScan model (Kulldorff, 2021)
- GscanStat model (Santafe et al., 2025)
has been included [here](https://github.com/spatialstatisticsupna/GscanStat_article/tree/main/R).

The [BigDM](https://github.com/spatialstatisticsupna/bigDM) library is used to fit spatio-temporal models. Additionally, the divide-and-conquer approach (Orozco-Acosta et al., 2023) implemented in the [BigDM](https://github.com/spatialstatisticsupna/bigDM) library is used when dealing with large risk maps.

# Example
An example to fit SaTScan and GscanStat models using a simulated spatio-temporal risk map for the autonomous community of Navarre (Spain) has been included [here](https://github.com/spatialstatisticsupna/GscanStat_article/tree/main/R). The [example_GscanStat.R](https://github.com/spatialstatisticsupna/GscanStat_article/blob/main/R/example_GscanStat.R) file also includes examples to fit spatio-temporal models using the divide-and-conquer approach.

# Acknowledgements
This work has been supported by the Spanish Ministry of Science and Innovation - State Research Agency (PID2020-113125RB-I00). 2021-2025

# References
[Kulldorff, M (2001). Prospective time-periodic geographical disease surveillance using a scan statistic. Journal of the Royal Statistical Society, A164:61-72.](https://www.jstor.org/stable/pdf/2680534)

[Orozco-Acosta, E., Adin, A., and Ugarte, M.D. (2023). Big problems in spatio-temporal disease mapping: methods and software. Computer Methods and Programs in Biomedicine, 231, 107403.](https://doi.org/10.1016/j.cmpb.2023.107403)

[Santafé, G., Adin, A., and Ugarte, M.D. (2025). Improving Disease Risk Estimation in Small Areas by Accounting for Spatio-Temporal Local Discontinuities.]  
