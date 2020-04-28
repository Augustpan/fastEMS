# fastEMS -- faster analsyis of metacommunity structure with NGS data

- author: Yuanfei Pan
- mail: yfpan16@fudan.edu.cn
- ref: Yuanfei Pan (2020). fastEMS: faster analsyis of metacommunity structure with NGS data. https://github.com/Augustpan/fastEMS

fastEMS.jl is a Julia re-implementation of the EMS analysis (elements of meta-community structure, Presley et al. 2010). This script is more than 100x faster than the original MATLAB implementation, and also much faster than the implementation in R package `metacom` (Dallas, 2020). As is known to me, this is the **only** software that is capable of doing EMS analysis using next-generation sequencing data.


## Reference:
Presley, S. J. , Higgins, C. L. , & Willig, M. R. . (2010). A comprehensive framework for the evaluation of metacommunity structure. Oikos, 119(6), 908-917.

Tad Dallas (2020). metacom: Analysis of the 'Elements of Metacommunity Structure'. R package version 1.5.3. https://CRAN.R-project.org/package=metacom
