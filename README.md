# Replication Code of "Exponential Random Graph Models for Dynamic Signed Networks: An Application to International Politics"

------------



Description: 

This is the accompanying code repository for the paper Fritz et al. (2024) - ''Exponential Random Graph Models for Dynamic Signed Networks: An Application to International Politics''.
The code folder contains the R code to replicate all computations of the paper together with an R package implementing the SERGM.
The replication package includes all data needed to reproduce the results of the paper. 
The sources of the data are described in Section 3.1 of the Supplemental Materials. The script includes the preprocessing starting with the raw data. 
For both applications covered in the publication, the scripts rely on two packages to gather the data, namely, 'peacesciencer' and 'signnet'. To guarantee full reproducability in the face of package changes, the current state () of the preprocessed data are saved in ".RDS" files. 

Preparations: 

At first, the following packages need to be installed: 
"RcppArmadillo", "ergm", "stringr", "lpSolveAPI","signnet", "Rcpp",
"ggraph", "igraph", "tidygraph", "ggpubr", "coda", "data.table",
"ggmcmc","gridExtra","grid","Hmisc". To install them, execute the following command:
install.packages("signnet") 

The entire replication package is written in R and to be executed in RStudio.
Before running any script, the package "ergm.sign" needs to be installed. 
The package is provided as a folder. First, set the working directory of the 
RStudio session to be the main folder of the replication package.
Then execute the following command:

install.packages("ergm.sign_1.0.tar.gz")

Scripts:

There are three separate folders for the two applications ('International Cooperation and Conflict' from Section 4 in the main mauscript and 'Enmity and Friendship among New Guinean Highland Tribes' from Section 4 in the Supplemental Material)
and the simulation study (detailed in Section 2 in the Supplemental Material and 
split into two scripts according to the Setting 1 and 2).
Accordingly, there are separate files to run the respective analyes. 
The "other_functions.R" script includes routines to ease plotting the results and other miscellaneous functions. 
Following, the runner scripts of each folder are named with the needed computational resources and approximated time for execution: 

1. "application_dca.R". Needed resources: Computer with around 30 GB of RAM that allows to parallelize with 10 cores. Estimated time: 3 hours 
2. "application_tribes.R". Needed resources: Computer with around 10 GB of RAM that allows to parallelize with 20 cores. Estimated time: 1 hour 
3.1. "simulation_mle" (Setting 1 in Section 2 of the Supplemental Material). Needed resources: Computer with around 100 GB of RAM that allows to parallelize with 10 cores. Estimated time: 2 days 
3.2. "simulation_cp" (Setting 2 in Section 2 of the Supplemental Material). Needed resources: Computer with around 50 GB of RAM that allows to parallelize with 50 cores. Estimated time: 1 hour

Since the full replication takes around 2 days, all results are saved as ".RDS" files. 

Executed on a server with the following info: 
> sessionInfo()
R version 4.3.3 (2024-02-29)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.9 (Ootpa)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] vctrs_0.6.3        cli_3.6.2          ggmcmc_1.5.1.1     rlang_1.1.3        GGally_2.1.2       purrr_1.0.2        generics_0.1.3    
 [8] glue_1.7.0         colorspace_2.1-0   plyr_1.8.9         scales_1.3.0       fansi_1.0.6        grid_4.3.3         munsell_0.5.0     
[15] tibble_3.2.1       lifecycle_1.0.4    compiler_4.3.3     dplyr_1.1.3        RColorBrewer_1.1-3 Rcpp_1.0.12        pkgconfig_2.0.3   
[22] tidyr_1.3.0        rstudioapi_0.15.0  R6_2.5.1           tidyselect_1.2.0   utf8_1.2.4         pillar_1.9.0       magrittr_2.0.3    
[29] tools_4.3.3        gtable_0.3.4       reshape_0.8.9      ggplot2_3.4.3  

