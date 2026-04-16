# Doubly valid and doubly sharp sensitivity analysis to unobserved confounding for survival outcomes

This repository contains the R implementation of [Doubly valid and doubly sharp sensitivity analysis to unobserved confounding for survival outcomes](https://arxiv.org/). **IMPORTANT: For licensing reasons, some parts of the code were changed or removed (random forests, quantile regression, survival function, and parallel computing). Therefore, the present code is NOT FUNCTIONAL. To get a working code, see instructions in src/utils.R.**

* [Description](#description)
* [Requirements and Licenses](#requirements-and-licenses)
* [Output of the `sessionInfo()` command](#output-of-the-sessioninfo-command)

## Description

This repository is related to the paper *Doubly valid and doubly sharp sensitivity analysis to unobserved confounding for survival outcomes*, Baitairian et al. (2026). It is organized as follows:

- **dvds_bounds_survival_outcomes.Rmd**: an R Markdown document to reproduce the figures and numerical results from the paper. Users need to modify **src/utils.R** before executing code from this file;
- **src/utils.R**: functions that must be redefined by the user to get a working implementation and other useful functions;
- **src/data_simulation_fun.R**: functions to create the simulated datasets;
- **src/data_importation.R**: a file to preprocess the real data;
- **src/mc_data_generation.R**: a file to generate several Monte Carlo datasets;
- **src/dvds_fun.R**: functions used to compute the proposed doubly valid and doubly sharp bounds;
- **src/lee_zsb_fun.R**: functions used to compute the bounds from Lee et al. (2024);
- **src/mc_experiments.R**: a file to execute the experiments from the paper on the simulated and real data;
- **/figures**: a folder to save the generated figures;
- **/data/simul**: a folder that contains the simulated Monte Carlo datasets;
- **/data/RHC**: a folder that contains the real data on Right Heart Catheterization (RHC) (Connors et al., 1996);
- **/results**: a folder that contains the results of the sensitivity analyses.

*Note:* the direct optimization bounds from Lee et al. (2024) were reimplemented.

### References

Connors, A. F., Speroff, T., Dawson, N. V., Thomas, C., Harrell, F. E., Wagner, D., ... & Knaus, W. A. (1996). The effectiveness of right heart catheterization in the initial care of critically III patients. Jama, 276(11), 889-897. [Link](https://doi.org/10.1001/jama.1996.03540110043030)

Lee, S., Park, J. H., & Lee, W. (2024). Sensitivity analysis for unmeasured confounding in estimating the difference in restricted mean survival time. Statistical Methods in Medical Research, 33(11-12), 1979-1992. [Link](https://doi.org/10.1177/09622802241280782)

## Requirements and Licenses

The following libraries are used in this repository with R version 4.3.2.

|Library|Version|License|
|---|---|---|
|`fastDummies`|1.7.3|MIT|
|`foreach`|1.5.2|Apache License 2.0|
|`ggplot2`|3.4.4|MIT|

## Output of the `sessionInfo()` command

```
R version 4.3.2 (2023-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 8.7 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.15.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] fastDummies_1.7.3 ggplot2_3.4.4     foreach_1.5.2    

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5       cli_3.6.2         knitr_1.45        rlang_1.1.2       xfun_0.41         generics_0.1.3    doSNOW_1.0.20     labeling_0.4.3   
 [9] glue_1.6.2        colorspace_2.1-0  htmltools_0.5.7   snow_0.4-4        fansi_1.0.6       scales_1.3.0      rmarkdown_2.25    grid_4.3.2       
[17] evaluate_0.23     munsell_0.5.0     tibble_3.2.1      fastmap_1.1.1     yaml_2.3.8        lifecycle_1.0.4   compiler_4.3.2    dplyr_1.1.4      
[25] codetools_0.2-19  pkgconfig_2.0.3   rstudioapi_0.15.0 farver_2.1.1      digest_0.6.34     R6_2.5.1          tidyselect_1.2.0  utf8_1.2.4       
[33] pillar_1.9.0      parallel_4.3.2    magrittr_2.0.3    withr_2.5.2       tools_4.3.2       gtable_0.3.4      iterators_1.0.14 
```