
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Semi-Distance Correlation and MV Index: Measure Dependence Between Categorical and Continuous Variables

<!-- badges: start -->
<!-- badges: end -->

The goal of package `semidist` is to provide an easy way to implement
the semi-distance methods (Zhong et al., 2022+) and MV index methods
(Cui, Li and Zhong, 2015; Cui and Zhong, 2019).

## Installation

To install `semidist`,

``` r
install.packages("semidist")
```

## Example

Here is a simple example showing how to use `semidist` to measure the
dependence between a categorical variable and a multivariate continuous
variable, and apply the measure on testing the independence and conduct
groupwise feature screening.

``` r
library(semidist)
X <- mtcars[, c("mpg", "disp", "drat", "wt")]
y <- factor(mtcars[, "am"])

sdcov(X, y)
#> [1] 31.78288
sdcor(X, y)
#> [1] 0.3489821

sd_test(X, y)
#> $method
#> [1] "Semi-Distance Independence Test (Permutation Test with K = 10000)"
#> 
#> $name_data
#> [1] "X and y"
#> 
#> $n
#> [1] 32
#> 
#> $test_type
#> [1] "perm"
#> 
#> $stat
#> [1] 940.344
#> 
#> $pvalue
#> [1] 0.0009999
#> 
#> $num_perm
#> [1] 10000
#> 
#> attr(,"class")
#> [1] "indtest"

sd_sis(X, y, d = 2)
#> $group_info
#> $group_info$`Grp mpg`
#> [1] "mpg"
#> 
#> $group_info$`Grp disp`
#> [1] "disp"
#> 
#> $group_info$`Grp drat`
#> [1] "drat"
#> 
#> $group_info$`Grp wt`
#> [1] "wt"
#> 
#> 
#> $measurement
#>   Grp mpg  Grp disp  Grp drat    Grp wt 
#> 0.3447938 0.3488447 0.5054821 0.5358834 
#> 
#> $selected
#> [1] "wt"   "drat"
#> 
#> $ordering
#> [1] "Grp wt"   "Grp drat" "Grp disp" "Grp mpg"

# Suppose we have prior information for the group structure as
# ("mpg", "drat"), ("disp", "hp") and ("wt", "qsec")
group_info <- list(
  mpg_drat = c("mpg", "drat"),
  disp_wt = c("disp", "wt")
)
sd_sis(X, y, group_info, d = 2)
#> $group_info
#> $group_info$mpg_drat
#> [1] "mpg"  "drat"
#> 
#> $group_info$disp_wt
#> [1] "disp" "wt"  
#> 
#> 
#> $measurement
#>  mpg_drat   disp_wt 
#> 0.3518051 0.3488598 
#> 
#> $selected
#> [1] "mpg"  "drat"
#> 
#> $ordering
#> [1] "mpg_drat" "disp_wt"
```

## References

1.  Wei Zhong, Zhuoxi Li, Wenwen Guo and Hengjian Cui. (2022+)
    “Semi-Distance Correlation and Its Applications.” *Journal of the
    American Statistical Association.* Under review.
2.  Hengjian Cui and Wei Zhong (2019). “A Distribution-Free Test of
    Independence Based on Mean Variance Index.” *Computational
    Statistics & Data Analysis*, 139, 117-133.
3.  Hengjian Cui, Runze Li and Wei Zhong (2015). “Model-Free Feature
    Screening for Ultrahigh Dimensional Discriminant Analysis.” *Journal
    of the American Statistical Association*, 110, 630-641.
