---
title: "NPBin"
author: "Anthony Aylward"
date: "2018-03-20"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



We will analyze a subset of one of the sample dataset for illustration
purposes.

```r
library(npbin)

minimum_coverage <- 5 # minimum total coverage allowed
n_cores <- 47 # the number of cores to be used, can ONLY be 1 if run on Windows.

dt.ct <- atac
colnames(dt.ct)
#>  [1] "chr"             "location"        "m"              
#>  [4] "xm"              "winning.chip"    "motif"          
#>  [7] "pval.mat.atSNP"  "pval.pat.atSNP"  "pval.rank.atSNP"
#> [10] "winnig.motif"    "potential_TP"    "potential_FP"
```

for illustration purpose, keep only the 2000 points, remove this line will
end up with analyzing the whole data set. It could be slow if only one core
is used.

```r
dt.ct <- dt.ct[1:2000, ]
```

NPBin

```r
n_breaks <- 11 # number of breaks
spline_order <- 4 # order of splines
breaks <- seq(0, 1, length.out = n_breaks)
pi_init <- hist(
  dt.ct[, p_hat],
  breaks = seq(0, 1, length.out = n_breaks + spline_order - 3),
  plot = FALSE
)[["density"]] # initialized the weights using the histogram of p_hat
#> Error in `[.data.frame`(dt.ct, , p_hat): object 'p_hat' not found
```

estimate the overall model

```r
overall_model_estimate <- emBinBspl(
  dt.ct[, xm],
  dt.ct[, m],
  breaks = breaks,
  k = spline_order,
  pi.init = pi_init,
  ncores = n_cores,
  err.max = 1e-3,
  iter.max = 200
)  
#> Error in `[.data.frame`(dt.ct, , xm): object 'xm' not found
```

estimate the null model

```r
null_model_estimate <- estNull(
  dt.ct[, xm],
  dt.ct[, m],
  overall_model_estimate,
  init = NULL,
  iter.max = 200,
  ncores = n_cores,
  ub = rep(log(1e4), 2),
  err.max = 1e-4
)
#> Error in estNull1(mod, pseq, ncores = ncores): object 'overall_model_estimate' not found
dt.ct[,
  fnp := null_model_estimate[["f"]]
][,
  f0np := null_model_estimate[["f0"]]
][,
  locfdrnp := null_model_estimate[["locfdr"]]
][,
  fdrnp := locfdr2FDR(locfdrnp)
][,
  ranknp := rank(locfdrnp, ties.method = "max")
]
#> Error in `:=`(fnp, null_model_estimate[["f"]]): could not find function ":="
names(null_model_estimate)
#> Error in eval(expr, envir, enclos): object 'null_model_estimate' not found
```


```r
null_model_estimate$coef.null
#> Error in eval(expr, envir, enclos): object 'null_model_estimate' not found
```

Empirical Bayes test using p_hat

```r
pct0 <- 0.45         
empirical_bayes_beta_hat <- ebBeta(
  dt.ct[, xm],
  dt.ct[, m],
  dt.ct[, p_hat],
  breaks = breaks,
  k = spline_order,
  pi.init = pi_init,
  pct0 = pct0,
  init = NULL,
  iter.max = 200,
  err.max = 1e-4,
  ncores = n_cores
)
#> Error in `[.data.frame`(dt.ct, , m): object 'm' not found
dt.ct[,
  fhat := empirical_bayes_beta_hat[["f"]]
][,
  f0hat := empirical_bayes_beta_hat[["f0"]]
][,
  locfdrhat := empirical_bayes_beta_hat[["locfdr"]]
][,
  fdrhat := locfdr2FDR(locfdrhat)
][,
  rankhat := rank(locfdrhat, ties.method = "max")
]
#> Error in `:=`(fhat, empirical_bayes_beta_hat[["f"]]): could not find function ":="
names(empirical_bayes_beta_hat)
#> Error in eval(expr, envir, enclos): object 'empirical_bayes_beta_hat' not found
```

null parameters of EBE

```r
empirical_bayes_beta_hat[["coef.null"]]
#> Error in eval(expr, envir, enclos): object 'empirical_bayes_beta_hat' not found
```

Binomial test

```r
p_binomial <- sapply(
  1:n,
  function(y) binom.test(dt.ct[y, xm], dt.ct[y, m])[["p.value"]]
)
#> Error in lapply(X = X, FUN = FUN, ...): object 'n' not found
dt.ct[,
  pvbin := p_binomial
][,
  fdrbin := p.adjust(pvbin, method = "BH")
][,
  rankbin := rank(pvbin, ties.method = "max")
]
#> Error in `:=`(pvbin, p_binomial): could not find function ":="
```

Evaluate the results using motifs

number of potential TP defined based on motif

```r
dt.ct[, sum(potential_TP)]
#> Error in `[.data.frame`(dt.ct, , sum(potential_TP)): object 'potential_TP' not found
```

number of potential FP defined based on motif

```r
dt.ct[, sum(potential_FP)]
#> Error in `[.data.frame`(dt.ct, , sum(potential_FP)): object 'potential_FP' not found
```

find the number of TP and FP in top ranked SNPs

```r
dt.ct[,
  tpnp := rank2nhit(ranknp ,potential_TP)
][,
  fpnp := rank2nhit(ranknp, potential_FP)
]
#> Error in `:=`(tpnp, rank2nhit(ranknp, potential_TP)): could not find function ":="
dt.ct[,
  tp_hat := rank2nhit(rankhat, potential_TP)
][,
  fp_hat := rank2nhit(rankhat, potential_FP)
]
#> Error in `:=`(tp_hat, rank2nhit(rankhat, potential_TP)): could not find function ":="
dt.ct[,
  tpbin := rank2nhit(rankbin, potential_TP)
][,
  fpbin := rank2nhit(rankbin, potential_FP)
]
#> Error in `:=`(tpbin, rank2nhit(rankbin, potential_TP)): could not find function ":="
```

plot the accuracy measure as in the main paper. 
We presented a zoom-in version in the main paper to the top 20%, 
because usually there are not many ALI SNPs.
Note that the default of the demo only select a subset of the data for
illustration purposes.
Thus the figure may not an exact replica of the one in the paper.
To reproduce the results in the paper, please use the whole dataset

```r
cbfpalette <- c(
  "#D55E00",
  "#0072B2",
  "#CC79A7",
  "#009E73",
  "#E69F00",
  "#56B4E9",
  "#F0E442"
)
plotidac <- c(' NPB','EBE','Binom')
```
