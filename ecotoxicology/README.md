Ecotoxicology exemple
================

# Example of Bayesian inference on a model in ecotoxicology

Our aim is to study the effect of a toxic substance suspected to pollute
lakes and rivers. We wish to model the effect of this substance on the
death rate of daphnids (small invertebrates of fresh water, called
“water feas”). An in vitro experiment has been conducted to observe the
effect of the pollutant concentration on the death rate of 20 organisms
after an exposition of 21 days. The data for this experiment:

The tested concentrations (µg.ml<sup>-1</sup>): 0.19 0.38 0.76 1.53 3.05
6.11 The number of survivors among the 20 organisms: 16 12 4 3 1 1

We wish to estimate **LC50**, the concentration under which 50% of the
organisms are dead after 21 days, with a log-logistic modeling of the
21-days survival probability, through the formula:
![Formule](./assets/formule.png)

## Model formalization

![DAG model](./assets/dag_ecotoxicology.png)

``` r
library(rjags)
```

    ## Le chargement a nécessité le package : coda

    ## Linked to JAGS 4.3.0

    ## Loaded modules: basemod,bugs

``` r
model1 <-
  "
  model {
  

  for (i in 1:N)
{
    N_alive[i] ~ dbin(p[i], n)
    p[i] <- 1 / (1 + (conc[i]/LC50)^b)
}

b ~ dunif(0, 10)
LC50 <- 10^logLC50
logLC50 ~ dunif(-1,1)

}
  "
```

## Data

``` r
  conc <- c(0.19, 0.38, 0.76, 1.53, 3.05, 6.11)
  N_alive <- c(16, 12, 4, 3, 1, 1)
  N <- 6
  n <- 20

data_jags <- list(
  conc = conc,
  N_alive = N_alive,
  N = N,
  n = n
             )
```

## Initial values

Start values need to be in the fixed interval of prior distribution:

``` r
init <- list(
  list(logLC50 = -0.99, b = 1),
  list(logLC50 = 0, b = 5),
  list(logLC50 = 0.99, b = 9))
```

## Implementation

### Simulations

``` r
m <- jags.model(file=textConnection(model1),
                data = data_jags,
                inits = init,
                n.chains = 3
                )
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 6
    ##    Unobserved stochastic nodes: 2
    ##    Total graph size: 45
    ## 
    ## Initializing model

``` r
update(m, 3000)
mcmc1 <- coda.samples(m, c("LC50", "b"), n.iter = 5000)
```

### Minimal check of convergence

``` r
plot(mcmc1)
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
require(lattice)
```

    ## Le chargement a nécessité le package : lattice

``` r
xyplot(mcmc1)
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
summary(mcmc1)
```

    ## 
    ## Iterations = 4001:9000
    ## Thinning interval = 1 
    ## Number of chains = 3 
    ## Sample size per chain = 5000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##        Mean      SD  Naive SE Time-series SE
    ## LC50 0.4382 0.07878 0.0006432      0.0009022
    ## b    1.4806 0.28452 0.0023231      0.0035046
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%    75% 97.5%
    ## LC50 0.2907 0.3846 0.4361 0.4876 0.604
    ## b    0.9600 1.2819 1.4667 1.6644 2.067

``` r
gelman.diag(mcmc1) # return only values of 1 for convergence
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## LC50          1          1
    ## b             1          1
    ## 
    ## Multivariate psrf
    ## 
    ## 1

``` r
gelman.plot(mcmc1)
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
\#\#\# Autocorrelation plot

``` r
autocorr.plot(mcmc1[[1]])
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
effectiveSize(mcmc1)
```

    ##     LC50        b 
    ## 7644.493 6595.004

``` r
raftery.diag(mcmc1)
```

    ## [[1]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                             
    ##       Burn-in  Total Lower bound  Dependence
    ##       (M)      (N)   (Nmin)       factor (I)
    ##  LC50 7        7263  3746         1.94      
    ##  b    8        11888 3746         3.17      
    ## 
    ## 
    ## [[2]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                             
    ##       Burn-in  Total Lower bound  Dependence
    ##       (M)      (N)   (Nmin)       factor (I)
    ##  LC50 6        6520  3746         1.74      
    ##  b    5        5916  3746         1.58      
    ## 
    ## 
    ## [[3]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                             
    ##       Burn-in  Total Lower bound  Dependence
    ##       (M)      (N)   (Nmin)       factor (I)
    ##  LC50 7        7968  3746         2.13      
    ##  b    7        7263  3746         1.94

## Evolution of MCMC quantiles over iterates

``` r
cumuplot(mcmc1[[1]])
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Simulations without observed data (Monte Carlo)

``` r
d0 <- list(conc = conc, n = n, N = N)
m0 <- jags.model(file = textConnection(model1), data = d0, n.chains = 1)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 0
    ##    Unobserved stochastic nodes: 8
    ##    Total graph size: 45
    ## 
    ## Initializing model

``` r
update(m0, 3000)
mc0 <- coda.samples(m0, c("b","LC50"), n.iter = 3000)
```

## Description of prior and posterior marginal distributions

It is a raw and easy-to-get representation, but poorly suited to uniform
distributions (e.g. b prior) which would be better represented using the
theoretical density rather than using the kernel density estimation from
a sample.

``` r
## Plot of prior densities
plot(mc0, trace = FALSE, density = TRUE)
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
summary(mc0)
```

    ## 
    ## Iterations = 3001:6000
    ## Thinning interval = 1 
    ## Number of chains = 1 
    ## Sample size per chain = 3000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##       Mean    SD Naive SE Time-series SE
    ## LC50 2.162 2.511  0.04585        0.04741
    ## b    4.992 2.880  0.05258        0.05258
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%   50%   75% 97.5%
    ## LC50 0.1109 0.3184 0.986 3.206 8.866
    ## b    0.2279 2.4891 5.002 7.432 9.753

``` r
## Plot of posterior densities
plot(mcmc1, trace = FALSE, density = TRUE)
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
summary(mcmc1)
```

    ## 
    ## Iterations = 4001:9000
    ## Thinning interval = 1 
    ## Number of chains = 3 
    ## Sample size per chain = 5000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##        Mean      SD  Naive SE Time-series SE
    ## LC50 0.4382 0.07878 0.0006432      0.0009022
    ## b    1.4806 0.28452 0.0023231      0.0035046
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%    75% 97.5%
    ## LC50 0.2907 0.3846 0.4361 0.4876 0.604
    ## b    0.9600 1.2819 1.4667 1.6644 2.067

## Joint posterior probability

``` r
mctot <- as.data.frame(as.matrix(mcmc1))
mcsample <- mctot[sample.int(nrow(mctot), size = 500), ]
## Plot of the joint posterior distribution as a scatter plot
pairs(mcsample)
```

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
b <- mctot[, "b"]
LC50 <- mctot[, "LC50"]

p <- 1/(1+(1.53/LC50)^b) # mediane et intervalle
```

## Confrontation of the model to the data - posterior predictive check

``` r
summary(mcsample)
```

    ##       LC50              b         
    ##  Min.   :0.2258   Min.   :0.7048  
    ##  1st Qu.:0.3875   1st Qu.:1.2964  
    ##  Median :0.4393   Median :1.4930  
    ##  Mean   :0.4439   Mean   :1.4770  
    ##  3rd Qu.:0.4946   3rd Qu.:1.6486  
    ##  Max.   :0.7400   Max.   :2.2801

``` r
dic.samples(m, n.iter = 5000, type = "pD")
```

    ## Mean deviance:  20.39 
    ## penalty 2.037 
    ## Penalized deviance: 22.43
