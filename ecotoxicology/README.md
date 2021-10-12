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

    library(rjags)

    ## Le chargement a nécessité le package : coda

    ## Linked to JAGS 4.3.0

    ## Loaded modules: basemod,bugs

    model1 <-
      "
      model {
      

      for (i in 1:N)
    {
        N_alive[i] ~ dbin(p[i], n)
        p[i] <- 1 / (1 + (conc[i]/LC50)^b)
    }

    b ~ dunif(0, 10)
    logLC50 ~ dunif(-1,1)
    LC50 <- 10^logLC50

    }
      "

## Data

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

## Initial values

Start values need to be in the fixed interval of prior distribution:

    init <- list(
      list(logLC50 = -0.99, b = 1),
      list(logLC50 = 0, b = 5),
      list(logLC50 = 0.99, b = 9))

## Implementation

### Simulations

    m <- jags.model(file=textConnection(model1),
                    data = data_jags,
                    inits = init,
                    n.chains = 3
                    )

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 6
    ##    Unobserved stochastic nodes: 2
    ##    Total graph size: 45
    ## 
    ## Initializing model

    update(m, 3000)
    mcmc1 <- coda.samples(m, c("logLC50", "b"), n.iter = 5000)

### Minimal check of convergence

    plot(mcmc1)

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    require(lattice)

    ## Le chargement a nécessité le package : lattice

    xyplot(mcmc1)

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-6-2.png)

    summary(mcmc1)

    ## 
    ## Iterations = 4001:9000
    ## Thinning interval = 1 
    ## Number of chains = 3 
    ## Sample size per chain = 5000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##            Mean      SD  Naive SE Time-series SE
    ## b        1.4852 0.28574 0.0023331      0.0033346
    ## logLC50 -0.3639 0.07995 0.0006528      0.0009172
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##            2.5%     25%     50%     75%   97.5%
    ## b        0.9687  1.2861  1.4734  1.6728  2.0747
    ## logLC50 -0.5308 -0.4141 -0.3589 -0.3092 -0.2195

    gelman.diag(mcmc1) # return only values of 1 for convergence

    ## Potential scale reduction factors:
    ## 
    ##         Point est. Upper C.I.
    ## b                1          1
    ## logLC50          1          1
    ## 
    ## Multivariate psrf
    ## 
    ## 1

    gelman.plot(mcmc1)

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-8-1.png)
\#\#\# Autocorrelation plot

    autocorr.plot(mcmc1[[1]])

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-9-1.png)

    effectiveSize(mcmc1)

    ##        b  logLC50 
    ## 7424.937 7627.755

    raftery.diag(mcmc1)

    ## [[1]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                                
    ##          Burn-in  Total Lower bound  Dependence
    ##          (M)      (N)   (Nmin)       factor (I)
    ##  b       8        10824 3746         2.89      
    ##  logLC50 6        6406  3746         1.71      
    ## 
    ## 
    ## [[2]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                                
    ##          Burn-in  Total Lower bound  Dependence
    ##          (M)      (N)   (Nmin)       factor (I)
    ##  b       5        6078  3746         1.62      
    ##  logLC50 6        6756  3746         1.80      
    ## 
    ## 
    ## [[3]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                                
    ##          Burn-in  Total Lower bound  Dependence
    ##          (M)      (N)   (Nmin)       factor (I)
    ##  b       5        5771  3746         1.54      
    ##  logLC50 7        7397  3746         1.97

## Evolution of MCMC quantiles over iterates

    cumuplot(mcmc1[[1]])

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-10-1.png)

### Simulations without observed data (Monte Carlo)

    d0 <- list(conc = conc, n = n, N = N)
    m0 <- jags.model(file = textConnection(model1), data = d0, n.chains = 1)

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 0
    ##    Unobserved stochastic nodes: 8
    ##    Total graph size: 45
    ## 
    ## Initializing model

    update(m0, 3000)
    mc0 <- coda.samples(m0, c("b","LC50"), n.iter = 3000)

## Description of prior and posterior marginal distributions

It is a raw and easy-to-get representation, but poorly suited to uniform
distributions (e.g. b prior) which would be better represented using the
theoretical density rather than using the kernel density estimation from
a sample.

    ## Plot of prior densities
    plot(mc0, trace = FALSE, density = TRUE)

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-12-1.png)

    summary(mc0)

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
    ## LC50 2.131 2.465   0.0450         0.0450
    ## b    4.967 2.859   0.0522         0.0522
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%   75% 97.5%
    ## LC50 0.1112 0.3112 0.9873 3.145 8.672
    ## b    0.2325 2.4843 4.9327 7.484 9.729

    ## Plot of posterior densities
    plot(mcmc1, trace = FALSE, density = TRUE)

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-14-1.png)

    summary(mcmc1)

    ## 
    ## Iterations = 4001:9000
    ## Thinning interval = 1 
    ## Number of chains = 3 
    ## Sample size per chain = 5000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##            Mean      SD  Naive SE Time-series SE
    ## b        1.4852 0.28574 0.0023331      0.0033346
    ## logLC50 -0.3639 0.07995 0.0006528      0.0009172
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##            2.5%     25%     50%     75%   97.5%
    ## b        0.9687  1.2861  1.4734  1.6728  2.0747
    ## logLC50 -0.5308 -0.4141 -0.3589 -0.3092 -0.2195

    mctot <- as.data.frame(as.matrix(mcmc1))
    mcsample <- mctot[sample.int(nrow(mctot), size = 500), ]
    ## Plot of the joint posterior distribution as a scatter plot
    pairs(mcsample)

![](/home/tommaso/M2/Semestre_3/Stat_Bay/repo_git/ecotoxicology/README_files/figure-markdown_strict/unnamed-chunk-16-1.png)
