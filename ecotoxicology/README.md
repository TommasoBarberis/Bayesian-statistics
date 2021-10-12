# Example of Bayesian inference on a model in ecotoxicology

## Context

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
$p = \\frac{1}{1+(\\frac{conc}{LC50})^b}$

## Data

``` r
conc <- c(0.19, 0.38, 0.76, 1.53, 3.05, 6.11)
Nalive <- c(16, 12, 4, 3, 1, 1)
n <- 20 # number of organisms
  
plot(Nalive ~ conc, pch = 16, xlab = "concentration of toxic substance", ylab = "number of survivors")
```

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)

## Implementation

### Model formalization

![DAG model](./assets/dag_ecotoxicology.png)

``` r
library(rjags)
```

    ## Le chargement a nécessité le package : coda

    ## Linked to JAGS 4.3.0

    ## Loaded modules: basemod,bugs

``` r
desc_model <-
  "
  model {
  

  for (i in 1:N)
{
    Nalive[i] ~ dbin(p[i], n)
    p[i] <- 1 / (1 + (conc[i]/LC50)^b)
}

b ~ dunif(0, 10)
LC50 <- 10^logLC50
logLC50 ~ dunif(-1,1)

}
  "
```

## MCMC Simulation

### Data

``` r
data_jags <- list(
  conc = conc,
  Nalive = Nalive,
  N = length(conc),
  n = n
)
```

### Initial values

Start values need to be in the fixed interval of prior distribution:

``` r
init <- list(
  list(logLC50 = -0.99, b = 1),
  list(logLC50 = 0, b = 5),
  list(logLC50 = 0.99, b = 9))
```

### Simulations of Markov chains

``` r
model <- jags.model(file=textConnection(desc_model),
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
update(model, 3000) # burn-in phase
mcmc1 <- coda.samples(model, c("LC50", "b"), n.iter = 5000)
```

### Minimal check of convergence

``` r
plot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
require(lattice)
```

    ## Le chargement a nécessité le package : lattice

``` r
xyplot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-2.png)

``` r
gelman.diag(mcmc1) # return only values of 1 for adequate convergence 
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

![](README_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
geweke.diag(mcmc1)
```

    ## [[1]]
    ## 
    ## Fraction in 1st window = 0.1
    ## Fraction in 2nd window = 0.5 
    ## 
    ##     LC50        b 
    ## -1.48297 -0.02203 
    ## 
    ## 
    ## [[2]]
    ## 
    ## Fraction in 1st window = 0.1
    ## Fraction in 2nd window = 0.5 
    ## 
    ##  LC50     b 
    ## 1.402 0.455 
    ## 
    ## 
    ## [[3]]
    ## 
    ## Fraction in 1st window = 0.1
    ## Fraction in 2nd window = 0.5 
    ## 
    ##   LC50      b 
    ## 0.7104 0.7070

#### Evolution of MCMC quantiles over iterates

``` r
cumuplot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-1.png)![](README_files/figure-markdown_github/unnamed-chunk-10-2.png)![](README_files/figure-markdown_github/unnamed-chunk-10-3.png)

### Autocorrelation plot

``` r
acfplot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
autocorr.plot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-11-2.png)![](README_files/figure-markdown_github/unnamed-chunk-11-3.png)![](README_files/figure-markdown_github/unnamed-chunk-11-4.png)

``` r
effectiveSize(mcmc1) # estimation of effective size in function of autocorrelation for the 3 chains (3 x 5000 = 15000 iterations)
```

    ##     LC50        b 
    ## 7838.172 6928.176

``` r
# if we want reduce as possible autocorrelation, we can set the parameter thin = 15000/7000 (7000 => effective size)

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
    ##  LC50 6        6520  3746         1.74      
    ##  b    5        5673  3746         1.51      
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
    ##  LC50 6        6878  3746         1.84      
    ##  b    6        6406  3746         1.71      
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
    ##  LC50 6        7131  3746         1.90      
    ##  b    6        8820  3746         2.35

``` r
# if the dependance factor is grater than 5, it means that a strong autocorrelation is present in the chains
```

#### New simulation with adequate *n.iter* and *thin* values

-   *n.iter*: 10000 iterations (see *raftery.diag* results);
-   *thin*: 2 (see *effectiveSize* results)

``` r
# new simulations
mcmc1 <- coda.samples(model, c("LC50", "b"), n.iter = 10000, thin = 2)

# new autocorrelation test
acfplot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
autocorr.plot(mcmc1) 
```

![](README_files/figure-markdown_github/unnamed-chunk-12-2.png)![](README_files/figure-markdown_github/unnamed-chunk-12-3.png)![](README_files/figure-markdown_github/unnamed-chunk-12-4.png)

``` r
xyplot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

### Parallelization of MCMC chains

Here, the parallelization is not necessary, but it is interesting when
**MCMC** require a huge number of iterations.

#### No-parallel version

``` r
# t1 <- Sys.time()
# model2 <- jags.model(
#   file = textConnection(desc_model),
#   data = data_jags, 
#   inits = init, 
#   n.chains = 3)
# 
# update(model, 5000)
# 
# mcmc <- coda.samples(model2, c("LC50", "b"), n.iter = 200000, thin = 40)
# t2 <- Sys.time()
```

``` r
# # computational time
# t2 - t1
```

#### Parallel version

``` r
# library(dclone)
# # cl <- makePSOCKcluster(3)
# # sous linux ou MAC ce sera plus efficace d utiliser makeForkCluster
# cl <- makeForkCluster(3)
# t1 <- Sys.time()
# 
# parJagsModel(
#   cl, 
#   name = "modelpar", 
#   file = "desc_model.txt",
#   data = data_jags, 
#   inits = init, 
#   n.chains = 3)
# 
# parUpdate(cl, object = "modelpar", n.iter = 5000)
# mcmcpar <- parCodaSamples(cl, model = "modelpar",
# variable.names = c("LC50", "b"),
# n.iter = 200000, thin = 40)
# t2 <- Sys.time()
# stopCluster(cl)
```

## Mcmc utilisation

### Description of posterior marginal distribution

``` r
summary(mcmc1)
```

    ## 
    ## Iterations = 9002:19000
    ## Thinning interval = 2 
    ## Number of chains = 3 
    ## Sample size per chain = 5000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##        Mean      SD  Naive SE Time-series SE
    ## LC50 0.4376 0.08005 0.0006536      0.0007372
    ## b    1.4787 0.28793 0.0023509      0.0027156
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%    75%  97.5%
    ## LC50 0.2891 0.3836 0.4335 0.4881 0.6056
    ## b    0.9576 1.2780 1.4652 1.6628 2.0859

``` r
densityplot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
HPDinterval(mcmc1)
```

    ## [[1]]
    ##          lower     upper
    ## LC50 0.2767340 0.5985109
    ## b    0.9220456 2.0521601
    ## attr(,"Probability")
    ## [1] 0.95
    ## 
    ## [[2]]
    ##          lower     upper
    ## LC50 0.2865038 0.5977355
    ## b    0.9365833 2.0534086
    ## attr(,"Probability")
    ## [1] 0.95
    ## 
    ## [[3]]
    ##          lower     upper
    ## LC50 0.2797155 0.5879739
    ## b    0.9723578 2.0745191
    ## attr(,"Probability")
    ## [1] 0.95

### Description of posterior join distribution

``` r
crosscorr(mcmc1)
```

    ##           LC50         b
    ## LC50 1.0000000 0.3168721
    ## b    0.3168721 1.0000000

``` r
crosscorr.plot(mcmc1)
```

![](README_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
mctot <- as.data.frame(as.matrix(mcmc1))
mctotsample <- mctot[sample.int(nrow(mctot), size = 300), ]
## Plot of the joint posterior distribution as a scatter plot
pairs(mctotsample)
```

![](README_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
require(IDPmisc)
```

    ## Le chargement a nécessité le package : IDPmisc

``` r
ipairs(mctotsample)
```

![](README_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
panel.hist <- function(x, col.hist = "grey", ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, col=col.hist, ...)
}
panel.dens <- function(x, col.dens = "red", lwd.dens = 2, ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
densx <- density(x)
vx <- densx$x
vy <- densx$y
lines(vx,vy/max(vy),col=col.dens,lwd=lwd.dens, ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r <- abs(cor(x, y,method="spearman"))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#text(0.5, 0.5, txt, cex = cex.cor * r)
text(0.5, 0.5, txt, cex = cex.cor * 0.5, ...)
}
panel.xy <- function(x, y, pch.xy = 1, col.xy = "black", cex.xy = 0.5, ...)
{
points(x,y,pch=pch.xy, col=col.xy,cex=cex.xy, ...)
}
```

``` r
pairs(
  mctotsample, 
  upper.panel = panel.xy,
  diag.panel = panel.dens,
  lower.panel = panel.cor)
```

![](README_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
pairs(mctotsample, upper.panel = panel.xy,
diag.panel = panel.hist,
lower.panel = panel.cor,
col.hist = "blue", pch.xy = 1, col.xy = "red", cex.xy = 1)
```

    ## Warning in plot.window(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in title(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...): "pch.xy"
    ## n'est pas un paramètre graphique

    ## Warning in rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...): "col.xy"
    ## n'est pas un paramètre graphique

    ## Warning in rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...): "cex.xy"
    ## n'est pas un paramètre graphique

    ## Warning in plot.window(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in title(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.hist" n'est
    ## pas un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "pch.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "cex.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.hist" n'est
    ## pas un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "pch.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "cex.xy" n'est pas
    ## un paramètre graphique

    ## Warning in text.default(0.5, 0.5, txt, cex = cex.cor * 0.5, ...): "col.hist"
    ## n'est pas un paramètre graphique

    ## Warning in text.default(0.5, 0.5, txt, cex = cex.cor * 0.5, ...): "pch.xy" n'est
    ## pas un paramètre graphique

    ## Warning in text.default(0.5, 0.5, txt, cex = cex.cor * 0.5, ...): "col.xy" n'est
    ## pas un paramètre graphique

    ## Warning in text.default(0.5, 0.5, txt, cex = cex.cor * 0.5, ...): "cex.xy" n'est
    ## pas un paramètre graphique

    ## Warning in plot.window(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in title(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.hist" n'est
    ## pas un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "pch.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "cex.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.hist" n'est
    ## pas un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "pch.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "col.xy" n'est pas
    ## un paramètre graphique

    ## Warning in axis(side = side, at = at, labels = labels, ...): "cex.xy" n'est pas
    ## un paramètre graphique

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "col.hist" n'est pas un
    ## paramètre graphique

    ## Warning in plot.window(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.window(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.hist" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "col.xy" n'est pas un paramètre graphique

    ## Warning in plot.xy(xy, type, ...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.hist" n'est pas un paramètre graphique

    ## Warning in title(...): "pch.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "col.xy" n'est pas un paramètre graphique

    ## Warning in title(...): "cex.xy" n'est pas un paramètre graphique

    ## Warning in rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...): "pch.xy"
    ## n'est pas un paramètre graphique

    ## Warning in rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...): "col.xy"
    ## n'est pas un paramètre graphique

    ## Warning in rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...): "cex.xy"
    ## n'est pas un paramètre graphique

![](README_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
require(GGally)
```

    ## Le chargement a nécessité le package : GGally

    ## Le chargement a nécessité le package : ggplot2

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

``` r
ggscatmat(mctotsample)
```

![](README_files/figure-markdown_github/unnamed-chunk-25-1.png)

### Comparaison between *prior* and *posterior* laws

#### Simulations without observed data (Monte Carlo)

``` r
d0 <- list(conc = conc, n = n, N = length(conc))
model0 <- jags.model(file = textConnection(desc_model), data = d0, n.chains = 1)
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
update(model0, 5000)
mcmc0 <- coda.samples(model0, c("b","LC50"), n.iter = 5000)
mcmctot0 <- as.data.frame(as.matrix(mcmc0))

par(mfrow = c(1, 3))
par(mar = c(5, 2, 1, 1))
for (i in 1:ncol(mcmctot0))
{
hist(mctot[,i], main = "", xlab = names(mcmctot0)[i], freq = FALSE)
lines(density(mcmctot0[,i]), col = "blue", lwd = 2)
}
```

![](README_files/figure-markdown_github/unnamed-chunk-26-1.png)

### Using *posterior* distribution to estimate parameter

``` r
b <- mctot[, "b"]
LC50 <- mctot[, "LC50"]

var_conc <- 1.53  # ug/ml
p <- 1/(1+(var_conc/LC50)^b)
quantile(p, probs = c(0.025, 0.5, 0.975))
```

    ##       2.5%        50%      97.5% 
    ## 0.06674147 0.13707394 0.23337481

``` r
N_alive <- rbinom(n=length(p),
                  size=20,
                  prob=p)
quantile(N_alive, probs = c(0.025, 0.5, 0.975))
```

    ##  2.5%   50% 97.5% 
    ##     0     3     7

## Model validation

### DIC

``` r
dic.samples(model, n.iter = 5000, type = "pD")
```

    ## Mean deviance:  20.44 
    ## penalty 2.079 
    ## Penalized deviance: 22.52

### WAIC

``` r
desc_model_loo <-
  desc_model <-
  "
  model {
  

  for (i in 1:N)
{
    Nalive[i] ~ dbin(p[i], n)
    p[i] <- 1 / (1 + (conc[i]/LC50)^b)
    
    loglik[i] <- log(dpois(p[i], Nalive[i]))
}

b ~ dunif(0, 10)
LC50 <- 10^logLC50
logLC50 ~ dunif(-1,1)

}
  "
```

``` r
# model_loo <- jags.model(file = textConnection(desc_model_loo), data = data_jags,
# inits = init, n.chains = 3)
# 
# update(model_loo, 5000)
# mcmc_loo <- coda.samples(model_loo, c("loglik"), n.iter = 5000)
# loglik_loo <- as.matrix(mcmc_loo)
```
