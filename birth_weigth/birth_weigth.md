# 1. Context

We wish to compare the birth weight between boys and girls. To do this,
we know the birth weight for a set of 48 children (file pnais48.txt).
Column POIDNAIS represents birth weight, and column SEXE represents sex
(0 : boy, 1 : girl). Variables are renamed as weight and sex. Moreover,
for ease of use, the sex variable is re-encoded by adding 1 at each
child (hence, the value 1 correspond to a boy, and 2 to a girl).

Hypothesysis: variability of the wight is the same between boys and
girls.

``` r
data_weight <- read.table("pnais48.txt", h = T, sep = ",")
sex <- data_weight$SEX + 1
weight <- data_weight$POIDNAIS
table(sex)
```

    ## sex
    ##  1  2 
    ## 29 19

The study is on 29 boys and 19 girls.

The weight of children is supposed to follow a normal distribution, with
different means between the two sexes (*μ*<sub>*b*</sub> and
*μ*<sub>*g*</sub>), but a shared standard deviation (*σ*). We set two
priors: a uniform distribution between `2500` and `5000` grams for the
mean weight for boys and girls and a uniform distribution between 200
and 800 grams for the standard deviation of the weight.

# 2. Model formalization

-   Set up the DAG associated to the model, specifying stochastic and
    deterministic links.

![](./assets/diag.png)

# 3. Implementation of the model

## Model

The model is implemented as a string, called `desc_model1`.

-   Beware: In ***JAGS***, the normal distribution is parameterized with
    the precision (1/*σ*2 = *τ*) and not the standard deviation. The
    birth weight of a child `i` follows a normal distribution (`dnorm`),
    centered on mean *μ*<sub>*i*</sub> (dependent on sex), and with
    standard deviation *σ* (or precision *τ* under ***JAGS***).  
    Data is made of N = 48 cases: a loop is used to define the
    distribution followed by each of these observations.

-   For ease of use, the mean weights of boys and girls will be stored
    in a vector called `means`, of length 2. Hence, `mean[1]` will
    correspond to the mean weight of boys, and `mean[2]` to the mean
    weight of girls. For the `i` th child, `mean[sex[i]]` will give the
    mean weight expected for the child. The mean used for the
    distribution associated to each observation depends on the sex:
    `mean[2]` for girls and `mean[1]` for boys. The link between the
    mean used for the distribution associated to an observation and the
    sex is deterministic (simple equality, with no random part).

``` r
library(rjags)
```

    ## Le chargement a nécessité le package : coda

    ## Linked to JAGS 4.3.0

    ## Loaded modules: basemod,bugs

``` r
desc_model1 <-
  "
  model {
  
  # Defining links
  
  for (i in 1:Nchild)
{
    mu[i] <- mean[sex[i]]   # mean: vector of two elements, mean for boys                                and means for girls
    w[i] ~ dnorm(mu[i], tau) # weight
}

  # Definition for a prior distribution  
  mean[1] ~ dunif(2500, 5000) 
  mean[2] ~ dunif(2500, 5000)
  tau <- 1/sd # precision
  sd ~ dunif(200, 800)
  }
  "
```

**\~**: stochastic link

**\<-**: deterministic link

## Data

-   Define the data required for this model (data). Beware : Do no
    forget to include in the loop the number of observations (N).

``` r
data_birth_w <- list(Nchild = length(data_weight$POIDNAIS),
             sex = sex,
             w = weight
             )
```

## Initial values

Start values need to be in the fixed interval of prior distribution:

``` r
init_birth_w <- list(
  list(mean = c(2600, 4000), sd = 500),
  list(mean = c(4500, 2700), sd = 750),
  list(mean = c(4000, 4000), sd = 250))
```

## Implementation

``` r
model_birth_w <- jags.model(file=textConnection(desc_model1),
                data = data_birth_w,
                inits = init_birth_w,
                n.chains = 3
                )
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 48
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 106
    ## 
    ## Initializing model

``` r
update(model_birth_w, 3000)
mcmc_birth_w <- coda.samples(model_birth_w, c("mean", "sd"), n.iter = 5000)
```

The table of chain `i` is obtained with the following command:

``` r
# mcmc_birth_w[[1]]
```

Each table contains the parameter as columns, with the iterations of the
***MCMC*** as lines.

Computation of the posterior mean of the average birth weight of boys
from the first ***MCMC***:

``` r
mean(mcmc_birth_w[[1]][, "mean[1]"])
```

    ## [1] 3727.995

⇒ On this example, it is not possible to get an explicit description of
the posterior distributions of the parameters. Yet, if *σ* is known,
*μ*<sub>*b*</sub> and *μ*<sub>*g*</sub> both follow a *normal*
distribution. Likewise, if *μ*<sub>*b*</sub> and *μ*<sub>*g*</sub> are
known, *σ* follow as posterior an inverted *gamma* distribution.

But since neither parameter is known, the MCMC algorithm will proceed in
an iterative way. From a starting value (given in `inits` parameter of
the model) for *σ*, values of *μ*<sub>*b*</sub> and *μ*<sub>*g*</sub>
are sampled in the conditional distribution given *σ*. Then a value of
*σ* is sample from its conditional distribution given the generated
values of *μ*<sub>*b*</sub> and *μ*<sub>*g*</sub>. This process is
repeated a large number of times.

This process is a **Gibbs sampling**, which is a specialized case of the
Metropolis algorithm. Gibbs sampling is relevant when the conditional
distributions of the parameters are explicitly given.

# 4. Analysis of convergence and autocorrelation

## Convergence

⇒ Check visually the convergence form the track of the 3 runs of MCMC:

``` r
plot(mcmc_birth_w)
```

![](birth_weigth_files/figure-markdown_github/unnamed-chunk-8-1.png)

Whatever the initial values of the parameters, the three chains seem to
sample in similar ranges. The chains overlap, and mix well. This is a
sign of convergence of the algorithm (convergence towards a
distribution, and not a value).

``` r
gelman.diag(mcmc_birth_w)
```

    ## Potential scale reduction factors:
    ## 
    ##         Point est. Upper C.I.
    ## mean[1]          1          1
    ## mean[2]          1          1
    ## sd               1          1
    ## 
    ## Multivariate psrf
    ## 
    ## 1

``` r
gelman.plot(mcmc_birth_w)
```

![](birth_weigth_files/figure-markdown_github/unnamed-chunk-10-1.png)

The variance reduction index is 1 for the 3 parameters on all the 5000
kept iterations. This index is defined as:

$\\sqrt{\\frac{var(total)}{var(wthin-chains)}}$

It means that between-chains variance is negligible in front of
within-chains variance, and so the three chains sample values from the
same posterior distribution).  
The index of variance reduction has been close to 1 since the first
iterations that have been kept. 3000 iterations burn-in time are enough
here to reach convergence.

## Autocorrelation

⇒ Check the autocorrelation of values inside the chains:

``` r
autocorr.plot(mcmc_birth_w)
```

![](birth_weigth_files/figure-markdown_github/unnamed-chunk-11-1.png)![](birth_weigth_files/figure-markdown_github/unnamed-chunk-11-2.png)![](birth_weigth_files/figure-markdown_github/unnamed-chunk-11-3.png)

Autocorrelation is weak or the three parameters. The autocorrelation is
computed for each parameter and each chain.

``` r
lapply(mcmc_birth_w, effectiveSize)
```

    ## [[1]]
    ##  mean[1]  mean[2]       sd 
    ## 3389.770 3084.136 1280.059 
    ## 
    ## [[2]]
    ##   mean[1]   mean[2]        sd 
    ## 3253.3262 3175.7315  941.8373 
    ## 
    ## [[3]]
    ##  mean[1]  mean[2]       sd 
    ## 3069.195 3104.783 1268.988

The efficiency size of the chains following the chains and parameter is
between 2000 and 3000, not too different from the number of iterations
per chain (5000).

``` r
raftery.diag(mcmc_birth_w)
```

    ## [[1]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                                
    ##          Burn-in  Total Lower bound  Dependence
    ##          (M)      (N)   (Nmin)       factor (I)
    ##  mean[1] 5        6078  3746         1.62      
    ##  mean[2] 5        5577  3746         1.49      
    ##  sd      8        9123  3746         2.44      
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
    ##  mean[1] 6        6406  3746         1.71      
    ##  mean[2] 6        6406  3746         1.71      
    ##  sd      12       12572 3746         3.36      
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
    ##  mean[1] 6        6637  3746         1.77      
    ##  mean[2] 5        6119  3746         1.63      
    ##  sd      9        9892  3746         2.64

Raftery’s diagnostic signals dependency factors close to 1, with a total
N close to the 5000 realized iterations. We can conclude that is not
necessary to perform more iterations.

## Number of sampled values

⇒ Check that are enough sampled values to predict correct posterior
distributions of the parameters (`cumuplot`):

``` r
cumuplot(mcmc_birth_w)
```

![](birth_weigth_files/figure-markdown_github/unnamed-chunk-14-1.png)![](birth_weigth_files/figure-markdown_github/unnamed-chunk-14-2.png)![](birth_weigth_files/figure-markdown_github/unnamed-chunk-14-3.png)

The estimates of means and quantiles at 2.5% and 97.5% seem stable after
2000 iterations for the 3 parameters. It means that are enough sampled
values to describe the posterior distribution (perhaps a few more would
be relevant for parameter *σ*). With this graph we cans also evaluate
convergence.

# 5. Exploit the results

⇒ Represent joint posterior distribution of the parameters (function
`pairs`).  
It is useful to concatenate the three chains to evaluate the *a
posteriori* distribution (`as.data.frame.(as.matrix(.))`).

``` r
mcmc_birth_w_tot <- as.data.frame(as.matrix(mcmc_birth_w))
pairs(mcmc_birth_w_tot)
```

![](birth_weigth_files/figure-markdown_github/unnamed-chunk-15-1.png)

Sampled values are a sample of the joint posterior distribution of the
parameters. So, the covariance between the samples values of the
parameters is an estimate of the real covariance between parameters. On
this example, the correlation between the sampled values seems very
weak, which shows absence of correlation between the parameters.

⇒ Compute a point estimate and a 95% - credibility interval of
*μ*<sub>*b*</sub> and *μ*<sub>*g*</sub>.

``` r
summary(mcmc_birth_w)
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
    ##           Mean     SD  Naive SE Time-series SE
    ## mean[1] 3728.1 5.3661 0.0438138       0.054531
    ## mean[2] 3379.3 6.6360 0.0541828       0.068581
    ## sd       799.9 0.1062 0.0008673       0.001816
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##           2.5%    25%    50%  75% 97.5%
    ## mean[1] 3717.7 3724.4 3728.0 3732  3739
    ## mean[2] 3366.4 3374.9 3379.2 3384  3392
    ## sd       799.6  799.9  799.9  800   800

Posterior mean and median are very close because the posterior
distributions are symmetric. 95% - credibility intervals can be obtained
from the 2.5% and 97.5% quantiles of the sampled values.

Girls: 3728.2 \[3718.0, 3739\]  
Boys: 3379.2 \[3366.6, 3392\]

⇒ Find out the above results “by hand”, without using the function
`summary` applied to MCMC. We have to compute the statistic on the
sampled values on the three chains together.

``` r
summary(mcmc_birth_w_tot)
```

    ##     mean[1]        mean[2]           sd       
    ##  Min.   :3707   Min.   :3354   Min.   :799.0  
    ##  1st Qu.:3724   1st Qu.:3375   1st Qu.:799.9  
    ##  Median :3728   Median :3379   Median :799.9  
    ##  Mean   :3728   Mean   :3379   Mean   :799.9  
    ##  3rd Qu.:3732   3rd Qu.:3384   3rd Qu.:800.0  
    ##  Max.   :3749   Max.   :3405   Max.   :800.0

``` r
quantile(mcmc_birth_w_tot[, "mean[2]"], prob = c(0.025, 0.975))
```

    ##     2.5%    97.5% 
    ## 3366.422 3392.326

``` r
quantile(mcmc_birth_w_tot[, "mean[1]"], prob = c(0.025, 0.975))
```

    ##     2.5%    97.5% 
    ## 3717.708 3738.621

``` r
quantile(mcmc_birth_w_tot[, "sd"], prob = c(0.025, 0.975))
```

    ##     2.5%    97.5% 
    ## 799.6003 799.9972

``` r
apply(mcmc_birth_w_tot, 2, quantile, prob = c(0.025, 0.975))
```

    ##        mean[1]  mean[2]       sd
    ## 2.5%  3717.708 3366.422 799.6003
    ## 97.5% 3738.621 3392.326 799.9972

⇒ Compute the probability that the mean weight of boys is larger than
the mean weight of girls.

An estimate of this probability can be calculated as the proportion,
among the couples of parameters sampled in the posterior distribution,
of the couples wheres the boys’ mean weight is larger than the girls’
mean weight.  
Graphically, it can be seen as the proportion of iterations for wichh
the couple (*μ*<sub>*b*</sub>, *μ*<sub>*g*</sub>) is over the diagonal.

``` r
sum(mcmc_birth_w_tot[, "mean[1]"] > mcmc_birth_w_tot[, "mean[2]"])/nrow(mcmc_birth_w_tot)
```

    ## [1] 1

The probability that the boys’ mean weight is larger than the girls’
mean weight is very high.

# 6. Evolution of the model

Children weight is positive. This constraint is not in the model, since
a normal distribution also cover negative values. A solutions is to
truncate the distribution of weights at 0. In **JAGS**, the truncature
operator is `T(a, b) : a` is the low troncature, `b`, the high
truncature (leaving it empty if we do not want truncature on a side).

⇒ Build a new model with this constraint, and compare the results with
those of the first model.

``` r
desc_model_birth_w2 <- 
"model {
  for (i in 1:Nchild) {
    w[i] ~ dnorm(mu[i], tau)T(0, )
    mu[i] <- mean[sex[i]]
  }
  mean[1] ~ dunif(2500, 5000)
  mean[2] ~ dunif(2500, 5000)
  tau <- 1/(sd^2)
  sd ~ dunif(200, 800)
}"
```

``` r
model_birth_w2 <- jags.model(textConnection(desc_model_birth_w2), 
                             data = data_birth_w, 
                             inits = init_birth_w, 
                             n.chains = 3)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 48
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 109
    ## 
    ## Initializing model

``` r
update(model_birth_w2, 3000)
mcmc_birth_w2 <- coda.samples(model_birth_w2, variable.names = c("mean[1]", "mean[2]", "sd"), n.iter = 5000)
```

``` r
summary(mcmc_birth_w2)
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
    ##           Mean     SD Naive SE Time-series SE
    ## mean[1] 3727.9  99.13   0.8094         1.0111
    ## mean[2] 3379.6 123.29   1.0067         1.2989
    ## sd       530.7  56.51   0.4614         0.6342
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##           2.5%    25%    50%    75%  97.5%
    ## mean[1] 3534.1 3662.0 3727.4 3793.4 3922.9
    ## mean[2] 3131.0 3298.4 3380.2 3461.8 3619.3
    ## sd       433.5  490.5  525.9  565.6  653.1

The results do not change. Indeed, considering the standard deviation of
the first model, the probability to get negative weight is quasi-null.
Likewise, we could have constrained the prior distribution on
*μ*<sub>*b*</sub> and *μ*<sub>*g*</sub> to be positive, but it would not
have any impact on the results in this example, since the *a posteriori*
distributions are totally negatives values.
