# Tutorial for the R package *bayou*
The purpose of *bayou* is to fit Bayesian models of adaptive evolution to phylogenetic comparative data. Specifically, *bayou* provides a flexible framework for fitting multi-optima Ornstein-Uhlenbeck models to phylogenetic comparative data. This tutorial demonstrates some of the options for running *bayou*.

## Reversible-jump MCMC over regime placement
In this example, we will fit a reversible-jump MCMC model to an included dataset (the Chelonia dataset of Jaffe et al. 2012). To start, we will specify a model in which no parameters are fixed. This will estimate the posterior of shift number, location and magnitude as well as all other parameters.

We begin by loading the package and data. We will also assume a constant standard error across all taxa. This can instead be a named vector with species-specific measurement error.

```r
require(bayou)
data(chelonia)
tree <- chelonia$phy
dat <- chelonia$dat
SE <- 0.05
```



### Defining a prior function
We now need to define prior function to set up our model. This can be the trickiest part. We will set half-cauchy priors for alpha and sigma^2, a normal prior for theta, a conditional Poisson for the number of shifts, and "dsb" controls how many shifts can be per branch (either 0, 1 or Inf) and the probability of a shift being on that branch. Since we have set bmax = 1 and prob = 1; we are specifying a model with a maximum of 1 shift per branch, and equal probability across all branches.

```r
prior <- make.prior(tree, dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
    dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), param = list(dalpha = list(scale = 1), 
    dsig2 = list(scale = 1), dk = list(lambda = 15, kmax = 200), dsb = list(bmax = 1, 
        prob = 1), dtheta = list(mean = mean(dat), sd = 2)))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

The figure produced gives a rough visual of the chosen prior distributions. 

### Running the MCMC
Now we are going to run the mcmc. We are going to output files to our working directory. If you want to specify another director, replace "getwd()" with the path of the directory. By default, *bayou outputs to the R temporary directory. We will run a relatively short chain of only 10,000 generations.

```r
par(mfrow = c(2, 3))
fit1 <- bayou.mcmc(tree, dat, SE = SE, model = "OU", prior, ngen = 10000, new.dir = getwd(), 
    plot.freq = 2000, ticker.freq = 1000)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1000			-130.43			-97.46			26.27			0.32			15			0.38			0.02			1			0.38			1			0.48			0.75			0.92			0.6			0.83			
```

```
## 2000			-116.05			-84.44			17.07			0.27			13			0.39			0.02			0.95			0.46			0.91			0.41			0.75			0.95			0.57			0.68			
## 3000			-111.7			-74.28			12.58			0.21			11			0.37			0.02			0.94			0.4			0.92			1			0.42			0.68			0.95			0.5			0.44			
```

```
## 4000			-109.56			-64.07			14.15			0.28			9			0.35			0.02			0.87			0.37			0.89			1			0.4			0.64			0.96			0.44			0.35			
## 5000			-113.35			-64.33			14.83			0.24			9			0.34			0.02			0.89			0.37			0.82			1			0.39			0.62			0.95			0.45			0.32			
```

```
## 6000			-119.61			-73.4			9.41			0.29			11			0.34			0.02			0.9			0.39			0.8			1			0.38			0.6			0.96			0.48			0.32			
## 7000			-111.14			-63.5			10.26			0.22			9			0.34			0.02			0.88			0.41			0.79			1			0.37			0.59			0.96			0.47			0.33			
```

```
## 8000			-105.83			-90.55			14.71			0.21			14			0.34			0.02			0.88			0.4			0.77			1			0.37			0.59			0.96			0.45			0.31			
## 9000			-114.33			-94.41			12.41			0.25			15			0.34			0.02			0.88			0.4			0.8			1			0.37			0.59			0.96			0.46			0.32			
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```
## 10000			-110.92			-74.86			13.86			0.25			11			0.35			0.02			0.89			0.41			0.79			1			0.36			0.59			0.96			0.46			0.33			
```


Most of the output is saved in a file, only an overview of the run is saved in *fit1*. 

```r
fit1
```

```
## bayou modelfit
## OU parameterization
## 
## Results are stored in directory
## /home/josef/repos/bayou/UQYZHRKAJW/bayou.* 
## To load results, use 'load.bayou(bayouFit)'
## 
## 10000  generations were run with the following acceptance probabilities:
##   .alpha  birth.k D0.slide D1.slide  death.k R1.slide    .sig2   .theta 
##     0.35     0.02     0.89     0.41     0.79     1.00     0.36     0.59 
## U0.slide U1.slide U2.slide 
##     0.96     0.46     0.33 
##  Total number of proposals of each type:
##   .alpha  birth.k D0.slide D1.slide  death.k R1.slide    .sig2   .theta 
##     1812     4523      171      280      112        1      853     1773 
## U0.slide U1.slide U2.slide 
##      178      157      140
```


We can load the actual chains by running the following code:

```r
chain <- load.bayou(fit1, save.Rdata = FALSE, cleanup = FALSE)
chain <- set.burnin(chain, 0.3)
```


We can return a summary of our MCMC results by summarizing the chain. Notice the very small effective sample sizes for parameters for only 10,000 generations, need to get more like 100 for each by running the MCMC for more generations.

```r
out <- summary(chain)
```

```
## bayou MCMC chain: 10000 generations
## 1000 samples, first 300 samples discarded as burnin
## 
## 
## Summary statistics for parameters:
##                 Mean       SD  Naive SE Time-series SE Effective Size
## lnL       -113.32076  4.25131 0.1605696       1.143409         13.824
## prior      -75.74338 13.59449 0.5134568       5.658736          5.771
## alpha        0.06657  0.02644 0.0009986       0.008482          9.717
## sig2         0.03117  0.01075 0.0004060       0.003640          8.721
## k           11.27817  2.62808 0.0992611       1.096298          5.747
## ntheta      12.27817  2.62808 0.0992611       1.096298          5.747
## root         3.63773  0.13500 0.0050989       0.057955          5.426
## all theta    3.85801  1.02167        NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##         pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 372 0.9857               4.388           0.010529      1.56470
## 408 0.9287               4.969           0.009236      1.84564
## 411 0.6748               3.228           0.003935      1.92874
## 45  0.4237               4.728           0.021360      0.23616
## 264 0.3381               3.000           0.033473      0.75655
## 238 0.2739               5.619           0.039559      0.56377
## 409 0.2568               3.694           0.057911      0.70005
## 51  0.2482               2.546           0.047140      2.89243
## 345 0.2068               4.040           0.048312     50.00154
## 43  0.1926               4.941           0.027902      0.16253
## 177 0.1755               3.208           0.092757      1.45866
## 275 0.1712               3.512           0.050490      1.24727
## 232 0.1698               5.519           0.083403      3.05968
## 439 0.1698               3.105           0.045447      0.40109
## 91  0.1583               4.086           0.045433      0.29825
## 307 0.1583               1.846           0.096121      0.51436
## 183 0.1569               2.731           0.030756      0.03843
## 441 0.1569               2.985           0.031698      1.00535
## 42  0.1526               5.031           0.042362      1.07982
## 240 0.1498               4.593           0.055112      0.55352
## 241 0.1498               4.166           0.024298      0.47713
## 191 0.1469               5.016           0.069491      0.40012
## 252 0.1455               4.294           0.018131      1.11714
## 412 0.1455               3.920           0.013670      0.81304
## 110 0.1384               1.858           0.057153      1.16477
## 19  0.1312               4.333           0.029933      0.01383
## 346 0.1284               4.900           0.074465      5.16316
## 415 0.1270               2.952           0.063830      2.83099
## 115 0.1170               4.305           0.044354     15.05473
## 224 0.1098               3.047           0.061372      1.36024
## 206 0.1070               2.694           0.022599      0.21588
## 407 0.1056               3.221           0.013790      0.56280
## 99  0.1041               3.535           0.043017      0.29114
```


We can visualize the chains by viewing traces for each parameter using the plotting utilities of the R package coda:

```r
plot(chain)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-71.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-72.png) 


We can view where there are shifts of high probability.

```r
par(mfrow = c(1, 1))
plotSimmap.mcmc(tree, chain, burnin = 0.3, type = "circles", fsize = 0.4, pts = FALSE)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


And we can view the density of phenotypic optima and location of highly supported shifts on the phenogram. Here we show all shifts with posterior probabilities greater than *pp.cutoff = 0.3*. 

```r
phenogram.density(tree, dat, chain = chain, burnin = 0.3, pp.cutoff = 0.3)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


### Diagnosing Convergence
It is useful to run two chains with independent starting positions. By default, *bayou* simulates starting parameters from the prior distribution. We will re-run the analysis here to obtain an independent chain.

```r
fit2 <- bayou.mcmc(tree, dat, SE = SE, model = "OU", prior, ngen = 10000, new.dir = getwd(), 
    plot.freq = NULL, ticker.freq = 1000)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1000			-142.72			-90.91			39.65			0.45			14			0.58			0.03			0.88			0.65			1			0.36			0.87			0.9			0.89			0.91			
## 2000			-131.82			-76.25			18.25			0.24			11			0.52			0.03			0.92			0.48			1			0.32			0.87			0.92			0.86			0.75			
## 3000			-131.58			-65.27			27.21			0.43			9			0.48			0.03			0.9			0.45			1			0.31			0.82			0.93			0.8			0.59			
## 4000			-125.7			-68.82			13.25			0.32			10			0.45			0.03			0.93			0.46			0.95			0.32			0.77			0.93			0.72			0.49			
## 5000			-113.82			-104.88			10.77			0.27			17			0.43			0.03			0.91			0.42			0.94			0.31			0.75			0.93			0.71			0.46			
## 6000			-114.99			-89.67			15.45			0.26			14			0.41			0.03			0.89			0.4			0.92			0.3			0.73			0.93			0.67			0.37			
## 7000			-116.72			-71.39			12.74			0.31			10			0.4			0.03			0.82			0.38			0.9			0.31			0.71			0.93			0.57			0.34			
## 8000			-103.55			-79.2			8.52			0.2			12			0.38			0.02			0.83			0.4			0.87			0.32			0.68			0.94			0.54			0.32			
## 9000			-99.95			-68.66			2.7			0.16			10			0.38			0.02			0.85			0.4			0.82			0.32			0.65			0.92			0.51			0.31			
## 10000			-104.93			-58.88			6.63			0.21			8			0.37			0.02			0.85			0.38			0.8			0.32			0.62			0.93			0.44			0.28			
```

```r
chain2 <- load.bayou(fit2, save.Rdata = FALSE, cleanup = FALSE)
chain2 <- set.burnin(chain2, 0.3)
```


Now we can compare chains to see if the parameters have converged using Gelman and Rubin's R statistic. Values close to 1 indicate chains have converged (these two chains will not have converged in only 10,000 generations).

```r
RlnL <- gelman.R("lnL", chain1 = chain, chain2 = chain2, plot = TRUE, type = "n", 
    ylim = c(0.9, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-111.png) 

```r
Ralpha <- gelman.R("alpha", chain1 = chain, chain2 = chain2, plot = TRUE, type = "n", 
    ylim = c(0.9, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-112.png) 

```r
Rsig2 <- gelman.R("sig2", chain1 = chain, chain2 = chain2, plot = TRUE, type = "n", 
    ylim = c(0.9, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-113.png) 


Of particular interest to inference are the branch posterior probabilities, we can plot these to see if the two chains converged on similar answers. If the runs have converged, points should fall along the y=x line.

```r
L1 <- Lposterior(chain, tree, burnin = 0.3)
L2 <- Lposterior(chain2, tree, burnin = 0.3)
plot(L1$pp, L2$pp, xlim = c(0, 1), ylim = c(0, 1), xlab = "Chain 1", ylab = "Chain 2")
curve(1 * x, add = TRUE, lty = 2)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 


We can also combine chains:

```r
chains <- combine.chains(chain, chain2, burnin.prop = 0.3)
chains <- set.burnin(chains, 0)
```


### Estimating marginal likelihoods
To compare models, we want to estimate the marginal likelihood, which is performed using the stepping stone algorithm in *bayou*. The method first estimates a reference function by fitting a series of curves to the posterior of an MCMC chain. The function *steppingstone* will output a graphic showing the best-fitting density functions overlaying the posterior distribution. Then the stepping stone MCMC's are run for every value in the vector *Bk*, corresponding to each value of the power posterior function ranging from the reference function (*Bk = 0*) to the posterior (*Bk = 1*). 
To speed things up, you can specify multiple cores. Default is set to 2 cores, but if you have more you can use them.

```r
ss <- steppingstone(Bk = seq(0, 1, length.out = 5), chains, tree, dat, SE = SE, 
    prior = prior, new.dir = getwd(), ngen = 10000, cores = 2)
```

```
## Making power posterior function from provided mcmc chain...
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

```
## Running mcmc chains...
## Loading mcmc chains...
```

```r
ss
```

```
## Stepping stone estimation of marginal likelihood
## Marginal Likelihood:
## [1] -144.7
## A total of 5 power posteriors were run along the sequence: 0		0.25		0.5		0.75		1
## lnr_k -39.78 -37.25 -35.63 -32
```

```r
ss <- set.burnin(ss, 0.3)
```


Plotting the stepping stone object will give you an idea of convergence for each of the MCMC's by providing a trace of the lnL, ln prior, and reference function as well as the estimated marginal likelihood at each step.

```r
plot(ss)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 


## Fitting models with fixed parameters
*bayou* also allows the user to fit models with fixed parameters, including non-reversible jump models. A helpful utility for setting up fixed hypotheses is the function *identify.branches*. First, we set up some starting parameters.

```r
startpar <- list(alpha = 1, sig2 = 1, k = 2, ntheta = 3, theta = c(4, 5, 6))
```


Then we use *identify.branches()* to specify the location of shifts by clicking on them. f

```r
fixed.hypothesis <- identify.branches(tree, startpar$k)
```






```r
startpar$sb <- fixed.hypothesis$sb
startpar$t2 <- 2:startpar$ntheta
startpar$loc <- fixed.hypothesis$loc
plotBayoupars(startpar, tree, fsize = 0.5)
```

```
## no colors provided. using the following legend:
##        1        2        3 
##  "black"    "red" "green3"
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 


This is a good time to explain how *bayou* stores parameters. *bayou* takes parameter values as a list of values corresponding to the model. For the standard OU model, these parameters are:

```r
startpar
```

```
## $alpha
## [1] 1
## 
## $sig2
## [1] 1
## 
## $k
## [1] 2
## 
## $ntheta
## [1] 3
## 
## $theta
## [1] 4 5 6
## 
## $sb
## [1] 411 409
## 
## $t2
## [1] 2 3
## 
## $loc
## [1] 24  8
```

*k* is the number of shifts and *ntheta* is the number of phenotypic optima. *theta* is a vector of optima values equal in length to *ntheta*, with *theta*[1] corresponding to the root value and optimum. 
*sb* specifies the branch locations of shifts, and is equal in length to *k*. Note that because this vector specifies branches numbers, *bayou* always works with postorder trees. The vectors *loc* and *t2* are equal in length to *sb*, and specify the location of the shift on the branch (in distance from the starting point of that branch) and the identity of the optima *after the shift* (which corresponds to the *i*th element of the *theta* vector). 
While convergent optima (shifts that lead to the same *theta*), cannot be used in the reversible-jump algorithm, we can fit convergent models using fixed shift locations by specifying multiple shifts to the same value of *t2*. For example:


```r
converge_pars <- startpar
converge_pars$k <- 3
converge_pars$ntheta <- 3
converge_pars$sb <- c(converge_pars$sb, 49)
converge_pars$loc <- 0
converge_pars$t2 <- 3
converge_pars
```

```
## $alpha
## [1] 1
## 
## $sig2
## [1] 1
## 
## $k
## [1] 3
## 
## $ntheta
## [1] 3
## 
## $theta
## [1] 4 5 6
## 
## $sb
## [1] 411 409  49
## 
## $t2
## [1] 3
## 
## $loc
## [1] 0
```

```r
plotBayoupars(converge_pars, tree, fsize = 0.5)
```

```
## Error: 'names' attribute [453] must be the same length as the vector [451]
```


Now that we know how to specify fixed models, lets set up the prior function:

```r
prior.fixed <- make.prior(tree, dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
    dsb = "fixed", dk = "fixed", dtheta = "dnorm", dloc = "dloc"), param = list(dalpha = list(scale = 1), 
    dsig2 = list(scale = 1), dtheta = list(mean = mean(dat), sd = 2)), fixed = fixed.hypothesis)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 


A few things to note. First, we have set *dsb* and *dk* to "fixed". We could have also set *dloc*, *dalpha* and/or *dsig2* to "fixed". Here, we will let the location of the shift on the branch vary (*dloc* puts a uniform prior on the location of the shift on the branch). Second, we are passing our model through the *fixed.hypothesis* object we created using *identify.branches()*. For example, since simple Brownian motion is a special case of multi-optima OU models, we can specify through the the following prior:

```r
prior.BM <- make.prior(tree, dists = list(dalpha = "fixed", dsig2 = "dhalfcauchy", 
    dsb = "fixed", dk = "fixed", dloc = "fixed"), param = list(dsig2 = list(scale = 1), 
    dtheta = list(mean = mean(dat), sd = 2)), fixed = list(alpha = 0, k = 0, 
    sb = numeric(0), loc = numeric(0), t2 = numeric(0)))
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 


Now we can run our MCMC chains using our priors. Right now, you have to specify the starting parameters. This will soon be unnecessary...

```r
fit.fixed <- bayou.mcmc(tree, dat, SE = SE, model = "OU", prior = prior.fixed, 
    startpar = startpar, ngen = 10000, new.dir = getwd(), plot.freq = NULL)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 1000			-141.25			-5.79			60.42			0.63			2			0.38			1			0.32			0.33			1			
## 2000			-143.16			-5.89			37.1			0.45			2			0.41			0.98			0.3			0.38			1			
## 3000			-140.39			-5.77			31.25			0.46			2			0.43			0.98			0.33			0.38			0.99			
## 4000			-140.44			-5.75			33.39			0.48			2			0.44			0.98			0.32			0.39			0.99			
## 5000			-141.67			-5.79			45.48			0.49			2			0.45			0.98			0.32			0.39			1			
## 6000			-139.69			-5.76			33.15			0.4			2			0.45			0.98			0.33			0.4			0.99			
## 7000			-144.43			-5.99			33.46			0.4			2			0.45			0.97			0.32			0.4			0.99			
## 8000			-143.7			-5.76			71.75			0.71			2			0.45			0.97			0.32			0.4			0.99			
## 9000			-141.88			-5.78			26.45			0.38			2			0.46			0.98			0.32			0.4			0.99			
## 10000			-138.48			-5.75			38.67			0.48			2			0.46			0.98			0.33			0.41			1			
```






```r
chain.fixed <- load.bayou(fit.fixed, save.Rdata = FALSE, cleanup = FALSE)
chain.fixed <- set.burnin(chain.fixed, 0.3)
plot(chain.fixed)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-261.png) ![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-262.png) 

```r
out.fixed <- summary(chain.fixed)
```

```
## bayou MCMC chain: 10000 generations
## 1000 samples, first 300 samples discarded as burnin
## 
## 
## Summary statistics for parameters:
##                 Mean       SD  Naive SE Time-series SE Effective Size
## lnL       -140.86585 1.754982 0.0662848      0.1476304         141.32
## prior       -5.79783 0.050636 0.0019125      0.0043990         132.50
## alpha        0.01939 0.004967 0.0001876      0.0004995          98.88
## sig2         0.01814 0.002846 0.0001075      0.0003536          64.79
## k            2.00000 0.000000 0.0000000      0.0000000           0.00
## ntheta       3.00000 0.000000 0.0000000      0.0000000           0.00
## root         3.76423 0.175547 0.0066303      0.0107514         266.60
## all theta    3.65028 0.354228        NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##     pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 409  1               3.635            0.01202        1.358
## 411  1               3.551            0.01778        2.176
```

```r
phenogram.density(tree, dat, chain = chain.fixed, burnin = 0.3, pp.cutoff = 0.5)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-263.png) 

```r
plotSimmap.mcmc(tree, chain.fixed, burnin = 0.3, type = "circles", pts = FALSE, 
    fsize = 0.5)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-264.png) 





As before, we can estimate the marginal likelihood:

```r
ss.fixed <- steppingstone(Bk = seq(0, 1, length.out = 5), chain.fixed, tree, 
    dat, SE = 0, startpar = startpar, prior = prior.fixed, ngen = 10000, cores = 5)
```

```
## Making power posterior function from provided mcmc chain...
```

```
## Running mcmc chains...
## Loading mcmc chains...
```

```r
ss.fixed <- set.burnin(ss.fixed, 0.3)
ss.fixed$lnr
```

```
## [1] -162.8
```

```r
plot(ss.fixed)
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-281.png) ![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-282.png) 


And calculate Bayes Factors in support of the Reversible-Jump model vs. the fixed model. 

```r
2 * (ss$lnr - ss.fixed$lnr)
```

```
## [1] 36.28
```


## Converting to and from *OUwie* format

*bayou* includes some utilities to switch between *OUwie* and *bayou* formatted models. 

```r
require(OUwie)
OUwieData <- bayou2OUwie(startpar, tree, dat)
OUwieRes <- OUwie(OUwieData$tree, OUwieData$dat, model = "OUM")
```

```
## Initializing... 
## Finished. Begin thorough search... 
## Finished. Summarizing results.
```


Here's a comparison of *bayou* to *OUwie* models. NB: At the moment they won't be that similar because the root state is treated differently in both, and I have not yet included the ability to use measurement error in *bayou2OUwie()*, which is on my to do list...

```r
OUwieRes
```

```
## 
## Fit
##    -lnL   AIC  AICc model ntax
##  -148.1 306.1 306.4   OUM  226
## 
## 
## Rates
##                1       2       3
## alpha    0.03233 0.03233 0.03233
## sigma.sq 0.02785 0.02785 0.02785
## 
## Optima
##               1      2      3
## estimate 3.7433 3.4769 3.5108
## se       0.1221 0.3874 0.1312
## 
## Arrived at a reliable solution
```

```r
out.fixed$statistics
```

```
##                 Mean       SD  Naive SE Time-series SE Effective Size
## lnL       -140.86585 1.754982 0.0662848      0.1476304         141.32
## prior       -5.79783 0.050636 0.0019125      0.0043990         132.50
## alpha        0.01939 0.004967 0.0001876      0.0004995          98.88
## sig2         0.01814 0.002846 0.0001075      0.0003536          64.79
## k            2.00000 0.000000 0.0000000      0.0000000           0.00
## ntheta       3.00000 0.000000 0.0000000      0.0000000           0.00
## root         3.76423 0.175547 0.0066303      0.0107514         266.60
## all theta    3.65028 0.354228        NA             NA             NA
```


Of course, we can also go the other direction, and convert back to *bayou* format:

```r
back_pars <- OUwie2bayou(OUwieData$tree, OUwieData$dat)
plotBayoupars(back_pars, tree)
```

```
## no colors provided. using the following legend:
##        1        2        3 
##  "black"    "red" "green3"
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32.png) 


## Some other useful things
### Some plotting features
It's easy to convert *bayou* formatted models into trees that can be visualized using *phytools*' *plotSimmap* or *phenogram* functions. For example, let's pull a few samples from the posterior out of our original chain:

```r
samp_pars <- pull.pars(500, chain, model = "OU")
samp_pars
```

```
## $alpha
## [1] 0.04674
## 
## $sig2
## [1] 0.02289
## 
## $k
## [1] 9
## 
## $ntheta
## [1] 10
## 
## $theta
##  [1] 3.697 4.904 4.388 3.102 5.515 2.776 1.660 3.458 4.921 2.415
## 
## $sb
## [1] 409 372 407  43 230 108 344  40 210
## 
## $loc
## [1] 5.877 1.935 5.336 5.832 5.777 0.113 1.057 2.140 4.429
## 
## $t2
## [1]  2  3  4  5  6  7  8  9 10
```


Now let's convert it to *phytools*' simmap format:

```r
sm_tree <- pars2simmap(samp_pars, tree)
plotSimmap(sm_tree$tree, col = sm_tree$col)
phenogram(sm_tree$tree, dat, col = sm_tree$col, ftype = "off")
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-341.png) ![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-342.png) 


We can also make some regime plots, which show the location of each optimum and the expected stationary variance.

```r
plot(c(0, 210), c(2, 6), type = "n", xlab = "time", ylab = "phenotype")
regime.plot(samp_pars, sm_tree$tree, type = "density", cols = sm_tree$col)
phenogram(sm_tree$tree, dat, col = sm_tree$col, ftype = "off", add = TRUE)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35.png) 


....Put in how to do OU ancestral state reconstructions here....


## How to specify some other models
Rather than restricting our interest to only models with a single shift allowed per branch, we can fit a number of models that allow multiple shifts per branch, have unequal probabilities, etc. These are all specified with the prior function. Here are some examples:

A model with an arbitrary number of shifts allowed per branch, and branches chosen with probability proportional to their length:

```r
tree <- reorder(tree, "postorder")
dat <- dat[tree$tip.label]
prior <- make.prior(tree, dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
    dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), param = list(dalpha = list(scale = 1), 
    dsig2 = list(scale = 1), dk = list(lambda = 15, kmax = 200), dsb = list(bmax = Inf, 
        prob = tree$edge.length), dtheta = list(mean = mean(dat), sd = 2)))
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-36.png) 


A model that disallows shifts on terminal branches:

```r
terminal.branches <- as.numeric(tree$edge[, 2] < length(tree$tip.label) + 1)
prior <- make.prior(tree, dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
    dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), param = list(dalpha = list(scale = 1), 
    dsig2 = list(scale = 1), dk = list(lambda = 15, kmax = 200), dsb = list(bmax = terminal.branches, 
        prob = 1), dtheta = list(mean = mean(dat), sd = 2)))
```

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-37.png) 



...Add QG, OUrepar models here...


