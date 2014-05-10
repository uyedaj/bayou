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
## 1000			-134.58			-63.84			30.35			0.4			9			0.49			0.03			1			0.6			1			0.38			0.7			1			0.47			0.36			
```

```
## 2000			-125.64			-64.41			17.62			0.28			9			0.42			0.02			0.93			0.52			1			0.37			0.71			1			0.45			0.37			
## 3000			-112.33			-70			11.11			0.23			10			0.4			0.03			0.92			0.46			0.95			0.33			0.71			0.96			0.43			0.33			
```

```
## 4000			-104.8			-116.63			5.9			0.2			19			0.38			0.03			0.93			0.46			0.88			0.33			0.68			0.93			0.43			0.3			
## 5000			-111.65			-94.73			9.49			0.21			15			0.37			0.03			0.92			0.45			0.86			0.31			0.67			0.93			0.42			0.29			
```

```
## 6000			-109.62			-100.71			17.25			0.23			16			0.37			0.03			0.91			0.43			0.84			0.31			0.65			0.93			0.43			0.29			
## 7000			-115.51			-87.4			11.29			0.2			13			0.35			0.03			0.9			0.46			0.82			0.31			0.64			0.92			0.43			0.33			
```

```
## 8000			-111.41			-79.88			15.76			0.25			12			0.36			0.03			0.92			0.46			0.79			0.31			0.63			0.9			0.41			0.32			
## 9000			-107.88			-79.45			8.46			0.27			12			0.35			0.02			0.9			0.45			0.77			0.32			0.62			0.9			0.4			0.3			
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```
## 10000			-103.63			-104.43			10.15			0.21			17			0.35			0.02			0.89			0.44			0.75			0.32			0.62			0.9			0.41			0.29			
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
## /home/josef/repos/bayou/WREQXOFCZA/bayou.* 
## To load results, use 'load.bayou(bayouFit)'
## 
## 10000  generations were run with the following acceptance probabilities:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     0.35     0.02     0.89     0.44     0.75     0.32     0.62     0.90 
## U1.slide U2.slide 
##     0.41     0.29 
##  Total number of proposals of each type:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     1821     4424      206      245      131      909     1797      182 
## U1.slide U2.slide 
##      139      146
```


We can load the actual chains by running the following code:

```r
chain <- load.bayou(fit1, save.Rdata = FALSE, cleanup = TRUE)
```

```
## deleting temporary directory /home/josef/repos/bayou/WREQXOFCZA/
```

```r
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
## lnL       -108.99424  4.51500 0.1705294       1.369693         10.866
## prior      -89.93537 14.11418 0.5330852       6.218478          5.152
## alpha        0.07312  0.02337 0.0008826       0.008971          6.785
## sig2         0.03145  0.00715 0.0002701       0.001870         14.622
## k           13.91726  2.70323 0.1020994       1.055387          6.561
## ntheta      14.91726  2.70323 0.1020994       1.055387          6.561
## root         3.69236  0.07101 0.0026821       0.022876          9.636
## all theta    3.87570  1.17245        NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##         pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 408 0.8260               4.869           0.009390      1.92984
## 372 0.6177               4.236           0.011101      1.91403
## 232 0.4779               5.295           0.033253      2.45346
## 346 0.4779               5.293           0.026893      4.75256
## 42  0.4736               5.522           0.028520      1.38170
## 368 0.3780               2.011           0.042089     20.23833
## 333 0.3652               3.669           0.053000      1.15205
## 411 0.3652               3.180           0.002330      1.78639
## 406 0.3452               3.189           0.002652      0.39133
## 439 0.3124               3.151           0.020047      0.34204
## 365 0.2967               4.282           0.016230      2.99965
## 441 0.2882               2.962           0.015469      1.25136
## 66  0.2853               5.105           0.048520      1.30531
## 43  0.2582               4.810           0.018755      0.32188
## 437 0.2397               3.037           0.025373      1.92896
## 45  0.2382               4.888           0.014709      0.26125
## 264 0.2354               4.825           0.098697      2.07276
## 370 0.2211               2.088           0.028785     75.40673
## 243 0.2126               4.306           0.093597      5.25943
## 407 0.2126               3.243           0.001712      1.27622
## 300 0.2054               4.171           0.030948      0.28542
## 65  0.1997               4.738           0.022647      5.18765
## 269 0.1997               3.594           0.060601      2.22674
## 336 0.1969               4.342           0.070638      0.71842
## 272 0.1940               1.117           0.108825      2.51353
## 271 0.1912               3.484           0.037183      0.33508
## 273 0.1897               4.387           0.056632      3.86264
## 252 0.1769               4.086           0.015752      1.14261
## 237 0.1526               2.815           0.093503      0.88509
## 77  0.1484               3.183           0.081758      0.04389
## 263 0.1455               2.999           0.057618      0.17478
## 367 0.1441               1.245           0.044919      7.12252
## 306 0.1369               2.158           0.020246      1.16421
## 101 0.1341               5.707           0.085489      0.02325
## 389 0.1298               3.796           0.024842      0.32272
## 178 0.1255               5.841           0.078406      0.93575
## 242 0.1227               3.864           0.033218      7.55929
## 37  0.1213               6.139           0.046494      0.13356
## 387 0.1127               3.689           0.046900      0.10509
## 36  0.1098               3.536           0.107129      0.25252
## 292 0.1084               5.056           0.070089     37.60859
## 74  0.1041               3.779           0.078586      0.01571
```


We can visualize the chains by viewing traces for each parameter using the plotting utilities of the R package coda:

```r
plot(chain)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-71.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-72.png) 


We can view where there are shifts of high probability.

```r
par(mfrow = c(1, 1))
plotSimmap.mcmc(tree, chain, burnin = 0.3, circle = TRUE, fsize = 0.4)
```

```
## Error: object 'L' not found
```


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
## 1000			-140.19			-58.64			26.48			0.39			8			0.61			0.02			0.84			0.53			1			0.49			0.81			1			0.8			0.69			
## 2000			-125.21			-53.61			29.39			0.31			7			0.51			0.01			0.86			0.5			0.88			0.42			0.77			1			0.71			0.5			
## 3000			-129.84			-50.45			44.81			0.47			6			0.46			0.01			0.89			0.4			0.92			0.39			0.74			1			0.66			0.44			
## 4000			-127.17			-74.58			31.69			0.43			11			0.43			0.02			0.89			0.41			0.87			0.38			0.72			0.98			0.65			0.42			
## 5000			-133.52			-72.64			29.82			0.39			10			0.41			0.02			0.89			0.4			0.87			0.36			0.7			0.99			0.63			0.42			
## 6000			-126.01			-86.75			27.92			0.26			13			0.4			0.02			0.86			0.38			0.88			0.34			0.71			0.96			0.59			0.42			
## 7000			-122			-81.27			21.49			0.31			12			0.39			0.02			0.86			0.38			0.87			0.34			0.72			0.96			0.56			0.4			
## 8000			-115.53			-79.76			14.83			0.3			12			0.39			0.02			0.86			0.37			0.85			0.33			0.72			0.97			0.53			0.38			
## 9000			-108.96			-73.57			10.96			0.27			11			0.38			0.02			0.88			0.38			0.82			0.33			0.7			0.97			0.54			0.37			
## 10000			-111.15			-90.69			13.04			0.27			14			0.37			0.02			0.88			0.38			0.83			0.33			0.69			0.97			0.49			0.35			
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
## Loading required package: foreach
## foreach: simple, scalable parallel programming from Revolution Analytics
## Use Revolution R for scalability, fault tolerance and more.
## http://www.revolutionanalytics.com
## Loading required package: doMC
## Loading required package: iterators
## Loading required package: parallel
```

```
## Making power posterior function from provided mcmc chain...
```

```
## Loading required package: fitdistrplus
## Loading required package: survival
## Loading required package: splines
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
## [1] -146.3
## A total of 5 power posteriors were run along the sequence: 0		0.25		0.5		0.75		1
## lnr_k -39.28 -37.2 -36.83 -32.95
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
## 1000			-142.11			-5.75			47.98			0.45			2			0.39			0.97			0.38			0.34			0.99			
## 2000			-140.24			-5.76			47.25			0.48			2			0.4			0.98			0.39			0.38			0.99			
## 3000			-140.45			-5.74			47.45			0.56			2			0.42			0.98			0.35			0.4			1			
## 4000			-139.51			-5.78			44.1			0.47			2			0.44			0.98			0.35			0.4			0.99			
## 5000			-140.32			-5.78			34.62			0.49			2			0.43			0.98			0.35			0.39			0.99			
## 6000			-139.69			-5.77			39.12			0.49			2			0.43			0.98			0.35			0.4			0.99			
## 7000			-138.8			-5.76			31.32			0.47			2			0.44			0.98			0.36			0.39			0.99			
## 8000			-139.9			-5.76			33.6			0.46			2			0.43			0.99			0.36			0.39			0.99			
## 9000			-143.17			-5.89			44.91			0.65			2			0.44			0.98			0.35			0.4			0.99			
## 10000			-143.85			-5.81			19.03			0.36			2			0.44			0.98			0.36			0.4			0.99			
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
## lnL       -140.71223 1.829723 0.0691077      0.2094088          76.35
## prior       -5.79445 0.073642 0.0027814      0.0096498          58.24
## alpha        0.01980 0.004796 0.0001811      0.0004733         102.70
## sig2         0.01842 0.002879 0.0001087      0.0002851         102.01
## k            2.00000 0.000000 0.0000000      0.0000000           0.00
## ntheta       3.00000 0.000000 0.0000000      0.0000000           0.00
## root         3.70618 0.155625 0.0058779      0.0175721          78.44
## all theta    3.61348 0.356982        NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##     pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 409  1               3.643            0.01124         1.37
## 411  1               3.492            0.01872         2.16
```

```r
phenogram.density(tree, dat, chain = chain.fixed, burnin = 0.3, pp.cutoff = 0.5)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-263.png) 

```r
plotSimmap.mcmc(tree, chain.fixed, burnin = 0.3, circle = TRUE, fsize = 0.5)
```

```
## Error: object 'L' not found
```





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
## [1] -163.8
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
## [1] 35.02
```


## Converting to and from *OUwie* format

*bayou* includes some utilities to switch between *OUwie* and *bayou* formatted models. 

```r
require(OUwie)
```

```
## Loading required package: OUwie
## Loading required package: nloptr
```

```r
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
## lnL       -140.71223 1.829723 0.0691077      0.2094088          76.35
## prior       -5.79445 0.073642 0.0027814      0.0096498          58.24
## alpha        0.01980 0.004796 0.0001811      0.0004733         102.70
## sig2         0.01842 0.002879 0.0001087      0.0002851         102.01
## k            2.00000 0.000000 0.0000000      0.0000000           0.00
## ntheta       3.00000 0.000000 0.0000000      0.0000000           0.00
## root         3.70618 0.155625 0.0058779      0.0175721          78.44
## all theta    3.61348 0.356982        NA             NA             NA
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
## [1] 0.07303
## 
## $sig2
## [1] 0.03083
## 
## $k
## [1] 15
## 
## $ntheta
## [1] 16
## 
## $theta
##  [1] 3.719 3.210 4.528 4.282 3.021 5.774 3.826 1.655 4.645 4.196 5.217
## [12] 4.615 3.753 3.348 4.326 3.025
## 
## $sb
##  [1] 406 408  66 441  42  80 367 365 149 232 336 429  22 346 284
## 
## $loc
##  [1]  4.10190 20.33575 65.33020  4.18944 22.54083  0.04220  3.51999
##  [8]  3.36531  5.41938  9.09958  0.05487 81.71603 17.09800  1.12861
## [15]  4.66872
## 
## $t2
##  [1]  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
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
```

```
## Loading required package: denstrip
```

```r
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


