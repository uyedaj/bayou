# Tutorial for the R package *bayou*
The purpose of *bayou* is to fit Bayesian models of adaptive evolution to phylogenetic comparative data. Specifically, *bayou* provides a flexible framework for fitting multi-optima Ornstein-Uhlenbeck models to phylogenetic comparative data. This tutorial demonstrates some of the options for running *bayou*.

## Reversible-jump MCMC over regime placement
In this example, we will fit a reversible-jump MCMC model to an included dataset (the Chelonia dataset of Jaffe et al. 2012). To start, we will specify a model in which no parameters are fixed. This will estimate the posterior of shift number, location and magnitude as well as all other parameters.

We begin by loading the package and data. We will also assume a constant standard error across all taxa. This can instead be a named vector with species-specific measurement error.

```r
require(devtools)
install_github("bayou", username="uyedaj")
```

```
## Installing github repo bayou/master from uyedaj
## Downloading master.zip from https://github.com/uyedaj/bayou/archive/master.zip
## Installing package from /tmp/RtmpFNOXjw/master.zip
## arguments 'minimized' and 'invisible' are for Windows only
## Installing bayou
## '/usr/lib/R/bin/R' --vanilla CMD INSTALL  \
##   '/tmp/RtmpFNOXjw/devtoolsfd7446881b7/bayou-master'  \
##   --library='/home/josef/R/x86_64-pc-linux-gnu-library/3.1'  \
##   --install-tests 
## 
## Reloading installed bayou
```

Now load the package and an example dataset


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
prior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=1,prob=1), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 
The figure produced gives a rough visual of the chosen prior distributions. 

### Running the MCMC
Now we are going to run the mcmc. We are going to output files to our working directory. If you want to specify another director, replace "getwd()" with the path of the directory. By default, *bayou outputs to the R temporary directory. We will run a relatively short chain of only 10,000 generations.

```r
par(mfrow=c(2,3))
fit1 <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior, ngen=10000, new.dir=TRUE, plot.freq=2000, ticker.freq=1000)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1000			-143.33			-78.41			48.52			0.58			11			0.5			0.03			0.88			0.59			0.93			0.28			0.78			1			1			0.67			
```

```
## 2000			-137.23			-109.13			31.26			0.37			17			0.46			0.03			0.81			0.49			0.96			0.25			0.83			0.96			0.94			0.84			
## 3000			-125.37			-72.55			32.22			0.39			10			0.43			0.03			0.83			0.47			0.98			0.3			0.82			0.96			0.75			0.81			
```

```
## 4000			-120.16			-69.6			19.39			0.31			10			0.43			0.03			0.85			0.44			0.96			0.32			0.8			0.96			0.65			0.64			
## 5000			-122.85			-68.65			21.42			0.28			10			0.42			0.03			0.86			0.44			0.97			0.31			0.78			0.97			0.63			0.65			
```

```
## 6000			-117.82			-48.37			14.15			0.27			6			0.42			0.02			0.86			0.42			0.95			0.31			0.74			0.97			0.61			0.47			
## 7000			-109			-68.65			8.5			0.19			10			0.4			0.02			0.87			0.39			0.94			0.3			0.72			0.98			0.54			0.41			
```

```
## 8000			-105.35			-96.35			12.28			0.2			15			0.4			0.02			0.88			0.39			0.88			0.29			0.69			0.97			0.54			0.38			
## 9000			-105.97			-84.47			8.42			0.23			13			0.39			0.02			0.86			0.38			0.83			0.3			0.68			0.98			0.53			0.36			
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```
## 10000			-81.2			-110.36			4.01			0.16			18			0.38			0.02			0.86			0.38			0.81			0.3			0.66			0.97			0.5			0.35			
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
## /tmp/RtmpFNOXjw/NZAPSOURJF/bayou.* 
## To load results, use 'load.bayou(bayouFit)'
## 
## 10000  generations were run with the following acceptance probabilities:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     0.38     0.02     0.86     0.38     0.81     0.30     0.66     0.97 
## U1.slide U2.slide 
##     0.50     0.35 
##  Total number of proposals of each type:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     1842     4405      194      292      123      877     1824      186 
## U1.slide U2.slide 
##      105      152
```

We can load the actual chains by running the following code:

```r
chain <- load.bayou(fit1, save.Rdata=FALSE, cleanup=TRUE)
```

```
## deleting temporary directory /tmp/RtmpFNOXjw/NZAPSOURJF/
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
## lnL       -113.48446 11.61540 0.4387076       5.522865          4.423
## prior      -75.10914 16.21010 0.6122470       7.943999          4.164
## alpha        0.06994  0.04543 0.0017158       0.020530          4.896
## sig2         0.03106  0.01288 0.0004864       0.005718          5.072
## k           11.14836  3.14201 0.1186722       1.535965          4.185
## ntheta      12.14836  3.14201 0.1186722       1.535965          4.185
## root         3.28963  0.14086 0.0053203       0.065534          4.620
## all theta    3.82311  1.04469        NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##         pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 408 1.0000               5.001           0.010483     1.211808
## 372 0.9743               4.369           0.011651     1.746000
## 412 0.7275               3.780           0.007811     0.379256
## 269 0.5720               4.727           0.052037     2.721977
## 116 0.4465               3.042           0.047721    17.049139
## 45  0.4151               4.615           0.009777     0.370536
## 232 0.2582               4.956           0.042973     3.487464
## 418 0.2582               3.361           0.034940     1.299152
## 272 0.2411               3.486           0.060388     1.342861
## 270 0.2254               3.849           0.048071     0.496541
## 2   0.2211               4.643           0.061307     0.147731
## 329 0.2183               4.268           0.086966     0.122718
## 252 0.2140               4.169           0.018954     1.117629
## 152 0.1997               2.897           0.056153     0.016277
## 43  0.1983               4.806           0.024630     0.194038
## 186 0.1983               1.179           0.103221     0.397393
## 438 0.1969               3.062           0.087867     3.884027
## 57  0.1940               3.827           0.032137     0.892684
## 435 0.1854               3.915           0.073428     0.766506
## 304 0.1655               3.078           0.058166     0.346133
## 391 0.1641               3.487           0.063434     1.851090
## 436 0.1469               3.340           0.052854     0.379028
## 426 0.1312               3.846           0.021095     2.295718
## 424 0.1213               4.313           0.022116     3.606174
## 417 0.1198               4.152           0.054690     0.529499
## 302 0.1141               4.455           0.098862     0.041822
## 188 0.1127               2.794           0.089926     0.015124
## 306 0.1127               2.275           0.028443     1.008554
## 76  0.1098               2.422           0.074394     0.004767
## 44  0.1013               4.666           0.041798     1.189356
```

We can visualize the chains by viewing traces for each parameter using the plotting utilities of the R package coda:

```r
plot(chain)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 

We can view where there are shifts of high probability.

```r
par(mfrow=c(1,1))
plotSimmap.mcmc(tree, chain, burnin=0.3, circle=TRUE,fsize=0.4)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

And we can view the density of phenotypic optima and location of highly supported shifts on the phenogram. Here we show all shifts with posterior probabilities greater than *pp.cutoff = 0.3*. 

```r
phenogram.density(tree, dat, chain=chain, burnin=0.3, pp.cutoff=0.3)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

### Diagnosing Convergence
It is useful to run two chains with independent starting positions. By default, *bayou* simulates starting parameters from the prior distribution. We will re-run the analysis here to obtain an independent chain.

```r
fit2 <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior, ngen=10000, new.dir=TRUE, plot.freq=NULL, ticker.freq=1000)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1000			-144.18			-73.76			9.67			0.34			11			0.37			0.05			0.82			0.67			0.95			0.39			0.66			1			0.46			0.64			
## 2000			-104.06			-94.72			6.85			0.17			15			0.34			0.04			0.85			0.53			0.84			0.35			0.64			0.96			0.46			0.45			
## 3000			-110.61			-69.66			9.7			0.22			10			0.33			0.03			0.9			0.44			0.81			0.34			0.62			0.96			0.39			0.37			
## 4000			-109.39			-74.06			11.79			0.22			11			0.34			0.03			0.92			0.43			0.76			0.32			0.62			0.97			0.4			0.35			
## 5000			-104.17			-69.17			10.46			0.2			10			0.33			0.03			0.93			0.4			0.75			0.32			0.6			0.96			0.38			0.29			
## 6000			-95.34			-89.43			4.38			0.16			14			0.33			0.03			0.93			0.36			0.69			0.33			0.59			0.9			0.37			0.29			
## 7000			-86.22			-94.72			3.55			0.14			15			0.33			0.02			0.9			0.34			0.7			0.34			0.55			0.88			0.35			0.26			
## 8000			-86.99			-131.3			6.02			0.16			22			0.33			0.03			0.88			0.33			0.67			0.34			0.54			0.88			0.33			0.27			
## 9000			-96.18			-104.83			3.8			0.19			17			0.33			0.03			0.85			0.32			0.69			0.33			0.54			0.88			0.33			0.28			
## 10000			-94.82			-90.42			7.5			0.21			14			0.34			0.03			0.85			0.3			0.69			0.33			0.53			0.88			0.3			0.27			
```

```r
chain2 <- load.bayou(fit2, save.Rdata=FALSE, cleanup=FALSE)
chain2 <- set.burnin(chain2, 0.3)
```

Now we can compare chains to see if the parameters have converged using Gelman and Rubin's R statistic. Values close to 1 indicate chains have converged (these two chains will not have converged in only 10,000 generations).

```r
RlnL <- gelman.R("lnL", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```r
Ralpha <- gelman.R("alpha", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 

```r
Rsig2 <- gelman.R("sig2", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) 

Of particular interest to inference are the branch posterior probabilities, we can plot these to see if the two chains converged on similar answers. If the runs have converged, points should fall along the y=x line.

```r
L1 <- Lposterior(chain,tree, burnin=0.3)
L2 <- Lposterior(chain2,tree, burnin=0.3)
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 

We can also combine chains:

```r
chains <- combine.chains(chain, chain2, burnin.prop=0.3)
chains <- set.burnin(chains, 0)
```

### Estimating marginal likelihoods
To compare models, we want to estimate the marginal likelihood, which is performed using the stepping stone algorithm in *bayou*. The method first estimates a reference function by fitting a series of curves to the posterior of an MCMC chain. The function *steppingstone* will output a graphic showing the best-fitting density functions overlaying the posterior distribution. Then the stepping stone MCMC's are run for every value in the vector *Bk*, corresponding to each value of the power posterior function ranging from the reference function (*Bk = 0*) to the posterior (*Bk = 1*). 
To speed things up, you can specify multiple cores by registering a parallel backend for the \code{foreach} package.

```r
#require(doParallel)
#registerDoParallel(cores=5)
ss <- steppingstone(Bk=seq(0,1,length.out=5), chains, tree, dat, SE=SE, prior=prior, new.dir=TRUE, ngen=10000, parallel=FALSE) #set parallel to TRUE if desired.
```

```
## Making power posterior function from provided mcmc chain...
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

```
## Running mcmc chains...
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0			1000			-217.24			-53.22			-46.38			20.46			0.7			7			0.78			0			1			1			1			0.62			0.72			1			1			1			
## 0			2000			-2189.9			-48.8			-40.25			6.19			0.1			6			0.77			0			1			1			1			0.67			0.76			1			1			1			
## 0			3000			-190.25			-58.8			-51.05			9.53			0.32			8			0.78			0.01			1			1			1			0.67			0.76			1			1			1			
## 0			4000			-523.46			-48.17			-38.2			6.8			0.15			6			0.79			0.01			1			1			1			0.67			0.78			1			1			1			
## 0			5000			-1976.51			-33.06			-27.31			6.54			0.09			3			0.8			0.01			1			1			1			0.69			0.78			1			1			1			
## 0			6000			-1065.3			-27.77			-22.16			5.42			0.06			2			0.79			0.01			1			1			1			0.7			0.77			1			1			1			
## 0			7000			-230.45			-27.58			-19.6			10.51			0.27			2			0.79			0.01			1			1			1			0.69			0.78			1			1			1			
## 0			8000			-191.84			-47.95			-36.96			9.36			0.19			6			0.79			0.01			1			1			1			0.7			0.78			1			1			1			
## 0			9000			-1357.35			-53.84			-46.48			6.57			0.27			7			0.79			0.01			1			1			1			0.71			0.77			1			1			1			
## 0			10000			-604.92			-33.72			-26.69			8.12			0.2			3			0.79			0.01			1			1			1			0.7			0.77			1			1			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.25			1000			-148.74			-55.97			-58.35			62.73			0.72			7			0.61			0.01			1			0.87			1			0.68			0.74			1			0.88			1			
## 0.25			2000			-150.25			-23.92			-20.16			18.89			0.36			1			0.62			0			0.97			0.82			1			0.59			0.69			1			0.93			1			
## 0.25			3000			-144.68			-32.65			-25.9			18.01			0.38			3			0.62			0.01			0.98			0.79			1			0.55			0.7			1			0.95			1			
## 0.25			4000			-149.72			-28.04			-22.75			20.17			0.47			2			0.6			0.01			0.99			0.81			1			0.57			0.69			1			0.96			1			
## 0.25			5000			-149.24			-43.28			-36.4			23.61			0.58			5			0.61			0.01			0.99			0.81			1			0.56			0.69			1			0.96			1			
## 0.25			6000			-150.19			-49.43			-47.18			30.88			0.5			6			0.6			0.01			0.99			0.82			1			0.56			0.69			1			0.97			1			
## 0.25			7000			-140.1			-48.3			-44.99			24.81			0.39			6			0.6			0.01			0.99			0.82			1			0.55			0.69			1			0.98			1			
## 0.25			8000			-138.58			-48.19			-40.49			18.43			0.4			6			0.6			0.01			0.99			0.82			1			0.55			0.68			1			0.98			0.97			
## 0.25			9000			-146.86			-28			-23.6			23.72			0.33			2			0.6			0.01			0.99			0.82			1			0.56			0.68			1			0.98			0.92			
## 0.25			10000			-145.6			-43.46			-36.87			26.23			0.45			5			0.6			0.01			1			0.82			1			0.55			0.68			1			0.98			0.92			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.5			1000			-146.22			-69.3			-71.71			87.88			0.92			10			0.59			0.02			1			0.73			0.9			0.47			0.77			1			0.62			0.8			
## 0.5			2000			-138.5			-33.97			-31.07			25.29			0.4			3			0.53			0.01			0.94			0.78			0.95			0.46			0.72			0.98			0.79			0.94			
## 0.5			3000			-140.98			-47.79			-40.44			31.34			0.52			6			0.55			0.01			0.95			0.7			0.96			0.45			0.75			0.97			0.79			0.97			
## 0.5			4000			-139.19			-43.6			-39.83			33.94			0.44			5			0.54			0.01			0.96			0.69			0.97			0.44			0.74			0.98			0.82			0.98			
## 0.5			5000			-140.8			-59.66			-56.63			31.75			0.56			8			0.53			0.01			0.95			0.7			0.97			0.45			0.74			0.98			0.87			0.98			
## 0.5			6000			-149.31			-38			-31.77			19.83			0.53			4			0.53			0.01			0.95			0.69			0.98			0.44			0.74			0.99			0.89			0.95			
## 0.5			7000			-140.58			-48.18			-42.09			36.37			0.53			6			0.53			0.01			0.95			0.73			0.98			0.44			0.74			0.98			0.9			0.96			
## 0.5			8000			-139.7			-48.17			-41.89			27.51			0.46			6			0.53			0.01			0.96			0.73			0.98			0.44			0.74			0.97			0.91			0.96			
## 0.5			9000			-136.02			-59.16			-51.38			13.6			0.25			8			0.52			0.01			0.96			0.72			0.97			0.45			0.73			0.97			0.9			0.91			
## 0.5			10000			-125.71			-38.15			-31.31			16.24			0.25			4			0.52			0.01			0.96			0.69			0.96			0.45			0.72			0.97			0.87			0.81			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.75			1000			-141.72			-43.49			-47.12			56.24			0.56			5			0.53			0.01			0.91			0.7			1			0.47			0.82			1			0.73			1			
## 0.75			2000			-137.66			-48.93			-42.44			22.85			0.42			6			0.52			0.01			0.94			0.78			1			0.39			0.76			1			0.55			0.65			
## 0.75			3000			-133.31			-32.66			-23.81			9.63			0.27			3			0.47			0.01			0.91			0.69			0.91			0.38			0.7			1			0.47			0.62			
## 0.75			4000			-132.2			-63.1			-50.74			12.38			0.35			9			0.47			0.01			0.91			0.6			0.89			0.38			0.68			1			0.45			0.52			
## 0.75			5000			-122.61			-68.18			-56.83			12.38			0.21			10			0.47			0.02			0.9			0.57			0.87			0.38			0.66			0.99			0.44			0.55			
## 0.75			6000			-123.53			-53.71			-47.83			25.12			0.35			7			0.46			0.02			0.91			0.53			0.88			0.37			0.64			0.94			0.42			0.51			
## 0.75			7000			-124.55			-68.72			-58.32			15.07			0.3			10			0.44			0.02			0.91			0.51			0.86			0.39			0.62			0.94			0.42			0.46			
## 0.75			8000			-118			-79.02			-67.62			11.8			0.29			12			0.44			0.02			0.91			0.5			0.85			0.4			0.62			0.94			0.42			0.48			
## 0.75			9000			-116.29			-73.71			-64.83			12.71			0.24			11			0.44			0.02			0.91			0.5			0.84			0.39			0.61			0.94			0.43			0.45			
## 0.75			10000			-111.52			-103.53			-89.18			7.98			0.2			17			0.43			0.02			0.89			0.48			0.83			0.38			0.62			0.95			0.41			0.46			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1			1000			-134.75			-116.05			-115.6			38.86			0.47			19			0.48			0.04			0.89			0.65			0.92			0.47			0.72			1			0.81			0.83			
## 1			2000			-136.51			-88.68			-83.03			26.27			0.37			14			0.41			0.03			0.94			0.52			0.88			0.41			0.74			1			0.76			0.5			
## 1			3000			-139.57			-86.01			-85.42			32.38			0.5			13			0.39			0.03			0.96			0.55			0.91			0.38			0.77			1			0.86			0.67			
## 1			4000			-126.16			-96.13			-94.12			30.71			0.35			15			0.39			0.03			0.93			0.59			0.9			0.37			0.79			1			0.88			0.67			
## 1			5000			-123.03			-85.87			-81.28			23.16			0.32			13			0.38			0.03			0.9			0.55			0.88			0.37			0.78			0.99			0.83			0.6			
## 1			6000			-121.26			-80.75			-81.14			22.89			0.31			12			0.36			0.03			0.92			0.54			0.88			0.37			0.78			0.99			0.81			0.54			
## 1			7000			-128.91			-54.62			-52.06			22.47			0.3			7			0.36			0.03			0.92			0.53			0.88			0.36			0.77			0.99			0.81			0.57			
## 1			8000			-105.32			-91.59			-84.55			11.14			0.17			14			0.37			0.03			0.92			0.51			0.88			0.36			0.75			0.99			0.8			0.51			
## 1			9000			-102.08			-75.41			-67.84			10.31			0.2			11			0.37			0.03			0.91			0.49			0.87			0.35			0.72			0.99			0.74			0.46			
## 1			10000			-94.23			-83.93			-72.57			10.2			0.2			13			0.36			0.03			0.89			0.48			0.86			0.35			0.7			0.99			0.71			0.42			
## Loading mcmc chains...
```

```r
ss
```

```
## Stepping stone estimation of marginal likelihood
## Marginal Likelihood:
## [1] -142.8
## A total of 5 power posteriors were run along the sequence: 0		0.25		0.5		0.75		1
## lnr_k -42.52 -36 -33.61 -30.63
```

```r
ss <- set.burnin(ss, 0.3)
```

Plotting the stepping stone object will give you an idea of convergence for each of the MCMC's by providing a trace of the lnL, ln prior, and reference function as well as the estimated marginal likelihood at each step.

```r
plot(ss)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 

## Fitting models with fixed parameters
*bayou* also allows the user to fit models with fixed parameters, including non-reversible jump models. A helpful utility for setting up fixed hypotheses is the function *identifyBranches*. First, we set up some starting parameters.

```r
startpar <- list(alpha=1, sig2=1, k=2, ntheta=3, theta=c(4,5,6))
```

Then we use *identifyBranches()* to specify the location of shifts by clicking on them. f

```r
fixed.hypothesis <- identifyBranches(tree, startpar$k)
```




```r
startpar$sb <- fixed.hypothesis$sb
startpar$t2 <- 2:startpar$ntheta
startpar$loc <- fixed.hypothesis$loc
plotBayoupars(startpar, tree, fsize=0.5)
```

```
## no colors provided. using the following legend:
##        1        2        3 
##  "black"    "red" "green3"
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 

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
converge_pars$loc <- c(converge_pars$loc, 0)
converge_pars$t2 <- c(converge_pars$t2, 3)
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
## [1] 2 3 3
## 
## $loc
## [1] 24  8  0
```

```r
plotBayoupars(converge_pars, tree, fsize=0.5)
```

```
## no colors provided. using the following legend:
##        1        2        3 
##  "black"    "red" "green3"
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 

Now that we know how to specify fixed models, lets set up the prior function:

```r
prior.fixed <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="fixed", dk="fixed", dtheta="dnorm", dloc="dloc"), 
                            param=list(dalpha=list(scale=1), dsig2=list(scale=1), dtheta=list(mean=mean(dat), sd=2)),
                                     fixed=fixed.hypothesis)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 

A few things to note. First, we have set *dsb* and *dk* to "fixed". We could have also set *dloc*, *dalpha* and/or *dsig2* to "fixed". Here, we will let the location of the shift on the branch vary (*dloc* puts a uniform prior on the location of the shift on the branch). Second, we are passing our model through the *fixed.hypothesis* object we created using *identifyBranches()*. For example, since simple Brownian motion is a special case of multi-optima OU models, we can specify through the the following prior:

```r
prior.BM <- make.prior(tree, dists=list(dalpha="fixed", dsig2="dhalfcauchy", dsb="fixed", dk="fixed", dloc="fixed"),
                            param=list(dsig2=list(scale=1), dtheta=list(mean=mean(dat), sd=2)),
                                fixed=list(alpha=0, k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0)))
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 

Now we can run our MCMC chains using our priors. Right now, you have to specify the starting parameters. This will soon be unnecessary...

```r
fit.fixed <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior=prior.fixed, startpar=startpar, ngen=10000, new.dir=TRUE, plot.freq=NULL)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 1000			-141.73			-5.82			24.94			0.39			2			0.4			0.99			0.42			0.35			1			
## 2000			-140.01			-5.75			37.85			0.54			2			0.4			0.95			0.44			0.39			1			
## 3000			-143.27			-5.85			29.08			0.51			2			0.41			0.96			0.4			0.41			0.99			
## 4000			-138.89			-5.79			32.79			0.42			2			0.41			0.97			0.39			0.4			0.99			
## 5000			-140.99			-5.74			34.61			0.45			2			0.41			0.97			0.38			0.42			1			
## 6000			-139.31			-5.76			39.79			0.5			2			0.41			0.98			0.37			0.41			1			
## 7000			-140.37			-5.83			46.7			0.58			2			0.42			0.97			0.37			0.41			0.99			
## 8000			-141.33			-5.75			41.26			0.6			2			0.43			0.97			0.36			0.41			0.99			
## 9000			-141.05			-5.75			47.34			0.44			2			0.43			0.97			0.37			0.42			0.99			
## 10000			-138.49			-5.75			32.1			0.44			2			0.44			0.97			0.36			0.42			0.99			
```




```r
chain.fixed <- load.bayou(fit.fixed, save.Rdata=FALSE, cleanup=FALSE)
chain.fixed <- set.burnin(chain.fixed, 0.3)
plot(chain.fixed)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-271.png) ![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-272.png) 

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
## lnL       -140.87327 1.731008 0.0653793      0.1895640          83.38
## prior       -5.79995 0.058820 0.0022216      0.0070591          69.43
## alpha        0.01889 0.004813 0.0001818      0.0006656          52.29
## sig2         0.01788 0.002863 0.0001081      0.0003333          73.80
## k            2.00000 0.000000 0.0000000      0.0000000           0.00
## ntheta       3.00000 0.000000 0.0000000      0.0000000           0.00
## root         3.69640 0.180655 0.0068232      0.0156335         133.53
## all theta    3.60502 0.380033        NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##     pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 409  1               3.681            0.01282        1.383
## 411  1               3.438            0.01865        2.208
```

```r
phenogram.density(tree, dat, chain=chain.fixed, burnin=0.3, pp.cutoff=0.5)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-273.png) 

```r
plotSimmap.mcmc(tree, chain.fixed, burnin=0.3, circle=TRUE, fsize=0.5)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-274.png) 



As before, we can estimate the marginal likelihood:

```r
ss.fixed <- steppingstone(Bk=seq(0,1,length.out=5), chain.fixed, tree, dat, SE=0, startpar=startpar, prior=prior.fixed, ngen=10000, parallel=FALSE)
```

```
## Making power posterior function from provided mcmc chain...
```

```
## Running mcmc chains...
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0			1000			-157.28			-5.78			8.99			31.7			0.39			2			0.59			1			0.52			0.52			1			
## 0			2000			-156.89			-5.77			7.79			30.31			0.44			2			0.58			1			0.51			0.5			1			
## 0			3000			-162.12			-5.8			7.94			43.19			0.49			2			0.61			1			0.48			0.5			1			
## 0			4000			-160.21			-5.78			6.18			44.81			0.69			2			0.62			1			0.49			0.51			1			
## 0			5000			-174.56			-6.11			1.78			38.04			0.42			2			0.63			1			0.49			0.51			1			
## 0			6000			-158.69			-5.77			7.67			50.39			0.66			2			0.64			1			0.48			0.52			1			
## 0			7000			-164.59			-6.03			3.61			36.56			0.47			2			0.64			1			0.47			0.51			1			
## 0			8000			-181.22			-5.98			2.59			33.82			0.38			2			0.64			1			0.46			0.51			1			
## 0			9000			-153.18			-5.79			8.39			36.01			0.53			2			0.65			1			0.45			0.51			1			
## 0			10000			-154.54			-5.78			9.1			38.04			0.51			2			0.65			1			0.46			0.51			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0.25			1000			-153.91			-5.8			8.39			35.47			0.5			2			0.57			1			0.49			0.43			0.99			
## 0.25			2000			-153.48			-5.79			8.25			34.17			0.5			2			0.59			1			0.48			0.43			0.99			
## 0.25			3000			-159.02			-5.95			5.2			44.45			0.62			2			0.58			1			0.44			0.44			1			
## 0.25			4000			-159.88			-5.76			7.47			52.13			0.71			2			0.58			1			0.42			0.43			0.99			
## 0.25			5000			-151.56			-5.76			7.77			36.01			0.59			2			0.58			1			0.42			0.43			0.99			
## 0.25			6000			-152.59			-5.76			9.17			37.43			0.54			2			0.58			1			0.42			0.44			0.99			
## 0.25			7000			-151.03			-5.75			8.01			26.86			0.4			2			0.58			1			0.41			0.44			0.99			
## 0.25			8000			-156.52			-5.79			6.28			40.08			0.66			2			0.59			1			0.41			0.44			0.99			
## 0.25			9000			-158.74			-5.81			6.11			52.2			0.74			2			0.59			1			0.4			0.44			0.99			
## 0.25			10000			-160.91			-5.9			6.18			31.31			0.39			2			0.58			1			0.4			0.44			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0.5			1000			-152.71			-5.78			5.78			26.67			0.44			2			0.56			0.98			0.42			0.45			1			
## 0.5			2000			-150.62			-5.75			8.18			28.56			0.44			2			0.54			0.98			0.37			0.45			0.99			
## 0.5			3000			-153.79			-5.84			4.8			33.51			0.54			2			0.53			0.99			0.4			0.44			1			
## 0.5			4000			-158.77			-5.79			5.14			50.48			0.9			2			0.53			0.99			0.39			0.44			0.99			
## 0.5			5000			-152.4			-5.76			8.3			30.77			0.43			2			0.53			0.99			0.39			0.42			1			
## 0.5			6000			-152.41			-5.79			7.26			33.12			0.53			2			0.53			0.99			0.38			0.42			1			
## 0.5			7000			-151.93			-5.75			6.36			33.87			0.6			2			0.53			0.99			0.39			0.41			0.99			
## 0.5			8000			-150.05			-5.78			5.56			26.95			0.44			2			0.53			0.99			0.38			0.41			1			
## 0.5			9000			-160.22			-6.04			3.06			44.08			0.66			2			0.53			0.99			0.37			0.42			1			
## 0.5			10000			-153			-5.82			6.69			32.86			0.53			2			0.53			0.99			0.37			0.42			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0.75			1000			-151.42			-5.77			8.58			32.36			0.48			2			0.44			1			0.31			0.32			1			
## 0.75			2000			-151.97			-5.79			-1.8			18.55			0.36			2			0.43			0.99			0.31			0.35			1			
## 0.75			3000			-150.25			-5.81			3.49			23.56			0.42			2			0.44			0.99			0.3			0.36			1			
## 0.75			4000			-156.21			-5.8			2.41			30.08			0.57			2			0.43			0.99			0.31			0.36			1			
## 0.75			5000			-153.26			-5.77			-0.51			23.81			0.55			2			0.45			0.99			0.32			0.37			1			
## 0.75			6000			-151.43			-5.78			4.58			24.26			0.44			2			0.45			0.99			0.32			0.38			1			
## 0.75			7000			-151.64			-5.8			5.02			27.99			0.51			2			0.44			0.98			0.32			0.37			1			
## 0.75			8000			-153.99			-5.81			5.81			38.59			0.69			2			0.44			0.98			0.33			0.37			1			
## 0.75			9000			-151.64			-5.8			3.13			24.17			0.44			2			0.44			0.98			0.33			0.37			1			
## 0.75			10000			-151.19			-5.76			7.13			26.32			0.42			2			0.44			0.98			0.33			0.37			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 1			1000			-150.75			-5.8			3.66			27.9			0.51			2			0.37			1			0.32			0.27			1			
## 1			2000			-150.69			-5.76			-7.28			15.98			0.38			2			0.35			0.99			0.31			0.31			0.99			
## 1			3000			-151.49			-5.79			-7.77			16.76			0.4			2			0.35			0.98			0.31			0.31			1			
## 1			4000			-148.37			-5.75			2.02			22.2			0.46			2			0.38			0.99			0.32			0.31			0.99			
## 1			5000			-152.07			-5.75			-16.58			13.18			0.38			2			0.37			0.99			0.31			0.32			1			
## 1			6000			-148.49			-5.75			0.14			20.3			0.43			2			0.38			0.99			0.32			0.32			0.99			
## 1			7000			-150.65			-5.77			5.72			28.61			0.48			2			0.39			0.99			0.31			0.33			0.99			
## 1			8000			-150.09			-5.75			4.22			21.38			0.37			2			0.39			0.99			0.3			0.34			0.99			
## 1			9000			-149.33			-5.76			-4.75			17.97			0.42			2			0.39			0.99			0.31			0.34			0.99			
## 1			10000			-151.66			-5.8			-1.86			23.55			0.56			2			0.39			0.98			0.31			0.34			0.99			
## Loading mcmc chains...
```

```r
ss.fixed <- set.burnin(ss.fixed, 0.3)
ss.fixed$lnr
```

```
## [1] -163.4
```

```r
plot(ss.fixed)
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-291.png) ![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-292.png) 

And calculate Bayes Factors in support of the Reversible-Jump model vs. the fixed model. 

```r
2*(ss$lnr-ss.fixed$lnr)
```

```
## [1] 41.3
```

## Converting to and from *OUwie* format

*bayou* includes some utilities to switch between *OUwie* and *bayou* formatted models. 

```r
require(OUwie)
OUwieData <- bayou2OUwie(startpar, tree, dat)
OUwieRes <- OUwie(OUwieData$tree, OUwieData$dat, model="OUM")
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
## lnL       -140.87327 1.731008 0.0653793      0.1895640          83.38
## prior       -5.79995 0.058820 0.0022216      0.0070591          69.43
## alpha        0.01889 0.004813 0.0001818      0.0006656          52.29
## sig2         0.01788 0.002863 0.0001081      0.0003333          73.80
## k            2.00000 0.000000 0.0000000      0.0000000           0.00
## ntheta       3.00000 0.000000 0.0000000      0.0000000           0.00
## root         3.69640 0.180655 0.0068232      0.0156335         133.53
## all theta    3.60502 0.380033        NA             NA             NA
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

![plot of chunk unnamed-chunk-33](figure/unnamed-chunk-33.png) 

## Some other useful things
### Some plotting features
It's easy to convert *bayou* formatted models into trees that can be visualized using *phytools*' *plotSimmap* or *phenogram* functions. For example, let's pull a few samples from the posterior out of our original chain:

```r
samp_pars <- pull.pars(500, chain, model="OU")
samp_pars
```

```
## $alpha
## [1] 0.03236
## 
## $sig2
## [1] 0.01831
## 
## $k
## [1] 10
## 
## $ntheta
## [1] 11
## 
## $theta
##  [1] 3.417 4.320 3.145 3.295 4.926 4.760 3.616 4.693 4.322 4.738 2.819
## 
## $sb
##  [1] 372 436 270 408 390 412   8 188 237 304
## 
## $loc
##  [1]  1.3559  6.2727  4.9575 15.0320  0.4805  4.1893 25.4025  0.8158
##  [9]  3.2784  3.0350
## 
## $t2
##  [1]  2  3  4  5  6  7  8  9 10 11
```

Now let's convert it to *phytools*' simmap format:

```r
sm_tree <- pars2simmap(samp_pars,tree)
plotSimmap(sm_tree$tree, col=sm_tree$col)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-351.png) 

```r
phenogram(sm_tree$tree, dat, col=sm_tree$col, ftype="off")
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-352.png) 

We can also make some regime plots, which show the location of each optimum and the expected stationary variance.

```r
plot(c(0, 210), c(2, 6), type="n", xlab="time", ylab="phenotype")
regime.plot(samp_pars, sm_tree$tree, type="density", cols=sm_tree$col)
```

```
## Error: could not find function "regime.plot"
```

```r
phenogram(sm_tree$tree, dat, col=sm_tree$col, ftype="off", add=TRUE)
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-36.png) 

....Put in how to do OU ancestral state reconstructions here....


## How to specify some other models
Rather than restricting our interest to only models with a single shift allowed per branch, we can fit a number of models that allow multiple shifts per branch, have unequal probabilities, etc. These are all specified with the prior function. Here are some examples:

A model with an arbitrary number of shifts allowed per branch, and branches chosen with probability proportional to their length:

```r
tree <- reorder(tree, "postorder")
dat <- dat[tree$tip.label]
prior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=Inf,prob=tree$edge.length), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-37.png) 

A model that disallows shifts on terminal branches:

```r
terminal.branches <- as.numeric(tree$edge[,2] < length(tree$tip.label)+1)
prior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=terminal.branches, prob=1), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-38.png) 


# Alternative parameterizations

```r
prior.repar <- make.prior(tree, model="OUrepar", dists=list(dhalflife="dlnorm", dVy="dlnorm",dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                          param=list(dhalflife=list(meanlog=1.5, sdlog=1.5), dVy=list(meanlog=-2, sdlog=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=1, prob=1), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-391.png) 

```r
fit.ourepar <- bayou.mcmc(tree, dat, SE=SE, model="OUrepar", prior=prior.repar, ngen=10000, new.dir=TRUE)
```

```
## gen			lnL			prior			half.life			Vy			K			birth.k			D0.slide			D1.slide			death.k			.halflife			.theta			U0.slide			U1.slide			U2.slide			.Vy			
## 1000			-157.38			-107.24			363.82			2.9			13			0.04			1			0.71			0.94			0.42			0.81			0.95			0.5			0.75			0.4			
```

```
## 2000			-142.96			-79.27			51.1			0.53			10			0.03			0.96			0.64			0.96			0.38			0.82			0.98			0.62			0.8			0.31			
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-392.png) 

```
## 3000			-123.82			-123.18			22.68			0.3			19			0.03			0.9			0.57			0.97			0.41			0.81			0.94			0.53			0.79			0.34			
```

```
## 4000			-102.1			-127.78			10.53			0.17			21			0.03			0.89			0.54			0.9			0.4			0.76			0.94			0.54			0.68			0.35			
```

```
## 5000			-112.01			-83.42			17.71			0.26			12			0.03			0.89			0.5			0.91			0.42			0.73			0.91			0.5			0.64			0.36			
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-393.png) 

```
## 6000			-109.28			-96.45			10.15			0.22			15			0.03			0.9			0.47			0.86			0.42			0.69			0.92			0.48			0.61			0.34			
```

```
## 7000			-112.43			-113.87			12.59			0.22			18			0.03			0.88			0.47			0.81			0.44			0.69			0.92			0.47			0.6			0.32			
```

```
## 8000			-99.39			-101.8			8.42			0.21			16			0.03			0.88			0.46			0.8			0.44			0.68			0.93			0.48			0.58			0.32			
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-394.png) 

```
## 9000			-97.46			-105.96			8.03			0.21			17			0.02			0.88			0.45			0.79			0.45			0.67			0.92			0.5			0.55			0.32			
```

```
## 10000			-86.82			-101			5.22			0.18			16			0.02			0.9			0.43			0.78			0.46			0.66			0.92			0.45			0.54			0.32			
```

```r
chain.ourepar <- load.bayou(fit.ourepar, save.Rdata=FALSE, cleanup=FALSE)
chain.ourepar <- set.burnin(chain.ourepar, 0.3)
plot(chain.ourepar)
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-395.png) ![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-396.png) ![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-397.png) 

```r
out.ourepar <- summary(chain.ourepar)
```

```
## bayou MCMC chain: 10000 generations
## 1000 samples, first 300 samples discarded as burnin
## 
## 
## Summary statistics for parameters:
##                Mean       SD Naive SE Time-series SE Effective Size
## lnL       -108.2104  8.71064 0.328996        3.61262          5.814
## prior     -103.2710 14.54810 0.549474        6.21921          5.472
## halflife    12.8217  4.52759 0.171005        1.32636         11.652
## Vy           0.2234  0.03874 0.001463        0.01203         10.374
## k           16.0642  2.72400 0.102884        1.15559          5.557
## ntheta      17.0642  2.72400 0.102884        1.15559          5.557
## root         3.6511  0.08047 0.003039        0.02926          7.566
## all theta    3.6192  1.13994       NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##         pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 372 1.0000               4.362           0.012182       0.3029
## 45  0.8787               4.737           0.008503       0.5251
## 407 0.5835               3.227           0.008806       0.4967
## 409 0.4893               4.479           0.030078       0.3394
## 58  0.4522               2.919           0.106104       0.7227
## 411 0.4451               4.246           0.046653       0.6646
## 53  0.3766               2.971           0.009821       0.4012
## 223 0.3595               2.409           0.039592       0.6958
## 252 0.3566               4.227           0.019732       0.4452
## 296 0.3338               2.747           0.029311       0.3459
## 441 0.3096               3.071           0.013560       0.7214
## 408 0.2967               4.831           0.017811       0.2426
## 405 0.2725               3.644           0.021396       0.4038
## 110 0.2682               5.248           0.057341       0.5575
## 404 0.2639               3.318           0.024308       0.3240
## 199 0.2482               2.586           0.098746       0.3057
## 424 0.2468               3.353           0.049496       0.6537
## 36  0.2268               4.890           0.073628       0.5971
## 181 0.2268               2.516           0.060687       0.2010
## 406 0.2225               3.162           0.003260       0.1511
## 234 0.2026               2.941           0.044137       0.3316
## 143 0.2011               1.951           0.040929       0.3864
## 226 0.2011               2.196           0.026369       0.4325
## 417 0.1969               6.012           0.030747       0.3110
## 136 0.1940               2.293           0.035075       0.5557
## 231 0.1883               2.933           0.188098       0.7416
## 293 0.1840               2.716           0.045446       0.5175
## 193 0.1769               3.102           0.048747       0.4558
## 28  0.1755               4.592           0.036175       0.5564
## 429 0.1755               3.340           0.019950       0.5747
## 266 0.1641               3.564           0.068064       0.5876
## 139 0.1612               2.716           0.023795       0.3924
## 54  0.1598               3.815           0.038943       0.2601
## 142 0.1598               3.836           0.053480       0.4774
## 325 0.1541               2.221           0.047904       0.5561
## 52  0.1526               3.827           0.062403       0.7291
## 42  0.1455               5.013           0.072816       0.3700
## 282 0.1455               3.466           0.097267       0.2854
## 333 0.1427               2.647           0.087986       0.3333
## 18  0.1412               4.943           0.027858       0.8437
## 225 0.1398               4.168           0.117400       0.8700
## 268 0.1255               3.330           0.051778       0.5499
## 427 0.1241               2.669           0.031077       0.3172
## 255 0.1198               2.981           0.019442       0.7165
## 306 0.1198               1.946           0.151294       0.2813
## 376 0.1184               4.792           0.097345       0.2698
## 336 0.1098               2.229           0.144106       0.7011
## 440 0.1070               2.537           0.016182       0.3523
## 334 0.1013               4.134           0.077875       0.8622
## 442 0.1013               4.657           0.047849       0.5450
```

```r
phenogram.density(tree, dat, chain=chain.ourepar, burnin=0.3, pp.cutoff=0.5)
plotSimmap.mcmc(tree, chain.ourepar, burnin=0.3, circle=TRUE, fsize=0.5)

prior.QG <- make.prior(tree, model="QG", dists=list(dh2="dbeta", dP="dlnorm", dw2="dlnorm", dNe="dlnorm", dk="cdpois", dtheta="dnorm"),
                       param = list(dh2=list(shape1=20, shape2=25), dP=list(meanlog=-2, sdlog=1), dw2=list(meanlog=5, sdlog=1.5), dNe=list(meanlog=10, sdlog=1), dk=list(lambda=15, kmax=200), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-398.png) ![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-399.png) 

```r
fit.QG <- bayou.mcmc(tree, dat, SE=SE, model="QG", prior=prior.QG, ngen=10000, new.dir=TRUE)
```

```
## gen			lnL			prior			half.life			Vy			K			birth.k			D0.slide			D1.slide			death.k			.h2			.Ne			.P			.theta			U0.slide			U1.slide			U2.slide			.w2			
## 1000			-153.75			-133.49			112.89			1.15			16			0.03			1			0.61			0.93			0.33			0.43			0.38			0.91			1			0.86			1			0.81			
```

```
## 2000			-141.5			-107.39			31.86			0.38			11			0.02			0.95			0.63			0.94			0.37			0.35			0.34			0			0.87			1			0.94			0.86			0.72			
```

```
## 3000			-141.97			-115.8			66.64			0.69			13			0.02			0.95			0.55			0.96			0.35			0.35			0.32			0			0.86			1			0.95			0.92			0.67			
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3910.png) 

```
## 4000			-133.39			-120.98			34.58			0.37			14			0.02			0.93			0.53			0.96			0.34			0.34			0.33			0			0.84			0.99			0.93			0.85			0.61			
```

```
## 5000			-132.75			-103.78			34.43			0.37			10			0.02			0.91			0.51			0.97			0.32			0.32			0.32			0			0.82			0.99			0.88			0.74			0.58			
```

```
## 6000			-123			-114.92			17.07			0.31			12			0.03			0.92			0.49			0.96			0.32			0.33			0.32			0			0.81			0.99			0.88			0.72			0.56			
```

```
## 7000			-129.62			-104.02			37.64			0.37			10			0.02			0.93			0.46			0.94			0.32			0.32			0.31			0			0.8			0.99			0.84			0.69			0.53			
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3911.png) 

```
## 8000			-131.5			-91.26			37.73			0.42			8			0.02			0.93			0.47			0.93			0.3			0.32			0.3			0			0.8			0.99			0.83			0.68			0.52			
```

```
## 9000			-138.36			-78.04			40.97			0.38			6			0.02			0.93			0.45			0.94			0.31			0.32			0.31			0			0.8			0.97			0.84			0.68			0.52			
```

```
## 10000			-129.22			-106.01			29.2			0.29			11			0.02			0.93			0.44			0.93			0.31			0.32			0.31			0			0.8			0.97			0.81			0.71			0.5			
```

```r
chain.QG <- load.bayou(fit.QG, save.Rdata=FALSE, cleanup=FALSE)
chain.QG <- set.burnin(chain.QG, 0.3)
plot(chain.QG)
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3912.png) ![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3913.png) ![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3914.png) ![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3915.png) 

```r
out.QG <- summary(chain.QG)
```

```
## bayou MCMC chain: 10000 generations
## 1000 samples, first 300 samples discarded as burnin
## 
## 
## Summary statistics for parameters:
##                Mean        SD Naive SE Time-series SE Effective Size
## lnL       -131.3957   4.24539 0.160346        1.46310          8.419
## prior     -102.2927  11.16789 0.421805        3.89408          8.225
## h2           0.4899   0.06674 0.002521        0.01032         41.844
## P            8.7589   5.40058 0.203977        3.16860          2.905
## w2          22.1001   6.79279 0.256560        1.06436         40.730
## Ne         259.0653 172.22944 6.505016       98.40297          3.063
## k           10.0471   2.12113 0.080114        0.75186          7.959
## ntheta      11.0471   2.12113 0.080114        0.75186          7.959
## root         3.5689   0.08736 0.003300        0.02466         12.551
## all theta    3.7366   1.45329       NA             NA             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##         pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 408 0.8003              5.2787            0.01564       0.6927
## 414 0.3310              3.6542            0.05406       0.5034
## 112 0.3224              5.2442            0.13870       0.7179
## 128 0.3138              5.7495            0.03006       0.5497
## 266 0.2853              3.0536            0.09111       0.6429
## 185 0.2625              4.7647            0.06854       0.3994
## 95  0.2439              3.1377            0.08176       0.2343
## 103 0.2211              4.7842            0.13542       0.6928
## 4   0.2140              4.3367            0.05545       0.5508
## 57  0.2054              1.9358            0.10002       0.4525
## 45  0.2011              4.8553            0.03623       0.4627
## 335 0.1669              4.2441            0.05691       0.6737
## 438 0.1655              2.7794            0.03033       0.4587
## 104 0.1626              5.1468            0.25094       0.3305
## 69  0.1569              2.6476            0.05379       0.5316
## 98  0.1541              3.6615            0.07427       0.5220
## 200 0.1541              2.2032            0.03900       0.7463
## 52  0.1512              4.3529            0.10189       0.5383
## 96  0.1455              3.3798            0.09976       0.4523
## 50  0.1398              3.8009            0.07325       0.6210
## 83  0.1398              3.3302            0.05991       0.3528
## 27  0.1341              1.3321            0.07844       0.4269
## 154 0.1341              0.6296            0.04201       0.5560
## 66  0.1298              4.8575            0.03744       0.1857
## 229 0.1270              2.1794            0.07742       0.2793
## 291 0.1255              3.7436            0.06933       0.4264
## 29  0.1184              1.8073            0.04313       0.7178
## 347 0.1184              4.8772            0.05367       0.9789
## 372 0.1127              4.9595            0.03532       0.9111
## 345 0.1098              3.8741            0.05433       0.3126
## 276 0.1056              4.1883            0.13502       0.5337
```

```r
phenogram.density(tree, dat, chain=chain.QG, burnin=0.3, pp.cutoff=0.5)
plotSimmap.mcmc(tree, chain.QG, burnin=0.3, circle=TRUE, fsize=0.5)
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-3916.png) 



