# Tutorial for the R package *bayou*
The purpose of *bayou* is to fit Bayesian models of adaptive evolution to phylogenetic comparative data. Specifically, *bayou* provides a flexible framework for fitting multi-optima Ornstein-Uhlenbeck models to phylogenetic comparative data. This tutorial demonstrates some of the options for running *bayou*.

## Reversible-jump MCMC over regime placement
In this example, we will fit a reversible-jump MCMC model to an included dataset (the Chelonia dataset of Jaffe et al. 2012). To start, we will specify a model in which no parameters are fixed. This will estimate the posterior of shift number, location and magnitude as well as all other parameters.

We begin by loading the package and data. We will also assume a constant standard error across all taxa. This can instead be a named vector with species-specific measurement error.

```r
require(bayou)
```

```
## Loading required package: bayou
## Loading required package: ape
## Loading required package: geiger
## Loading required package: phytools
## Loading required package: maps
## 
##  # ATTENTION: maps v3.0 has an updated 'world' map.        #
##  # Many country borders and names have changed since 1990. #
##  # Type '?world' or 'news(package="maps")'. See README_v3. #
## 
## 
## Loading required package: coda
```

```r
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

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 
The figure produced gives a rough visual of the chosen prior distributions. 

### Running the MCMC
Now we are going to run the mcmc. We are going to output files to our working directory. If you want to specify another director, replace "getwd()" with the path of the directory. By default, *bayou outputs to the R temporary directory. We will run a relatively short chain of only 10,000 generations.

```r
par(mfrow=c(2,3))
fit1 <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior, ngen=10000, new.dir=getwd(), plot.freq=2000, ticker.freq=1000)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1000			-125.53			-94.65			15.72			0.26			15			0.38			0.03			1			0.54			0.8			0.33			0.61			1			0.86			0.56			
```

```
## 2000			-123.52			-79.5			19.47			0.3			12			0.35			0.02			1			0.55			0.9			0.32			0.65			0.98			0.67			0.68			
## 3000			-120.82			-74.76			20.93			0.31			11			0.34			0.02			0.95			0.43			0.77			0.31			0.66			0.97			0.54			0.63			
```

```
## 4000			-124.49			-88.66			24.12			0.3			13			0.35			0.02			0.97			0.39			0.77			0.32			0.66			0.97			0.53			0.56			
## 5000			-115.62			-60.63			17.48			0.28			8			0.33			0.02			0.93			0.36			0.76			0.33			0.68			0.95			0.47			0.49			
```

```
## 6000			-113.94			-69.07			10.18			0.23			10			0.32			0.02			0.92			0.35			0.75			0.35			0.65			0.94			0.47			0.45			
## 7000			-115.61			-75.1			20.46			0.3			11			0.32			0.02			0.93			0.36			0.75			0.34			0.63			0.94			0.41			0.43			
```

```
## 8000			-115.48			-84.32			19.43			0.27			13			0.32			0.02			0.93			0.35			0.74			0.32			0.64			0.95			0.39			0.44			
## 9000			-103.94			-84.85			6.37			0.18			13			0.32			0.02			0.93			0.36			0.74			0.31			0.64			0.95			0.41			0.42			
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

```
## 10000			-101.56			-120.75			10.91			0.19			20			0.32			0.02			0.93			0.36			0.72			0.31			0.64			0.95			0.41			0.42			
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
## /home/josef/repos/bayou/SCYQNABVOL/bayou.* 
## To load results, use 'load.bayou(bayouFit)'
## 
## 10000  generations were run with the following acceptance probabilities:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     0.32     0.02     0.93     0.36     0.72     0.31     0.64     0.95 
## U1.slide U2.slide 
##     0.41     0.42 
##  Total number of proposals of each type:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     1801     4470      208      247      121      904     1818      192 
## U1.slide U2.slide 
##      107      132
```

We can load the actual chains by running the following code:

```r
chain <- load.bayou(fit1, save.Rdata=FALSE, cleanup=TRUE)
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
##                    Mean           SD     Naive SE Time-series SE
## lnL       -113.33199782  6.289098163 0.2375359593    2.536990156
## prior      -78.87443068 13.383969008 0.5055055328    4.958035901
## alpha        0.05298113  0.019942887 0.0007532325    0.007365484
## sig2         0.02548311  0.006813172 0.0002573300    0.001870560
## k           11.63908702  2.627244351 0.0992296496    0.977743352
## ntheta      12.63908702  2.627244351 0.0992296496    0.977743352
## root         3.45906124  0.138862071 0.0052447480    0.092966946
## all theta    3.85972343  1.361042120           NA             NA
##           Effective Size
## lnL             6.145245
## prior           7.287029
## alpha           7.331174
## sig2           13.266468
## k               7.220233
## ntheta          7.220233
## root            2.231054
## all theta             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##            pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 408 1.0000000           4.9947711        0.008341552    1.4108448
## 45  0.8815977           4.7492117        0.008426393    0.3397488
## 372 0.7817404           4.4678498        0.015082857    3.0136051
## 369 0.3480742           2.8055353        0.071980790   50.9970587
## 68  0.3281027           5.9475910        0.070164062    0.1042369
## 309 0.2796006           3.3397527        0.052709755    1.3202665
## 354 0.2282454           6.2253372        0.081139679    2.6090662
## 371 0.2282454           2.3537493        0.033680473   44.7883898
## 272 0.2182596           2.7595869        0.055108280    1.5059489
## 13  0.2111270           6.7840454        0.094898801    1.3370873
## 126 0.2068474           2.2358429        0.064994897    2.4649247
## 361 0.2068474           5.2429577        0.055309867   33.6945292
## 232 0.2011412           5.1732440        0.102656563    2.3930744
## 6   0.1997147           4.1778471        0.033393608    0.6757617
## 355 0.1654779           5.7713349        0.103384498    9.2292230
## 99  0.1640514           2.7083938        0.144858201    0.1787315
## 411 0.1626248           3.1496650        0.002056607    1.3398900
## 323 0.1554922           2.8572982        0.037390134    0.7247873
## 221 0.1469330           3.5745479        0.046490003    1.4317451
## 144 0.1455064           2.1073281        0.025260435    3.6989363
## 373 0.1426534           4.0804670        0.004812412    2.1976420
## 291 0.1412268           3.8948639        0.034436891   11.5895658
## 31  0.1369472           3.6580475        0.053312372    0.7661552
## 360 0.1369472           5.2044572        0.069612570   10.2968016
## 47  0.1340942           3.8178294        0.017288693    7.3322474
## 324 0.1340942           3.4189660        0.042819305    4.3501053
## 57  0.1312411           3.6610952        0.053943142    1.2088844
## 100 0.1298146           3.2230482        0.084184088    0.3390053
## 427 0.1212553           2.7410386        0.024575094    1.1181451
## 441 0.1098431           3.1857820        0.005805380    0.7178953
## 306 0.1055635           1.4392780        0.052907009    0.8005667
## 106 0.1041369           0.9331014        0.037549577    0.3712552
```

We can visualize the chains by viewing traces for each parameter using the plotting utilities of the R package coda:

```r
plot(chain)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png) 

We can view where there are shifts of high probability.

```r
par(mfrow=c(1,1))
plotSimmap.mcmc(chain, burnin=0.3, lwd=2, edge.type="theta", pal=colorRampPalette(c("black", "gray90")), show.tip.label=FALSE)
```

```
## Warning in plotSimmap.mcmc(chain, burnin = 0.3, lwd = 2, edge.type =
## "theta", : Length of post-burnin sample less than the requested parameter
## sample, using entire post-burnin chain instead
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

And we can view the density of phenotypic optima and location of highly supported shifts on the phenogram. Here we show all shifts with posterior probabilities greater than *pp.cutoff = 0.3*. 

```r
phenogram.density(tree, dat, chain=chain, burnin=0.3, pp.cutoff=0.3)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

### Diagnosing Convergence
It is useful to run two chains with independent starting positions. By default, *bayou* simulates starting parameters from the prior distribution. We will re-run the analysis here to obtain an independent chain.

```r
fit2 <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior, ngen=10000, new.dir=getwd(), plot.freq=NULL, ticker.freq=1000)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1000			-149.17			-85.47			2.07			0.27			13			0.38			0.03			0.59			0.27			0.85			0.39			0.56			1			0.5			0.4			
## 2000			-95.94			-88.94			6.31			0.18			14			0.34			0.02			0.62			0.36			0.63			0.35			0.52			1			0.33			0.3			
## 3000			-99.67			-90.56			6.66			0.22			14			0.32			0.02			0.67			0.33			0.6			0.29			0.53			0.98			0.32			0.2			
## 4000			-103.42			-80.13			6.56			0.22			12			0.32			0.02			0.71			0.29			0.65			0.31			0.52			0.96			0.33			0.17			
## 5000			-112.44			-89.45			15.71			0.25			14			0.32			0.02			0.74			0.3			0.64			0.29			0.53			0.94			0.36			0.17			
## 6000			-115.16			-94.94			15.88			0.24			15			0.32			0.02			0.75			0.31			0.65			0.3			0.56			0.95			0.4			0.21			
## 7000			-105.54			-79.19			10.59			0.19			12			0.32			0.02			0.77			0.32			0.67			0.31			0.56			0.95			0.36			0.24			
## 8000			-105.1			-84.14			4.87			0.18			13			0.33			0.02			0.79			0.32			0.68			0.31			0.56			0.96			0.39			0.22			
## 9000			-111			-64.09			10.45			0.25			9			0.34			0.02			0.79			0.32			0.68			0.31			0.55			0.96			0.39			0.21			
## 10000			-110.37			-70.06			6.88			0.21			10			0.34			0.02			0.77			0.31			0.65			0.32			0.54			0.97			0.37			0.21			
```

```r
chain2 <- load.bayou(fit2, save.Rdata=FALSE, cleanup=FALSE)
chain2 <- set.burnin(chain2, 0.3)
```

Now we can compare chains to see if the parameters have converged using Gelman and Rubin's R statistic. Values close to 1 indicate chains have converged (these two chains will not have converged in only 10,000 generations).

```r
RlnL <- gelman.R("lnL", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

```r
Ralpha <- gelman.R("alpha", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-2.png) 

```r
Rsig2 <- gelman.R("sig2", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-3.png) 

Of particular interest to inference are the branch posterior probabilities, we can plot these to see if the two chains converged on similar answers. If the runs have converged, points should fall along the y=x line.

```r
L1 <- Lposterior(chain,tree, burnin=0.3)
L2 <- Lposterior(chain2,tree, burnin=0.3)
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 

We can also combine chains:

```r
chains <- combine.chains(chain, chain2, burnin.prop=0.3)
chains <- set.burnin(chains, 0)
```

### Estimating marginal likelihoods
To compare models, we want to estimate the marginal likelihood, which is performed using the stepping stone algorithm in *bayou*. The method first estimates a reference function by fitting a series of curves to the posterior of an MCMC chain. The function *steppingstone* will output a graphic showing the best-fitting density functions overlaying the posterior distribution. Then the stepping stone MCMC's are run for every value in the vector *Bk*, corresponding to each value of the power posterior function ranging from the reference function (*Bk = 0*) to the posterior (*Bk = 1*). 
To speed things up, you can specify multiple cores. Default is set to 2 cores, but if you have more you can use them.

```r
ss <- steppingstone(Bk=seq(0,1,length.out=5), chains, tree, dat, SE=SE, prior=prior, new.dir=getwd(), ngen=10000)
```

```
## Making power posterior function from provided mcmc chain...
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 

```
## Running mcmc chains...
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0			1000			-1224.05			-70.59			-64.21			10.16			0.34			10			0.67			0.01			1			1			1			0.67			0.8			1			1			1			
## 0			2000			-144.12			-37.94			-28.99			27.39			0.52			4			0.69			0.01			1			1			1			0.67			0.79			1			1			1			
## 0			3000			-164.28			-43.42			-32.27			12.58			0.26			5			0.7			0.01			1			1			1			0.67			0.78			1			1			1			
## 0			4000			-273.39			-48.39			-38.13			7.98			0.2			6			0.72			0.01			1			1			1			0.66			0.78			1			1			1			
## 0			5000			-309.91			-42.89			-31.38			8.79			0.21			5			0.71			0.01			1			1			1			0.67			0.79			1			1			1			
## 0			6000			-166.03			-32.94			-24.53			19.41			0.26			3			0.72			0.01			1			1			1			0.67			0.79			1			1			1			
## 0			7000			-768.61			-32.83			-22.8			12.8			0.2			3			0.72			0.01			1			1			1			0.68			0.8			1			1			1			
## 0			8000			-139.05			-37.59			-28			26.21			0.45			4			0.73			0.01			1			1			1			0.67			0.8			1			1			1			
## 0			9000			-240.26			-17.53			-8.18			8.1			0.17			0			0.73			0.01			1			1			1			0.67			0.8			1			1			1			
## 0			10000			-260.72			-48.34			-37.09			12.86			0.29			6			0.72			0.01			1			1			1			0.68			0.81			1			1			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.25			1000			-151.93			-66.43			-59.87			21.79			0.34			9			0.54			0			1			0.78			1			0.59			0.8			1			0.92			0.82			
## 0.25			2000			-143.3			-80.03			-71			21.94			0.4			12			0.53			0.01			0.98			0.66			1			0.56			0.79			1			0.95			0.89			
## 0.25			3000			-142.16			-53.66			-45.23			22.45			0.45			7			0.55			0.01			0.96			0.65			1			0.53			0.77			1			0.93			0.92			
## 0.25			4000			-155.59			-42.67			-31.28			16.5			0.26			5			0.55			0.01			0.96			0.67			1			0.55			0.74			0.99			0.94			0.94			
## 0.25			5000			-144.68			-48.19			-38.92			19.86			0.46			6			0.55			0.01			0.97			0.72			1			0.56			0.73			0.99			0.91			0.95			
## 0.25			6000			-142.9			-42.85			-33.53			29.79			0.51			5			0.55			0.01			0.97			0.77			1			0.55			0.72			0.99			0.92			0.95			
## 0.25			7000			-140.24			-43.47			-38.03			38.31			0.52			5			0.55			0.01			0.97			0.79			1			0.55			0.72			0.99			0.92			0.95			
## 0.25			8000			-142.23			-43.61			-33.95			22.54			0.46			5			0.54			0.01			0.97			0.79			1			0.54			0.71			0.99			0.92			0.96			
## 0.25			9000			-140.53			-32.6			-23.9			27.77			0.46			3			0.55			0.01			0.97			0.8			1			0.55			0.71			1			0.93			0.96			
## 0.25			10000			-145.18			-32.76			-24.99			34.48			0.62			3			0.54			0.01			0.96			0.79			1			0.55			0.7			1			0.93			0.97			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.5			1000			-138.61			-69.31			-60.82			22.72			0.38			10			0.48			0.02			0.92			0.7			1			0.58			0.71			0.94			0.93			1			
## 0.5			2000			-142.23			-43.51			-37.52			33.15			0.49			5			0.48			0.01			0.88			0.8			1			0.5			0.74			0.98			0.95			1			
## 0.5			3000			-139.42			-27.73			-20.49			31.67			0.56			2			0.47			0.01			0.91			0.82			1			0.45			0.73			0.96			0.78			0.89			
## 0.5			4000			-141.44			-38.03			-31.67			36.88			0.6			4			0.48			0.01			0.93			0.73			1			0.43			0.72			0.97			0.83			0.87			
## 0.5			5000			-140.47			-27.83			-26.43			45.65			0.47			2			0.49			0.01			0.94			0.72			1			0.43			0.72			0.98			0.84			0.88			
## 0.5			6000			-140.44			-48.37			-40.15			31.87			0.53			6			0.49			0.01			0.95			0.69			1			0.41			0.71			0.98			0.85			0.9			
## 0.5			7000			-132.8			-49.83			-42.83			26.09			0.4			6			0.5			0.01			0.95			0.71			0.98			0.41			0.71			0.98			0.81			0.84			
## 0.5			8000			-132.24			-53.78			-43.67			18.97			0.32			7			0.49			0.01			0.91			0.7			0.98			0.41			0.71			0.98			0.8			0.82			
## 0.5			9000			-142.7			-32.86			-22.22			18.51			0.36			3			0.48			0.01			0.92			0.7			0.98			0.42			0.7			0.98			0.81			0.83			
## 0.5			10000			-138.7			-47.96			-37.06			22.78			0.42			6			0.48			0.01			0.93			0.7			0.98			0.42			0.7			0.99			0.79			0.82			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.75			1000			-142.72			-58.67			-49.34			25.93			0.47			8			0.45			0.02			1			0.66			1			0.46			0.8			1			0.67			0.6			
## 0.75			2000			-134.66			-70.26			-64.21			27.21			0.38			10			0.41			0.02			0.97			0.69			1			0.42			0.78			1			0.87			0.7			
## 0.75			3000			-127.63			-69.07			-59.12			22.65			0.33			10			0.4			0.02			0.98			0.71			1			0.41			0.77			1			0.75			0.73			
## 0.75			4000			-133.33			-44.34			-40.12			33.3			0.39			5			0.41			0.02			0.97			0.67			0.97			0.41			0.77			0.99			0.67			0.73			
## 0.75			5000			-140.64			-43.25			-33.07			19.03			0.42			5			0.42			0.02			0.98			0.65			0.97			0.4			0.77			0.99			0.73			0.75			
## 0.75			6000			-132.27			-64.33			-57.39			21.07			0.31			9			0.42			0.02			0.98			0.65			0.98			0.39			0.76			0.99			0.73			0.73			
## 0.75			7000			-142.23			-48.53			-40.35			20.86			0.33			6			0.42			0.02			0.97			0.65			0.98			0.39			0.76			0.99			0.74			0.71			
## 0.75			8000			-136.54			-66.42			-64.03			32.4			0.46			9			0.43			0.02			0.97			0.64			0.98			0.39			0.76			0.99			0.74			0.74			
## 0.75			9000			-132.39			-69.42			-63.38			30.69			0.42			10			0.43			0.02			0.97			0.62			0.98			0.39			0.76			0.99			0.77			0.76			
## 0.75			10000			-136.29			-63.84			-58.02			33.74			0.48			9			0.43			0.02			0.97			0.61			0.99			0.38			0.76			0.98			0.76			0.74			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1			1000			-137.31			-78.54			-76.79			39.49			0.53			11			0.46			0.03			0.95			0.69			1			0.49			0.85			1			0.87			0.8			
## 1			2000			-117.67			-93.93			-81.47			17.91			0.25			15			0.44			0.03			0.93			0.6			1			0.42			0.81			0.94			0.81			0.63			
## 1			3000			-114.26			-89.96			-80.97			19.42			0.25			14			0.39			0.03			0.9			0.49			0.9			0.41			0.76			0.93			0.6			0.52			
## 1			4000			-103.91			-68.65			-55.65			12.08			0.21			10			0.37			0.02			0.91			0.47			0.86			0.4			0.71			0.94			0.44			0.43			
## 1			5000			-111.11			-88.62			-75.44			16.57			0.23			14			0.36			0.02			0.92			0.44			0.81			0.4			0.66			0.94			0.38			0.4			
## 1			6000			-108.86			-89.71			-79.07			13.32			0.22			14			0.36			0.02			0.91			0.41			0.81			0.39			0.65			0.93			0.35			0.37			
## 1			7000			-120.83			-83.75			-69.65			15.52			0.33			13			0.36			0.02			0.9			0.42			0.81			0.39			0.66			0.94			0.36			0.38			
## 1			8000			-112.54			-63.66			-51.49			12.65			0.23			9			0.36			0.02			0.87			0.41			0.75			0.37			0.65			0.95			0.35			0.4			
## 1			9000			-110.95			-48.37			-36.72			10.25			0.23			6			0.37			0.02			0.88			0.42			0.74			0.35			0.62			0.95			0.35			0.4			
## 1			10000			-101.2			-84.52			-73.73			7.18			0.19			13			0.37			0.02			0.88			0.41			0.73			0.36			0.61			0.94			0.38			0.4			
## Loading mcmc chains...
```

```r
ss
```

```
## Stepping stone estimation of marginal likelihood
## Marginal Likelihood:
## [1] -145.9177
## A total of 5 power posteriors were run along the sequence: 0		0.25		0.5		0.75		1
## lnr_k -39.89 -36.87 -35.4 -33.75
```

```r
ss <- set.burnin(ss, 0.3)
```

Plotting the stepping stone object will give you an idea of convergence for each of the MCMC's by providing a trace of the lnL, ln prior, and reference function as well as the estimated marginal likelihood at each step.

```r
plot(ss)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png) 

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
plotBayoupars(startpar, tree, cex=0.5)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png) 

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
plotBayoupars(converge_pars, tree, cex=0.5)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png) 

Now that we know how to specify fixed models, lets set up the prior function:

```r
prior.fixed <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="fixed", dk="fixed", dtheta="dnorm", dloc="dloc"), 
                            param=list(dalpha=list(scale=1), dsig2=list(scale=1), dtheta=list(mean=mean(dat), sd=2)),
                                     fixed=fixed.hypothesis)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png) 

A few things to note. First, we have set *dsb* and *dk* to "fixed". We could have also set *dloc*, *dalpha* and/or *dsig2* to "fixed". Here, we will let the location of the shift on the branch vary (*dloc* puts a uniform prior on the location of the shift on the branch). Second, we are passing our model through the *fixed.hypothesis* object we created using *identifyBranches()*. For example, since simple Brownian motion is a special case of multi-optima OU models, we can specify through the the following prior:

```r
prior.BM <- make.prior(tree, dists=list(dalpha="fixed", dsig2="dhalfcauchy", dsb="fixed", dk="fixed", dloc="fixed"),
                            param=list(dsig2=list(scale=1), dtheta=list(mean=mean(dat), sd=2)),
                                fixed=list(alpha=0, k=0, sb=numeric(0), loc=numeric(0), t2=numeric(0)))
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png) 

Now we can run our MCMC chains using our priors. Right now, you have to specify the starting parameters. This will soon be unnecessary...

```r
fit.fixed <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior=prior.fixed, startpar=startpar, ngen=10000, new.dir=getwd(), plot.freq=NULL)
```

```
## gen			lnL			prior			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 1000			-138.77			-5.76			29.52			0.42			2			0.41			1			0.39			0.34			1			
## 2000			-142.96			-5.81			42.26			0.44			2			0.43			0.98			0.35			0.39			1			
## 3000			-140.13			-5.78			31.99			0.5			2			0.43			0.97			0.35			0.4			1			
## 4000			-143.94			-5.98			39.94			0.43			2			0.44			0.97			0.35			0.4			1			
## 5000			-141.29			-5.87			50.82			0.61			2			0.45			0.96			0.36			0.42			1			
## 6000			-141.58			-5.86			42.75			0.51			2			0.45			0.96			0.35			0.42			1			
## 7000			-140.81			-5.75			45.74			0.55			2			0.45			0.97			0.35			0.41			1			
## 8000			-139.25			-5.8			41.98			0.49			2			0.45			0.96			0.35			0.41			1			
## 9000			-138.91			-5.75			35.04			0.46			2			0.45			0.97			0.35			0.41			1			
## 10000			-139.94			-5.79			26.92			0.38			2			0.45			0.97			0.34			0.41			1			
```




```r
chain.fixed <- load.bayou(fit.fixed, save.Rdata=FALSE, cleanup=FALSE)
chain.fixed <- set.burnin(chain.fixed, 0.3)
plot(chain.fixed)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-1.png) ![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-2.png) 

```r
out.fixed <- summary(chain.fixed)
```

```
## bayou MCMC chain: 10000 generations
## 1000 samples, first 300 samples discarded as burnin
## 
## 
## Summary statistics for parameters:
##                    Mean          SD     Naive SE Time-series SE
## lnL       -140.80255750 1.883634869 0.0711439071    0.207740763
## prior       -5.80331910 0.080421745 0.0030374874    0.009848098
## alpha        0.01923291 0.004933909 0.0001863512    0.000606389
## sig2         0.01815236 0.002846870 0.0001075248    0.000309264
## k            2.00000000 0.000000000 0.0000000000    0.000000000
## ntheta       3.00000000 0.000000000 0.0000000000    0.000000000
## root         3.70446204 0.181409093 0.0068517269    0.013179465
## all theta    3.62297756 0.385489234           NA             NA
##           Effective Size
## lnL             82.21480
## prior           66.68717
## alpha           66.20330
## sig2            84.73770
## k                0.00000
## ntheta           0.00000
## root           189.46220
## all theta             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##     pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 409  1            3.658931         0.01292093     1.530067
## 411  1            3.505539         0.01978969     2.148033
```

```r
phenogram.density(tree, dat, chain=chain.fixed, burnin=0.3, pp.cutoff=0.5)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-3.png) 

```r
plotSimmap.mcmc(chain.fixed, burnin=0.3, lwd=2, edge.type="theta", pal=colorRampPalette(c("black", "gray90")), show.tip.label=FALSE)
```

```
## Warning in plotSimmap.mcmc(chain.fixed, burnin = 0.3, lwd = 2, edge.type =
## "theta", : Length of post-burnin sample less than the requested parameter
## sample, using entire post-burnin chain instead
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-4.png) 




As before, we can estimate the marginal likelihood:

```r
ss.fixed <- steppingstone(Bk=seq(0,1,length.out=5), chain.fixed, tree, dat, SE=0, startpar=startpar, prior=prior.fixed, ngen=10000)
```

```
## Making power posterior function from provided mcmc chain...
```

```
## Running mcmc chains...
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0			1000			-150.42			-5.76			6.82			29.79			0.51			2			0.64			1			0.45			0.48			1			
## 0			2000			-158.38			-5.76			9.46			51.34			0.64			2			0.61			1			0.46			0.47			1			
## 0			3000			-163.59			-5.74			8.46			27			0.33			2			0.63			1			0.45			0.45			1			
## 0			4000			-171.63			-5.82			5.35			44.91			0.46			2			0.63			1			0.46			0.47			1			
## 0			5000			-171.51			-5.92			3.6			45.9			0.45			2			0.64			1			0.45			0.47			1			
## 0			6000			-165.5			-5.82			6.65			27.35			0.35			2			0.64			1			0.46			0.47			1			
## 0			7000			-156.27			-5.75			7.79			32.09			0.51			2			0.64			1			0.47			0.47			1			
## 0			8000			-169.43			-5.93			5.07			49.52			0.56			2			0.64			1			0.47			0.47			1			
## 0			9000			-167.9			-5.75			8.29			48.5			0.47			2			0.64			1			0.47			0.47			1			
## 0			10000			-160.09			-5.79			8.04			61.37			0.77			2			0.64			1			0.47			0.48			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0.25			1000			-152.27			-5.76			4.9			24.27			0.43			2			0.59			0.99			0.4			0.5			1			
## 0.25			2000			-156.21			-5.8			8.48			32.84			0.45			2			0.59			0.99			0.4			0.52			1			
## 0.25			3000			-155.8			-5.8			7.98			38.49			0.54			2			0.59			0.99			0.39			0.51			1			
## 0.25			4000			-159.03			-5.86			5.82			56.48			0.8			2			0.58			0.99			0.38			0.49			1			
## 0.25			5000			-155.48			-5.76			6.88			44.25			0.79			2			0.58			1			0.37			0.48			1			
## 0.25			6000			-150.57			-5.77			7.52			28.45			0.46			2			0.58			0.99			0.37			0.47			1			
## 0.25			7000			-155.06			-5.93			1.09			33.56			0.62			2			0.57			0.99			0.39			0.47			1			
## 0.25			8000			-154.43			-5.75			9.13			41.01			0.56			2			0.57			0.99			0.39			0.47			0.99			
## 0.25			9000			-160.33			-5.81			8.65			34.79			0.42			2			0.57			0.99			0.39			0.47			0.99			
## 0.25			10000			-159.94			-5.79			6.08			62.6			0.97			2			0.57			1			0.39			0.47			0.99			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0.5			1000			-151.69			-5.78			7.18			33.46			0.58			2			0.5			0.96			0.35			0.42			1			
## 0.5			2000			-152.04			-5.87			3.29			24.3			0.42			2			0.47			0.98			0.36			0.38			1			
## 0.5			3000			-154.45			-5.75			7.07			42.2			0.73			2			0.48			0.97			0.37			0.38			1			
## 0.5			4000			-150.9			-5.8			4.52			29.89			0.53			2			0.5			0.97			0.37			0.39			0.99			
## 0.5			5000			-158.36			-5.83			7.03			52.89			0.73			2			0.51			0.98			0.37			0.4			0.99			
## 0.5			6000			-155.2			-5.76			9.52			32.93			0.42			2			0.52			0.98			0.37			0.4			0.99			
## 0.5			7000			-152.91			-5.77			8.27			36.61			0.52			2			0.52			0.98			0.37			0.4			0.99			
## 0.5			8000			-157.17			-5.77			9.47			44.92			0.55			2			0.52			0.98			0.36			0.4			0.99			
## 0.5			9000			-153.72			-5.77			9.46			36.99			0.51			2			0.52			0.98			0.36			0.4			0.99			
## 0.5			10000			-157.07			-5.92			5.11			40.76			0.63			2			0.52			0.98			0.36			0.4			0.99			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 0.75			1000			-154.67			-5.79			8.09			33.8			0.5			2			0.45			1			0.4			0.43			1			
## 0.75			2000			-150.39			-5.75			6.55			29.04			0.51			2			0.45			0.98			0.36			0.42			1			
## 0.75			3000			-149.29			-5.76			0.75			19.15			0.38			2			0.44			0.98			0.37			0.4			0.99			
## 0.75			4000			-153.4			-5.79			6.42			23.36			0.36			2			0.45			0.98			0.36			0.4			0.99			
## 0.75			5000			-152.57			-5.8			5.12			30.6			0.57			2			0.45			0.98			0.35			0.4			0.99			
## 0.75			6000			-154.06			-5.79			8.52			38.34			0.56			2			0.45			0.98			0.34			0.4			0.99			
## 0.75			7000			-149.74			-5.75			5.64			27.48			0.5			2			0.45			0.99			0.34			0.39			0.99			
## 0.75			8000			-153.1			-5.74			8.35			33.12			0.5			2			0.45			0.99			0.33			0.39			1			
## 0.75			9000			-154.88			-5.75			9.77			44.51			0.61			2			0.45			0.99			0.32			0.39			0.99			
## 0.75			10000			-152.93			-5.81			5.45			31.37			0.57			2			0.44			0.99			0.32			0.39			0.99			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			D0.slide			.sig2			.theta			U0.slide			
## 1			1000			-149.87			-5.76			5.3			27.67			0.49			2			0.37			1			0.35			0.32			1			
## 1			2000			-151.8			-5.81			3.7			24.52			0.46			2			0.37			0.99			0.32			0.34			1			
## 1			3000			-148.37			-5.75			0.94			21.28			0.46			2			0.36			0.99			0.31			0.32			1			
## 1			4000			-153.67			-5.92			4.1			33.86			0.59			2			0.36			0.99			0.31			0.33			1			
## 1			5000			-150.36			-5.8			5.74			26.18			0.45			2			0.37			0.99			0.28			0.32			0.99			
## 1			6000			-151			-5.77			-6.51			19.23			0.5			2			0.38			0.99			0.28			0.32			1			
## 1			7000			-152.07			-5.8			-9.94			16.69			0.45			2			0.38			0.98			0.29			0.33			1			
## 1			8000			-148.83			-5.76			0.98			19.92			0.41			2			0.38			0.99			0.3			0.33			1			
## 1			9000			-149.67			-5.77			4.41			22.98			0.43			2			0.39			0.98			0.3			0.33			1			
## 1			10000			-148.51			-5.75			-1.66			19.07			0.43			2			0.39			0.98			0.3			0.33			1			
## Loading mcmc chains...
```

```r
ss.fixed <- set.burnin(ss.fixed, 0.3)
ss.fixed$lnr
```

```
## [1] -163.7467
```

```r
plot(ss.fixed)
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28-1.png) ![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28-2.png) 

And calculate Bayes Factors in support of the Reversible-Jump model vs. the fixed model. 

```r
2*(ss$lnr-ss.fixed$lnr)
```

```
## [1] 35.65802
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
##       -lnL      AIC     AICc model ntax
##  -148.0593 306.1186 306.3913   OUM  226
## 
## 
## Rates
##                   1          2          3
## alpha    0.03232904 0.03232904 0.03232904
## sigma.sq 0.02785184 0.02785184 0.02785184
## 
## Optima
##                  1         2         3
## estimate 3.7433452 3.4768778 3.5108183
## se       0.1221091 0.3873866 0.1312257
## 
## Arrived at a reliable solution
```

```r
out.fixed$statistics
```

```
##                    Mean          SD     Naive SE Time-series SE
## lnL       -140.80255750 1.883634869 0.0711439071    0.207740763
## prior       -5.80331910 0.080421745 0.0030374874    0.009848098
## alpha        0.01923291 0.004933909 0.0001863512    0.000606389
## sig2         0.01815236 0.002846870 0.0001075248    0.000309264
## k            2.00000000 0.000000000 0.0000000000    0.000000000
## ntheta       3.00000000 0.000000000 0.0000000000    0.000000000
## root         3.70446204 0.181409093 0.0068517269    0.013179465
## all theta    3.62297756 0.385489234           NA             NA
##           Effective Size
## lnL             82.21480
## prior           66.68717
## alpha           66.20330
## sig2            84.73770
## k                0.00000
## ntheta           0.00000
## root           189.46220
## all theta             NA
```

Of course, we can also go the other direction, and convert back to *bayou* format:

```r
back_pars <- OUwie2bayou(OUwieData$tree, OUwieData$dat)
plotBayoupars(back_pars, tree)
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32-1.png) 

## Some other useful things
### Some plotting features
It's easy to convert *bayou* formatted models into trees that can be visualized using *phytools*' *plotSimmap* or *phenogram* functions. For example, let's pull a few samples from the posterior out of our original chain:

```r
samp_pars <- pull.pars(500, chain, model="OU")
samp_pars
```

```
## $alpha
## [1] 0.03965192
## 
## $sig2
## [1] 0.02253413
## 
## $k
## [1] 8
## 
## $ntheta
## [1] 9
## 
## $theta
## [1] 3.511195 5.261785 4.510632 5.472907 6.921682 4.075489 4.802339 3.105937
## [9] 1.976940
## 
## $sb
## [1] 408  45 355  68 232 372  77 107
## 
## $loc
## [1]  0.42805112  0.54618188  4.86392249  1.44505267 19.17510788  0.45299981
## [7]  0.07134491  0.07874794
## 
## $t2
## [1] 2 3 4 5 6 7 8 9
```

Now let's convert it to *phytools*' simmap format:

```r
sm_tree <- pars2simmap(samp_pars,tree)
plotSimmap(sm_tree$tree, col=sm_tree$col)
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-34-1.png) 

```r
phenogram(sm_tree$tree, dat, col=sm_tree$col, ftype="off", spread.labels=FALSE)
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-34-2.png) 

We can also make some regime plots, which show the location of each optimum and the expected stationary variance.

```r
plot(c(0, 210), c(2, 6), type="n", xlab="time", ylab="phenotype")
regime.plot(samp_pars, sm_tree$tree, type="density", cols=sm_tree$col)
```

```
## Error in eval(expr, envir, enclos): could not find function "regime.plot"
```

```r
phenogram(sm_tree$tree, dat, col=sm_tree$col, ftype="off", spread.labels=FALSE, add=TRUE)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35-1.png) 

....Put in how to do OU ancestral state reconstructions here....


## How to specify some other models
Rather than restricting our interest to only models with a single shift allowed per branch, we can fit a number of models that allow multiple shifts per branch, have unequal probabilities, etc. These are all specified with the prior function. Here are some examples:

A model with an arbitrary number of shifts allowed per branch, and branches chosen with probability proportional to their length:

```r
tree <- reorder(tree, "postorder")
dat <- dat[tree$tip.label]
prior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=Inf,prob=tree$edge.length), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-36-1.png) 

A model that disallows shifts on terminal branches:

```r
terminal.branches <- as.numeric(tree$edge[,2] < length(tree$tip.label)+1)
prior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=terminal.branches, prob=1), dtheta=list(mean=mean(dat), sd=2)))
```

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-37-1.png) 


...Add QG, OUrepar models here...


