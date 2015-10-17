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
## 1000			-127.15			-131.41			23.83			0.36			22			0.39			0.04			0.96			0.57			0.91			0.4			0.72			1			0.75			1			
```

```
## 2000			-122.87			-122.92			25.64			0.35			20			0.39			0.04			0.96			0.55			0.87			0.39			0.72			0.98			0.72			0.85			
## 3000			-115.56			-87.19			20.51			0.28			13			0.36			0.03			0.91			0.45			0.87			0.36			0.7			0.97			0.64			0.57			
```

```
## 4000			-112.89			-97.34			12.09			0.23			15			0.34			0.03			0.88			0.43			0.82			0.33			0.69			0.94			0.63			0.47			
## 5000			-92.88			-84.19			3.47			0.17			13			0.34			0.03			0.88			0.38			0.71			0.33			0.65			0.93			0.57			0.35			
```

```
## 6000			-95.13			-63.35			5.45			0.19			9			0.33			0.02			0.88			0.37			0.68			0.32			0.58			0.92			0.51			0.34			
## 7000			-94.52			-89.24			3.36			0.15			14			0.32			0.02			0.88			0.37			0.65			0.31			0.56			0.92			0.47			0.31			
```

```
## 8000			-100.38			-78.37			5.15			0.18			12			0.32			0.02			0.88			0.36			0.64			0.3			0.54			0.91			0.41			0.28			
## 9000			-99.33			-101.64			5.33			0.17			16			0.32			0.02			0.9			0.36			0.61			0.3			0.52			0.91			0.38			0.26			
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

```
## 10000			-96.05			-105.73			7.47			0.17			17			0.32			0.02			0.89			0.36			0.6			0.31			0.52			0.92			0.37			0.26			
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
## /home/josef/repos/bayou/TIGOSCKUHQ/bayou.* 
## To load results, use 'load.bayou(bayouFit)'
## 
## 10000  generations were run with the following acceptance probabilities:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     0.32     0.02     0.89     0.36     0.60     0.31     0.52     0.92 
## U1.slide U2.slide 
##     0.37     0.26 
##  Total number of proposals of each type:
##   .alpha  birth.k D0.slide D1.slide  death.k    .sig2   .theta U0.slide 
##     1843     4394      214      298      154      918     1733      179 
## U1.slide U2.slide 
##      140      127
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
##                   Mean         SD     Naive SE Time-series SE
## lnL       -98.05993062 8.84059387 0.3339046222    4.201966714
## prior     -87.25864762 9.29079604 0.3509085234    2.951837097
## alpha       0.13538866 0.06421844 0.0024254971    0.029744210
## sig2        0.04667301 0.01726971 0.0006522681    0.005863489
## k          13.48787447 1.75545962 0.0663027948    0.548204030
## ntheta     14.48787447 1.75545962 0.0663027948    0.548204030
## root        3.81809438 0.10409297 0.0039315372    0.039602598
## all theta   3.83687413 1.06778927           NA             NA
##           Effective Size
## lnL             4.426471
## prior           9.906519
## alpha           4.661382
## sig2            8.674772
## k              10.254093
## ntheta         10.254093
## root            6.908686
## all theta             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##            pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 372 1.0000000            4.341598        0.009389245   1.74699306
## 408 1.0000000            4.784165        0.005593138   1.54791526
## 252 0.8544936            4.063965        0.008918080   1.10638733
## 411 0.8245364            3.196444        0.003664430   2.76477671
## 43  0.6191155            4.922585        0.011481126   0.18772610
## 447 0.4964337            3.181774        0.010084098   1.72385808
## 316 0.4550642            2.542532        0.011416657   0.55026495
## 300 0.4079886            4.056284        0.023423004   0.28486991
## 45  0.3737518            4.922602        0.020314293   0.50060100
## 38  0.3537803            1.538912        0.075036144   0.06617060
## 444 0.2967190            3.296101        0.019580102   0.93632249
## 292 0.2853067            4.466551        0.030641533  18.64408483
## 440 0.2653352            2.936161        0.021351139   0.11620704
## 74  0.2410842            4.726965        0.066385532   0.02876886
## 402 0.2282454            4.670124        0.014033745   0.35279129
## 232 0.2211127            4.777008        0.048418261   2.29792701
## 1   0.1825963            3.196189        0.029903884   0.40379285
## 79  0.1825963            2.950661        0.238835072   0.04224958
## 306 0.1825963            2.398368        0.026053528   0.59511140
## 224 0.1797432            4.113987        0.076110746   1.64415549
## 409 0.1754636            3.126799        0.006293935   0.63002164
## 369 0.1669044            3.613279        0.059612264  73.57378274
## 294 0.1626248            5.710486        0.039749525   1.28785483
## 249 0.1554922            3.343114        0.068220006   1.81728662
## 219 0.1455064            2.945274        0.000000000   0.43802732
## 205 0.1426534            5.345859        0.070637701   0.26520713
## 99  0.1412268            4.845992        0.022030032   0.47704137
## 323 0.1269615            2.845389        0.022750747   1.08041135
## 304 0.1226819            1.857859        0.048737061   0.23549460
## 358 0.1226819            4.471553        0.041324821   2.80043443
## 293 0.1198288            2.817435        0.047507971  28.13344398
## 305 0.1198288            3.310183        0.065021663   0.20244831
## 359 0.1169757            3.370916        0.046094457  14.84878305
## 360 0.1141227            4.736807        0.073219161   4.81832529
## 442 0.1112696            3.687124        0.037485580   1.91125179
## 446 0.1112696            3.088867        0.015153627   0.31920117
## 403 0.1084165            4.151132        0.049782686   1.65001833
## 340 0.1041369            4.404482        0.079771446   2.28552331
## 200 0.1027104            3.558359        0.058252255   0.59785771
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
## 1000			-116.96			-102.03			18.55			0.32			16			0.43			0.03			0.88			0.64			0.76			0.4			0.69			0.93			0.57			0.62			
## 2000			-118.19			-89.95			14.47			0.29			14			0.39			0.02			0.83			0.56			0.67			0.36			0.72			0.89			0.58			0.55			
## 3000			-110.47			-90			7.6			0.2			14			0.37			0.02			0.82			0.49			0.66			0.34			0.7			0.91			0.55			0.44			
## 4000			-106.05			-105.48			11.51			0.21			17			0.38			0.02			0.83			0.44			0.65			0.32			0.67			0.94			0.51			0.39			
## 5000			-113.26			-102.18			13.36			0.22			16			0.39			0.02			0.83			0.4			0.67			0.31			0.66			0.95			0.5			0.38			
## 6000			-114.44			-89.96			15.31			0.22			14			0.38			0.02			0.84			0.38			0.66			0.31			0.66			0.96			0.46			0.35			
## 7000			-114.17			-94.8			12.98			0.23			15			0.39			0.02			0.84			0.36			0.69			0.31			0.67			0.96			0.46			0.36			
## 8000			-119.44			-83.76			12.43			0.29			13			0.38			0.02			0.84			0.36			0.69			0.31			0.66			0.96			0.46			0.37			
## 9000			-99.78			-89.44			8.48			0.2			14			0.38			0.02			0.86			0.34			0.68			0.31			0.67			0.96			0.45			0.34			
## 10000			-96.6			-85.31			5.12			0.21			13			0.38			0.02			0.85			0.34			0.67			0.3			0.65			0.96			0.41			0.34			
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
## 0			1000			-381.7			-59.23			-51.05			15.96			0.35			8			0.79			0.01			1			1			1			0.71			0.74			1			1			1			
## 0			2000			-192.78			-17.53			-9.22			10.73			0.2			0			0.8			0.01			1			1			1			0.75			0.77			1			1			1			
## 0			3000			-151.5			-43.12			-33.88			22.04			0.57			5			0.78			0.01			1			1			1			0.75			0.78			1			1			1			
## 0			4000			-373.28			-48.32			-39.3			11.61			0.34			6			0.8			0.01			1			1			1			0.75			0.79			1			1			1			
## 0			5000			-405.73			-28.27			-21.03			19.45			0.29			2			0.79			0.01			1			1			1			0.74			0.79			1			1			1			
## 0			6000			-882.03			-28.36			-21.75			15.96			0.31			2			0.79			0.01			1			1			1			0.74			0.79			1			1			1			
## 0			7000			-445.02			-33.29			-24.47			11.58			0.19			3			0.8			0.01			1			1			1			0.74			0.78			1			1			1			
## 0			8000			-493.64			-27.81			-19.66			14.35			0.31			2			0.8			0.01			1			1			1			0.74			0.78			1			1			1			
## 0			9000			-154.41			-17.56			-10.7			26.01			0.63			0			0.8			0.01			1			1			1			0.74			0.79			1			1			1			
## 0			10000			-175.76			-22.8			-15.79			11.15			0.34			1			0.79			0.01			1			1			1			0.74			0.79			1			1			1			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.25			1000			-139.71			-68.14			-58.65			30.14			0.44			10			0.56			0.01			1			0.81			1			0.52			0.75			1			1			1			
## 0.25			2000			-145.05			-43.48			-36.58			29.06			0.49			5			0.59			0.01			1			0.9			1			0.55			0.75			0.98			1			1			
## 0.25			3000			-147.46			-38.06			-28.19			14.8			0.32			4			0.57			0.01			0.98			0.9			1			0.54			0.71			0.98			1			1			
## 0.25			4000			-157.66			-43.65			-36.01			19.36			0.61			5			0.58			0.01			0.98			0.9			1			0.55			0.71			0.99			1			1			
## 0.25			5000			-135.84			-43.22			-33.84			16.04			0.33			5			0.57			0.01			0.96			0.83			1			0.55			0.69			0.98			0.9			0.95			
## 0.25			6000			-143.42			-37.56			-27.3			20.16			0.42			4			0.58			0.01			0.97			0.82			1			0.54			0.68			0.98			0.86			0.91			
## 0.25			7000			-146.65			-38.94			-32.46			22.54			0.36			4			0.58			0.01			0.97			0.82			1			0.54			0.69			0.98			0.87			0.92			
## 0.25			8000			-141.83			-27.85			-20.4			26.79			0.48			2			0.57			0.01			0.98			0.81			1			0.55			0.68			0.99			0.88			0.93			
## 0.25			9000			-139.1			-39.23			-34.9			29.98			0.39			4			0.58			0.01			0.97			0.8			1			0.55			0.68			0.98			0.89			0.92			
## 0.25			10000			-150.63			-33.04			-24.53			16.21			0.49			3			0.58			0.01			0.97			0.81			1			0.54			0.68			0.98			0.9			0.92			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.5			1000			-142.13			-74.04			-67.92			46.83			0.64			11			0.47			0.01			0.89			0.63			1			0.55			0.75			1			0.85			1			
## 0.5			2000			-137.71			-64.83			-61.88			43.86			0.45			9			0.5			0.01			0.91			0.66			1			0.53			0.75			1			0.78			0.79			
## 0.5			3000			-146.42			-53.63			-44.15			24.67			0.46			7			0.52			0.01			0.92			0.66			1			0.53			0.77			1			0.79			0.83			
## 0.5			4000			-141.42			-40.85			-38.56			34.74			0.47			4			0.53			0.01			0.94			0.67			1			0.51			0.76			1			0.83			0.84			
## 0.5			5000			-142.3			-43.52			-34.25			18.94			0.4			5			0.52			0.01			0.94			0.69			1			0.51			0.75			1			0.85			0.85			
## 0.5			6000			-138.4			-42.84			-33.54			27.62			0.39			5			0.52			0.01			0.95			0.69			1			0.5			0.73			1			0.85			0.84			
## 0.5			7000			-140.55			-27.88			-22.34			34.95			0.44			2			0.52			0.01			0.95			0.69			1			0.49			0.73			0.99			0.87			0.86			
## 0.5			8000			-138.25			-37.9			-32.02			41.61			0.5			4			0.52			0.01			0.95			0.7			1			0.48			0.72			0.99			0.88			0.87			
## 0.5			9000			-141.85			-48.32			-41.2			28.62			0.38			6			0.52			0.01			0.96			0.7			1			0.48			0.72			0.99			0.89			0.88			
## 0.5			10000			-139.78			-42.76			-35.35			39.55			0.55			5			0.52			0.01			0.96			0.71			1			0.48			0.72			0.98			0.87			0.89			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 0.75			1000			-133.49			-90.14			-82			28.08			0.42			14			0.41			0.03			0.65			0.62			1			0.38			0.74			0.86			0.83			0.7			
## 0.75			2000			-135.08			-63.24			-54.7			31.52			0.44			9			0.43			0.02			0.82			0.59			1			0.4			0.74			0.91			0.75			0.69			
## 0.75			3000			-122.19			-68.92			-58.85			16.69			0.28			10			0.44			0.02			0.89			0.57			1			0.39			0.72			0.94			0.62			0.55			
## 0.75			4000			-127.64			-54.49			-45.74			17.96			0.35			7			0.43			0.02			0.88			0.51			0.95			0.39			0.67			0.96			0.47			0.43			
## 0.75			5000			-140.43			-75.48			-74.09			44.38			0.44			11			0.44			0.02			0.87			0.55			0.96			0.38			0.69			0.95			0.47			0.49			
## 0.75			6000			-133.82			-74.62			-67.37			24			0.39			11			0.44			0.02			0.89			0.57			0.95			0.4			0.71			0.94			0.53			0.56			
## 0.75			7000			-132.57			-64.52			-56.6			23.8			0.28			9			0.45			0.02			0.91			0.58			0.95			0.4			0.72			0.95			0.59			0.55			
## 0.75			8000			-128.41			-73.96			-63.51			20.4			0.3			11			0.45			0.02			0.91			0.57			0.95			0.4			0.72			0.96			0.56			0.53			
## 0.75			9000			-134.45			-69.73			-61.84			26.37			0.33			10			0.44			0.02			0.92			0.57			0.94			0.39			0.72			0.95			0.55			0.52			
## 0.75			10000			-124.47			-83.91			-71.9			22.58			0.32			13			0.44			0.02			0.92			0.59			0.93			0.39			0.73			0.96			0.57			0.55			
## Bk_k			gen			lnL			prior			ref			half.life			Vy			K			.alpha			birth.k			D0.slide			D1.slide			death.k			.sig2			.theta			U0.slide			U1.slide			U2.slide			
## 1			1000			-139.04			-64.62			-62.17			47.18			0.48			9			0.42			0.02			0.93			0.47			1			0.39			0.74			1			0.92			1			
## 1			2000			-124.38			-76.01			-70.75			17.39			0.29			11			0.44			0.02			0.88			0.48			1			0.33			0.79			1			0.82			1			
## 1			3000			-125.02			-86.15			-81.92			29.92			0.38			13			0.39			0.02			0.84			0.42			0.9			0.35			0.78			0.91			0.72			0.85			
## 1			4000			-126.61			-87.38			-83.81			24.82			0.27			13			0.39			0.02			0.8			0.41			0.9			0.35			0.78			0.92			0.66			0.73			
## 1			5000			-115.07			-113.47			-105.87			18.83			0.29			18			0.37			0.03			0.81			0.41			0.87			0.33			0.77			0.93			0.61			0.72			
## 1			6000			-119.65			-90			-79.35			19.14			0.34			14			0.37			0.03			0.81			0.41			0.86			0.33			0.75			0.93			0.56			0.6			
## 1			7000			-107.56			-120.15			-107.77			13.99			0.24			20			0.35			0.03			0.82			0.4			0.86			0.32			0.73			0.93			0.54			0.59			
## 1			8000			-107.72			-103.06			-96			12.16			0.23			16			0.35			0.03			0.82			0.41			0.85			0.32			0.72			0.93			0.54			0.53			
## 1			9000			-113.72			-94.66			-82.92			10.07			0.26			15			0.35			0.03			0.83			0.43			0.82			0.31			0.71			0.94			0.54			0.51			
## 1			10000			-93.12			-94.32			-82.22			7.26			0.19			15			0.35			0.03			0.84			0.41			0.81			0.31			0.69			0.94			0.53			0.48			
## Loading mcmc chains...
```

```r
ss
```

```
## Stepping stone estimation of marginal likelihood
## Marginal Likelihood:
## [1] -145.548
## A total of 5 power posteriors were run along the sequence: 0		0.25		0.5		0.75		1
## lnr_k -40.54 -35.86 -35.8 -33.34
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
plotBayoupars(converge_pars, tree, cex=0.5)
```

```
## Error in names(segs) <- tmp.o: 'names' attribute [453] must be the same length as the vector [451]
```

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
## 1000			-148.29			-5.87			106.31			1.24			2			0.39			1			0.34			0.32			1			
## 2000			-142.79			-5.85			64.77			0.68			2			0.4			0.97			0.34			0.38			1			
## 3000			-138.95			-5.77			41.61			0.48			2			0.41			0.96			0.34			0.4			1			
## 4000			-139.69			-5.8			35.79			0.5			2			0.42			0.96			0.34			0.4			1			
## 5000			-139.38			-5.76			40.23			0.43			2			0.43			0.97			0.34			0.41			1			
## 6000			-140.57			-5.82			41.15			0.5			2			0.43			0.97			0.33			0.41			1			
## 7000			-138.12			-5.76			32.89			0.43			2			0.44			0.97			0.34			0.42			1			
## 8000			-142.49			-5.81			53.29			0.71			2			0.43			0.97			0.34			0.42			1			
## 9000			-140.18			-5.79			29.85			0.49			2			0.43			0.97			0.34			0.42			1			
## 10000			-140.24			-5.76			28.5			0.36			2			0.43			0.97			0.34			0.41			0.99			
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
##                    Mean         SD     Naive SE Time-series SE
## lnL       -140.42679836 1.76082659 0.0665055026   0.1963593169
## prior       -5.78371158 0.04370744 0.0016508071   0.0042180309
## alpha        0.02038045 0.00448490 0.0001693923   0.0004409626
## sig2         0.01871892 0.00265805 0.0001003932   0.0002646474
## k            2.00000000 0.00000000 0.0000000000   0.0000000000
## ntheta       3.00000000 0.00000000 0.0000000000   0.0000000000
## root         3.72587465 0.17000835 0.0064211268   0.0151826013
## all theta    3.60023378 0.31942313           NA             NA
##           Effective Size
## lnL             80.41372
## prior          107.37203
## alpha          103.44320
## sig2           100.87678
## k                0.00000
## ntheta           0.00000
## root           125.38574
## all theta             NA
## 
## 
## Branches with posterior probabilities higher than 0.1:
##     pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
## 409  1            3.613110         0.01032140     1.479101
## 411  1            3.461717         0.01546314     2.148915
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
ss.fixed <- steppingstone(Bk=seq(0,1,length.out=5), chain.fixed, tree, dat, SE=0, startpar=startpar, prior=prior.fixed, ngen=10000, cores=5)
```

```
## Making power posterior function from provided mcmc chain...
```

```
## Running mcmc chains...
```

```
## Error in .steppingstone.mcmc(k = k, Bk = Bk, tree = tree, dat = dat, SE = SE, : unused argument (cores = 5)
```

```r
ss.fixed <- set.burnin(ss.fixed, 0.3)
```

```
## Error in set.burnin(ss.fixed, 0.3): object 'ss.fixed' not found
```

```r
ss.fixed$lnr
```

```
## Error in eval(expr, envir, enclos): object 'ss.fixed' not found
```

```r
plot(ss.fixed)
```

```
## Error in plot(ss.fixed): object 'ss.fixed' not found
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28-1.png) 

And calculate Bayes Factors in support of the Reversible-Jump model vs. the fixed model. 

```r
2*(ss$lnr-ss.fixed$lnr)
```

```
## Error in eval(expr, envir, enclos): object 'ss.fixed' not found
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
##                    Mean         SD     Naive SE Time-series SE
## lnL       -140.42679836 1.76082659 0.0665055026   0.1963593169
## prior       -5.78371158 0.04370744 0.0016508071   0.0042180309
## alpha        0.02038045 0.00448490 0.0001693923   0.0004409626
## sig2         0.01871892 0.00265805 0.0001003932   0.0002646474
## k            2.00000000 0.00000000 0.0000000000   0.0000000000
## ntheta       3.00000000 0.00000000 0.0000000000   0.0000000000
## root         3.72587465 0.17000835 0.0064211268   0.0151826013
## all theta    3.60023378 0.31942313           NA             NA
##           Effective Size
## lnL             80.41372
## prior          107.37203
## alpha          103.44320
## sig2           100.87678
## k                0.00000
## ntheta           0.00000
## root           125.38574
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
## [1] 0.1998618
## 
## $sig2
## [1] 0.06799417
## 
## $k
## [1] 13
## 
## $ntheta
## [1] 14
## 
## $theta
##  [1] 3.800425 4.596740 4.754152 4.175607 4.499691 3.229389 3.144065
##  [8] 2.569950 4.263980 2.999493 4.754563 3.346101 1.352504 4.255794
## 
## $sb
##  [1]  45 408 372 361 411 440 316 252 200 402   1  79 300
## 
## $loc
##  [1] 20.2589703  5.4853744  2.8849834  3.0989981 43.0596109  1.4234151
##  [7]  5.3358854  4.2820953  1.6291464 13.3524716  0.4098723  0.1897309
## [13]  8.6149056
## 
## $t2
##  [1]  2  3  4  5  6  7  8  9 10 11 12 13 14
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


