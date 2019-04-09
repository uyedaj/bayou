This tutorial guides the user through analyses using bayou versions 2.0 and above. There have been some changes from bayou 1.0, particularly in how the MCMC analyses are built and deployed. 

Load packages
```{r}
install.packages("bayou")
library(bayou)

```

# Simulating a multi-optima OU process on a phylogeny

Let's begin by simulating a multi-optimum Ornstein-Uhlenbeck process on a phylogeny so that we can get a feel for how these models work to model adaptation on phylognies. First let's simulate a tree and rescale it to 100 million years. The second step is optional, but it will help us make sure all our trees and parameters are on a common scale that will probably be typical of the trees you'll analyze. We will also reorder the tree into "postorder" format. **bayou** will automatically reorder the tree to postorder format, but it helps to begin by reordering your tree so that the branch numbers used by **bayou** can be easily matched to the tree. 

```{r}

  tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 50, seed = 1)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))*100
  tree <- reorder(tree, "postorder")

  plot(tree, cex = 0.5)
  
```

Now let's simulate some data and use **bayou** to estimate adaptive shifts on phylogenies. Let's first define a set of parameter values to use to simulate our data. 

```{r}

set.seed(1)
truepars <- list(alpha=0.1, sig2=0.1, 
                    k=4, ntheta=5, theta=rnorm(5, 0, 4))
trueshiftlocations <- list(sb = c(94, 71, 45, 12), loc = c(6.9, 2.2, 0.8, 3.9)) # I encourage you to use identifyBranches instead.
truepars <- c(truepars, trueshiftlocations)
truepars$t2 <- 2:5

plotBayoupars(truepars, tree, cex = 0.5)

```

Now using the function *dataSim*, we can simulate trait data. 

```{r}

dat <- dataSim(truepars, tree, model="OU")$dat

```

To add realism, let's add some measurement error to the data. This is a good reminder to *always try to use measurement error in your analyses*. OU models especially are affected by measurement error. This is because OU models have the effect of "erasing" evolutionary history with increasing *alpha*. If you don't account for measurement error, then that measurement error will be transferred to the evolutionary process model. You can make a Brownian Motion model look very OU like if there is a lot of measurement error.

```{r}

MEvar <- 0.1
dat <- dat + rnorm(length(dat), 0, sqrt(MEvar))

```


We can now define our prior for our model. The prior function is going to take our parameters and output the *prior probability* of our parameter values. It represents our initial degree of belief in what values the parameters will take. 

```{r}

priorOU <- make.prior(tree, 
                      dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                 dk="cdpois", dtheta="dnorm"),
                      param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                 dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(dat), sd=1.5*sd(dat)))
)

```


*make.prior* tries to be reasonable and not make you type everything out, but **do not be comfortable with defaults**. One trick to make sure your prior functions are reasonable is to simulate a bunch of values and take the quantiles of the distribution. We are using a half-Cauchy distribution for *alpha* and *sig2*, which is a good weakly informative prior for scale parameters. 

To run our MCMC, we have to initiate the MCMC chain with some starting values. It's good to run multiple chains from multiple different starting values. Let's simulate some values from the prior distribution and make sure our prior functions works.

```{r}

startpars <- priorSim(priorOU, tree, plot=TRUE)$pars[[1]]
priorOU(startpars)

```


We're now going to take what we have and put it into the function *bayou.makeMCMC*. This function does not immediately initiate the MCMC, but it makes an object that we can use to manage our MCMC analysis. When *bayou* runs an MCMC, it writes the output to a set of files that need to put somewhere. This ensures that the memory doesn't get full of increasingly long chains. Here, I've specified the new.dir=TRUE, which means it will save the output into your R temporary directory. You can specify a particular directory if you wish. 


```{r}
set.seed(1)
mcmcOU <- bayou.makeMCMC(tree, dat, SE=MEvar, prior=priorOU, 
                         new.dir=TRUE, outname="modelOU_r001", plot.freq=NULL) # Set up the MCMC
mcmcOU$run(10000) # Run the MCMC

```

The full MCMC results are written to a set of files. We can load them back in to R as follows.

```{r}

chainOU <- mcmcOU$load()

```

Let's take a look at the results. We can set a "burnin" parameter that tells the package **coda** to discard the first bit of the chain.

```{r}

chainOU <- set.burnin(chainOU, 0.3)
summary(chainOU)
plot(chainOU, auto.layout=FALSE)

```

Our traces will probably look bad, 10,000 generations isn't long enough to obtain convergence. Also, note the small effective sample sizes in our summary (the NA's for the *all theta* row are expected, this is because these aren't a single parameter, but a variable number of optima that are coming in and out of existence throughout the chain). 

Let's visualize what we have so far. First, we will plot the truth, then 3 alternative ways of visualizing our chain.

```{r}

par(mfrow=c(2,2))
plotBayoupars(truepars, tree, main = "True parameters")
plotSimmap.mcmc(chainOU, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(tree, chainOU, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(tree, dat, burnin = 0.3, chainOU, pp.cutoff = 0.3)

```

Even though we haven't gotten convergence yet, we're probably picking up the major shifts pretty well. 

# Alternative parameterizations

Two alternative parameterizations of the OU model are built into **bayou**. First, is a parameterization where priors can be placed directly on phylogenetic half-life (*halflife*) and stationary variance (*Vy*), rather than *alpha* and *sig2*. For example, let's say we want to have a mildly informative prior on the phylogenetic half-life--say a log-normal distribution:

```{r}

par.halflife <- list(meanlog=2.5, sdlog=2.5)
#Draw a bunch of samples from this distribution:
samp <- rlnorm(10000, par.halflife$meanlog, par.halflife$sdlog)
hist(log(samp,10), breaks=100, main="Prior density of halflife")
abline(v=log(c(1,max(branching.times(tree))),10), col="red", lwd=2, lty=2)

```

Notice that there is about equal density of prior probability on the half-life being greater than tree height (rightmost red line) as there is below 1 million years (leftmost red line). 

Now a prior on the stationary variance:
```{r}

par.Vy <- list(meanlog=log(0.1), sdlog=0.25)
hist(rlnorm(10000, par.Vy$meanlog, par.Vy$sdlog), main="Prior density of Vy")

```

Let's make the prior, MCMC object and run the chain. Note that we have specified that *model = "OUrepar"* in *make.prior* and *bayou.makeMCMC*, which means we are using the *halflife*, *Vy* parameterization instead of *sig2* and *alpha*.

```{r}

priorBB <- make.prior(tree, 
                      dists=list(dhalflife="dlnorm", dVy="dlnorm", 
                                 dk="cdpois", dsb="dsb", dtheta="dnorm"),
                      param=list(dhalflife=par.halflife,
                                 dVy=par.Vy,
                                 dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                      model="OUrepar"
)


set.seed(1)
mcmcBB <- bayou.makeMCMC(tree, dat, SE=MEvar, model="OUrepar", prior=priorBB, new.dir=TRUE, outname="modelBB_r001", plot.freq=NULL)
mcmcBB$run(10000)

chainBB <- mcmcBB$load()
chainBB <- set.burnin(chainBB, 0.3)
summary(chainBB)
plot(chainBB)

```


```{r}

par(mfrow=c(2,2))
plotBayoupars(truepars, tree, main = "True parameters")
plotSimmap.mcmc(chainBB, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(tree, chainBB, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(tree, dat, burnin = 0.3, chainBB, pp.cutoff = 0.3)

```

Likely, we will have more shifts because we made the prior on *Vy* so narrow. Let's compare the posteriors from the two models.

```{r}

quantile(chainOU$sig2/(2*chainOU$alpha), quantiles)
quantile(chainBB$Vy, quantiles)

```

# Quantitative Genetics Model
We can also fit a model that follows the Quantitative Genetics parameterization. You can check the specified prior distributions on your own, but these are informative priors that specify moderate to high heritability, realistic phenotypic variances, large uncertainty regarding the strength of selection and reasonable effective population sizes for entire species. 

```{r}

par.h2 <- list(shape1=10, shape2=10)
par.P <- list(meanlog=log(0.12), sdlog=0.2)
par.w2 <- list(meanlog=log(100), sdlog=2.5)
par.Ne <- list(meanlog=log(500000), sdlog=2.5)

```

We should rescale the branch lengths to correspond roughly to generation time. However, here we will assume that for
most of the history of birds and mammals, the generation time has been around 2 year/gen.

```{r}

QGtree <- tree
QGtree$edge.length <- QGtree$edge.length/2

```


```{r}

priorQG <- make.prior(QGtree, plot.prior=FALSE,
                      dists=list(dh2="dbeta", dP="dlnorm",
                                 dw2="dlnorm", dNe="dlnorm",
                                 dk="cdpois", dtheta="dnorm"),
                      param=list(dh2=par.h2,
                                 dP=par.P,
                                 dw2=par.w2,
                                 dNe=par.Ne,
                                 dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                      model="QG"
)

```

Note that this model has difficulty fitting if the starting point is a poor fit. So rather than drawing from the prior distribution, we will start with shifts chosen by previous analyses:

```{r}

set.seed(1)
mcmcQG <- bayou.makeMCMC(QGtree, dat, SE=MEvar, model="QG", startpar=NULL, prior=priorQG, new.dir=TRUE, outname="modelQG_r001", plot.freq=NULL)
mcmcQG$run(10000)

```


```{r}

chainQG <- mcmcQG$load()
chainQG <- set.burnin(chainQG, 0.3)
summary(chainQG)
plot(chainQG, auto.layout=FALSE)

```

```{r}

par(mfrow=c(2,2))
plotBayoupars(truepars, tree, main = "True parameters")
plotSimmap.mcmc(chainQG, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(tree, chainQG, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(tree, dat, burnin = 0.3, chainQG, pp.cutoff = 0.3)

```

If you kept the seed the same, you should see that there are many more shifts recovered. This is because the QG model predicts such small stationary variances that shifts must occur costantly to explain the variation among species. 

# Model Comparison
Alternative parameterizations, shift locations, and priors can be compared using Bayes Factors. This requires estimation of the marginal likelihood, which can be difficult. **bayou** uses stepping-stone sampling to estimate the marginal likelihoods. To estimate marginal likelihoods, using the '$steppingstone' function in the mcmc object. For this exercise, we will do a much shorter run than is recommended. If you have multiple cores available on your machine, you can make use of these to run the stepping stone analysis in parallel and conduct the analysis much faster. 

While I have the complete code to do all 3 runs here, I suggest you partner with your neighbor and divide up the computational burden among you. 

```{r}
library(doParallel)
registerDoParallel(cores=2) #Using more cores will make it go faster if you have them.
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
ssOU <- mcmcOU$steppingstone(10000, chainOU, Bk, burnin=0.3, plot=FALSE)
ssBB <- mcmcBB$steppingstone(10000, chainBB, Bk, burnin=0.3, plot=FALSE)
ssQG <- mcmcQG$steppingstone(10000, chainQG, Bk, burnin=0.3, plot=FALSE)

mlnL <- c("OU"=ssOU$lnr, "BB"=ssBB$lnr, "QG"=ssQG$lnr)
mlnL

```

If you get a couple errors it's probably OK, the algorithm takes the posterior and tries to fit various distributions to the parameters, and if it fails to optimize them it will throw an error or two. However, as long as one of them fits OK it will run. Again, we have not run these for long enough or for enough steps (we prefer more like 50!), but you get the idea for how you would proceed. Obviously, having more cores makes this go a LOT faster, and this is a computationally intensive procedure!

```{r}

plot(ssOU)
plot(ssBB)
plot(ssQG)

```

# Fixed models
While using the reversible-jump MCMC of bayou is useful for exploratory analyses, it is likely that you will also have specific hypotheses regarding adaptive regimes. Like other approaches (OUwie, ouch, etc.) you can implement fixed models in bayou. We will set up two alternative hypotheses. First, we will set a prior with the shift locations fixed to be the true shift locations. Then we will specify an alternative prior with different shift locations. Finally, we will compare the two models using marginal likelihoods estimated using stepping stone sampling. 

```{r}

trueFixedPrior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         dk="fixed", dsb="fixed", 
                                         dtheta="dnorm"),
                                   param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                        dk="fixed", dsb="fixed", 
                                        dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                                   fixed =list(k = truepars$k, sb = truepars$sb)
                             )

```

Choose an alternative arrangement of shift locations. For fun, let's make this one have one extra shift to a convergent regime.

```{r}

altlocations <- list(sb = c(89, 70, 85, 47, 50), loc = c(1.5, 6.6, 5.2, 3.7, 11.1)) # You can also use identifyBranches here
altpars <- truepars
altpars$k <- 5
altpars$sb <- altlocations$sb
altpars$loc <- altlocations$loc
altpars$t2 <- c(2, 3, 4, 5, 3) # Shifts on branches 89 and 50 both lead to regime #3


alternativeFixedPrior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                              dk="fixed", dsb="fixed", 
                                              dtheta="dnorm"),
                             param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                        dk="fixed", dsb="fixed", 
                                        dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                             fixed=list(k = 5, ntheta = 5, sb = altpars$sb, loc = altpars$loc, t2 = altpars$t2)
)

par(mfrow=c(1,2))
plotBayoupars(truepars, tree, main="True Pars")
plotBayoupars(altpars, tree, main="Alternative Hypothesis")

```

Run both and load back into bayou. 

```{r}

mcmcFixed1 <- bayou.makeMCMC(tree, dat, SE=MEvar, prior=trueFixedPrior, new.dir=TRUE, outname="modelTrueFixed_r001", plot.freq=NULL)
mcmcFixed1$run(10000)

mcmcFixed2 <- bayou.makeMCMC(tree, dat, SE=MEvar, prior=alternativeFixedPrior, new.dir=TRUE, outname="modelAltFixed_r001", plot.freq=NULL)
mcmcFixed2$run(10000)

chainFixed1 <- mcmcFixed1$load()
chainFixed2 <- mcmcFixed2$load()

```

Again, we can estimate marginal likelihoods. Again, I suggest you divide up the tasks with your neighbor.

```{r}

## Run the stepping stone estimation of marginal likelihoods.
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
ssFixed1 <- mcmcFixed1$steppingstone(10000, chainFixed1, Bk)
ssFixed2 <- mcmcFixed2$steppingstone(10000, chainFixed2, Bk)

ssFixed1
ssFixed2

```


Most likely, we should see the true model with a much higher marginal likelihood than the alternative parameterization (even though again, these marginal likelihood estimates are likely bad).

***

# Customized/Allometric models
What if there is a known (or unknown) relationship between the trait of interest and another predictor variable? For example, we may be interested in a relationship between trait known to vary with body size, but consider the possibility that the relationship with body size itself varies over macroevolutionary time. Here, instead of having a single optimum that changes upon a regime shift, it is possible to have both the slope and intercept of the relationship change at once. bayou v2.0 allows you to include these additional predictors and test for shifts in the scaling between a trait and its predictors. 

Let's simulate a dataset where the slope and intercept shift at different points in the tree. We're going to use the same shift 
locations as before, but add in a covariate with body size that changes in different parts of the tree. We also need to simulate the predictor data, in this case, let's use Brownian Motion.

```{r}

set.seed(1)
tree <- sim.bdtree(b=1, d=0, stop="taxa", n=50, seed=1)
tree$edge.length <- tree$edge.length/max(branching.times(tree))*100
tree <- reorder(tree, "postorder")
truepars <- list(alpha = 0.5, sig2 = 0.05,
                 k = 3, ntheta = 4, 
                 beta_lnMass = c(0.75, 0.6, 0.9, 0.67), 
                 theta = c(-1, 1.25, 0.5, 0),
                 sb = c(94, 71, 50),
                 loc = c(0, 0, 0),
                 t2 = 2:4)

pred <- cbind("lnMass" = sim.char(tree, 0.2, model="BM", root=3)[,,1])
phytools::phenogram(tree, setNames(pred[,1], tree$tip.label), spread.labels=FALSE, main="Predictor: Body Size (lnMass)")

dat <- dataSim(truepars, tree, model="OU")$dat + truepars$beta_lnMass[bayou:::.tipregime(truepars, tree)] * pred[,1]

```

But our old visualization doesn't give the whole picture, because the trait covaries with body size:

```{r}

par(mfrow=c(1,2))
## Plot the regime locations
plotRegimes(pars2simmap(truepars,tree)$tr,  col=pars2simmap(truepars,tree)$col)
## Plot the allometry
plot(pred[,1], dat, pch=21, bg=bayou:::.tipregime(truepars, tree), xlab="lnMass", "ylab"="Trait")
## Add the regression lines
dum <- lapply(1:truepars$ntheta, function(x) abline(truepars$theta[x], truepars$beta_lnMass[x],  lty=2, col=x))

```

We are going to test 3 models in this analysis: Global intercepts & slopes (11), Separate intercepts & global slope (N1), and separate intercepts & slopes (NN). However, we're going to have to build these models to run them (they aren't built into **bayou**). As a convention, we're going to name our regression coefficients (other than the familiar intercept, *theta*) "beta_" followed by the predictor name (e.g. *beta_lnMass*). Here we imagine we have some fairly informative prior belief about what the allometry with body mass should be (normal distribution around 0.7).


```{r}

prior.11 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                  dsb="fixed", dk="fixed", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_lnMass=list(mean=0.7, sd=0.15),
                                  dtheta=list(mean=0, sd=1)),
                       fixed=list(k=0, sb=numeric(0))
)

prior.N1 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dbeta_lnMass=list(mean=0.7, sd=0.15),
                                  dk=list(lambda=10, kmax=50),
                                  dtheta=list(mean=0, sd=1))
)


prior.NN <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_lnMass=list(mean=0.7, sd=0.15),
                                  dk=list(lambda=10, kmax=50), 
                                  dtheta=list(mean=0, sd=1))
)

```

Manually set tuning parameters, and make the models. There is a bit of art to tuning the parameters, which may require making multiple runs and trying to get the acceptance ratios in the right region (0.2-0.4). But these should work well for these models and data. If the acceptance ratio for a certain parameter is too high, increase the tuning parameter for that variable. If the acceptance ratio is too low, decrease it. The scale of the regression coefficient, for example, should give you some idea of what these parameters should be. 

```{r}

D11 = list(alpha=2, sig2=2, beta_lnMass=0.1, k=1, theta=0.5, slide=1)
DN1 = list(alpha=2, sig2=2, beta_lnMass=0.1, k=1, theta=2, slide=1)
DNN = list(alpha=2, sig2=2, beta_lnMass=0.3, k=c(1,1), theta=2, slide=1)

```

Now we use the function *makeBayouModel* to create a **bayou** model object that specifies all the components **bayou** needs to drop a new model into the analysis. Note that if you are interested in developing in **bayou** these are intended to be easy to make for customized models and there is a lot more possible than what is shown here. Note that each model only differs in the number of reversible-jump parameters (0, 1, and 2), the prior and the tuning parameters. By default, the starting regime map is plotted.

```{r}

set.seed(1)
model.11 <- makeBayouModel(dat ~ lnMass, rjpars = c(), 
                           tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.11, D=D11)
model.N1 <- makeBayouModel(dat ~ lnMass, rjpars = c("theta"),  
                           tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.N1, D=DN1)
model.NN <- makeBayouModel(dat ~ lnMass, rjpars = c("theta", "lnMass"),  
                           tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.NN, D=DNN)

```

We can now drop these model object into the analysis, along with the generated starting values (replacing the out of the box options of "OU", "OUrepar" and "QG"). 

```{r}

## Make MCMC objects:
mcmc.11 <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.11$model, prior=prior.11, startpar=model.11$startpar, new.dir=TRUE, outname="model11_r001", plot.freq=NULL)
mcmc.N1 <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, new.dir=TRUE, outname="modelN1_r001", plot.freq=NULL)
mcmc.NN <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.NN$model, prior=prior.NN, startpar=model.NN$startpar, new.dir=TRUE, outname="modelNN_r001", plot.freq=NULL)

```

Run the models and load them in.

```{r}

set.seed(1)
mcmc.11$run(10000)
mcmc.N1$run(10000)
mcmc.NN$run(10000)

chain.11 <- set.burnin(mcmc.11$load(), 0.3)
chain.N1 <- set.burnin(mcmc.N1$load(), 0.3)
chain.NN <- set.burnin(mcmc.NN$load(), 0.3)

```

A particularly useful way to plot these is to use the *shiftSummaries* and *plotShiftSummaries* functions. Like other plotting functions, we define a posterior probability cutoff and only plot those shifts (*pp.cutoff*). Note that the global allometry (*11*), has no shifts and is not plotted here. 

```{r}

shiftsumsN1 <- shiftSummaries(chain.N1, mcmc.N1, pp.cutoff=0.5, burnin=0.3)
shiftsumsNN <- shiftSummaries(chain.NN, mcmc.NN, pp.cutoff=0.5, burnin=0.3)

plotShiftSummaries(shiftsumsN1, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsNN, lwd=2, single.plot=TRUE, label.pts=FALSE)

```

As before, we can compare different models by estimating marginal likelihoods. Divide and conquer. 

```{r}

registerDoParallel(cores=2)
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
set.seed(1)
ss.11 <- mcmc.11$steppingstone(10000, chain.11, Bk, burnin=0.3, plot=FALSE)
ss.N1 <- mcmc.N1$steppingstone(10000, chain.N1, Bk, burnin=0.3, plot=FALSE)
ss.NN <- mcmc.NN$steppingstone(10000, chain.NN, Bk, burnin=0.3, plot=FALSE)

mlnL <- c("11"=ss.11$lnr, "N1"=ss.N1$lnr, "NN"=ss.NN$lnr)
mlnL

```

That's it for now!

