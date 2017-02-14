

devtools::install_github("joeguinness/aldo")

# a short vignette demonstrating how to use the functions
library(aldo)

# grid size for data locations
gsize <- 40
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
locs <- simulateGrid(nvec,jittersize=0)
plot(locs[,1],locs[,2])

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, sig2noise = 1)

# simulate some data and plot a map of it
# uses Cholesky method so don't try this for large n
X <- cbind( rep(1,n), locs[,1], locs[,2] )
beta <- c(1,6,0)
y <- X %*% beta + simulateData(locs,covparms,covfun)
image( matrix(y,nvec) )

# generate an ordering and plot the first n/8
ord <- orderMaxMinLocal(locs)
n0 <- round(n/8)
plot( locs[ord[1:n0],1],locs[ord[1:n0],2] )

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]
Xord <- X[ord,]

# find the ordered m nearest neighbors
m <- 30
NNarray <- findOrderedNNfast(locsord,m)

# automatically group the observations
NNlist <- groupNN(NNarray)
NNlist[1:4]  # list elements are subsets of rows of NNarray



# compute exact loglik, ungrouped, and grouped ordered composite logliks
# system.time(  ll0 <- mvnMargLik(covparms,covfun,y,locs) # only do this if n is small 
system.time(  ll1 <- orderedCompLik(covparms,covfun,yord,locsord,NNarray)      )
system.time(  ll2 <- orderedGroupCompLik(covparms,covfun,yord,locsord,NNlist)  )



# an attempt to write a wrapper function to do all of this stuff:
# ordering, finding neighbors, maximizing parameters
# only covariance function implemented is isotropic matern 
# interesting thing: block independent likelihood is ridiculously
# good at finding parameter estimates that nearly maximize
# vecchia's likelihood approximation
# This function uses block independent likelihood to quickly get starting values
# for an optimization of Vecchia's approximation
system.time( result <- fitmodel(y,X,locs,maternIsotropic,numneighbors=30,fixedparameters=c(1,NA,NA,1))  )

result









