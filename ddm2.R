# Do the wiener fit:
# Some code copied from Kruske's textbook, other code taken from the Wabersich's examples
setwd("~/Documents/stats")

require(rjags)
load.module("wiener")

#--------------------------------- Read in data:
ADM_data <- read.csv("dataADM_PLSresults.csv")
rt <- ADM_data["rt"]
Ndata <- nrow(rt)
x <- ADM_data[,6:10]
nPredictors <- ncol(x)

# Specify data, as a list.
dataList <- list(
  rt <- rt[,] ,
  Ndata <- as.numeric(Ndata) ,
  x <- x[,] ,
  nPredictors <- as.numeric(nPredictors) ,
  Rc <- structure(c(1,0,0,1), .Dim=c(2,2))
) # Error reading in the data - nPredictors doesn't exist?

# -------------------------------- Set up the model

modelstring = "
model {

  for(i in 1 : Ndata){
    rt[i] ~ dwiener(alphadelta[i,1], tau, beta, alphadelta[i,2])

    alphadelta[i,] <- exp(k[i,])
    k[i,] ~ dmnorm( m[i,], C)

    mu[i,1] <- b01 + inprod(b1[], x[i,] )
    mu[i,2] <- b02 + inprod(b2[], x[i,])
  }
  
  C ~ dwish(Rc, 2)
  tau ~ dunif(0 , 1)
  beta <- 0.5
  
  b01 ~ dnorm(0, 0.000001)
  b02 ~ dnorm(0, 0.000001)

  for(j in 1 : nPredictors){
    b1[j] ~ dnorm(0, 0.000001)
    b2[j] ~ dnorm(0, 0.000001)
  }

}
"
# close quote for modelstring - copied from Kruske
writeLines(modelstring,con="model.txt")

# -------------------------------- Initialize chains

#initsList = list(
#  tau = 0.5 ,
#  alphadelta = 1
#)

# -------------------------------- Run the model

parameters = c("b0","b","C","tau","alpahdelta")
burnInSteps = 500            # Number of steps to "burn-in" the samplers.
nsteps = 5000

jagsModel = jags.model("model.txt", data = dataList)#, inits = initsList)

cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nsteps , thin=1 )

# -------------------------------- Examine results

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract chain values:
pad = mcmcChain[, "alphadelta" ]
ptau = mcmcChain[, "tau" ]
pC = mcmcChain[,"C"]


# Stupidly plot them:
hist(pad[,1])
mean(pad[,1])
hist(pad[,2])
mean(pad[,2])
hist(pC)
mean(pC)

# --------------------------------


# --------------------------------








