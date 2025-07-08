library(nimble)
library(coda)
source("sim.RT.R")
source("NimbleModelRT Poisson DA2 Marginal.R")
source("NimbleFunctionsRT Poisson DA2 Marginal.R")
source("init.RT.R")
source("sSampler Poisson Marginal.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
N <- 78
#detection parameters
lam0 <- 0.5
sigma <- 0.5

K <- 10 #number of occasions
buff <- 3 #state space buffer
X <- expand.grid(3:11,3:11) #make a trapping array
#ID thinning parameter
theta.thin <- 0.25

data <- sim.RT(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,theta.thin=theta.thin,obstype="poisson")

#What is the observed data?
str(data$y.ID) #the observed ID detections
str(data$y.noID) #trap by occasion counts of unidentified detections

#To fit models faster, we will sum data over occasions
data$y.ID <- apply(data$y.ID,c(1,2),sum)
data$y.noID <- rowSums(data$y.noID)

##Fit model in Nimble##
#data augmentation level
M <- 175
X <- data$X
J <- nrow(X) #number of detectors

inits <- list(lam0=1,sigma=1) #ballpark inits to build data

#augment, initialize data, check starting obsmod logProb
nimbuild <- init.RT(data,inits,M=M)
plot(NA,xlim=data$xlim,ylim=data$ylim)
points(X,pch=4)
points(nimbuild$s,pch=16) #initialized activity centers
for(i in 1:data$n.cap){
  trapcaps <- which(data$y.ID[i,]>0)
  traps <-  rbind(X[trapcaps,])
  s <- nimbuild$s[i,]
  points(s[1],s[2],col="goldenrod",pch=16)
  if(nrow(traps)>0){
    for(j in 1:nrow(traps)){
      lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
    }
  }
}

#inits for nimble
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #z and N inits must be consistent
                 lambda.N=sum(nimbuild$z), #converges faster if you set N and lambda.N at similar values
                 s=nimbuild$s,sigma=inits$sigma,lam0=inits$lam0,theta.thin=0.5)

#constants for Nimble
constants <- list(M=M,J=J,xlim=data$xlim,ylim=data$ylim)

# Supply data to Nimble
Nimdata <- list(y.ID=nimbuild$y.ID,X=as.matrix(X),K1D=data$K1D,y.noID=data$y.noID)

# set parameters to monitor
parameters <- c('lam0','sigma','theta.thin','N','lambda.N')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c('lam0','sigma','theta.thin','lambda.N')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,nodes=config.nodes)

#Add N/z update
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal.
#nodes used for update
y.ID.nodes <- Rmodel$expandNodeNames(paste("y.ID[1:",M,",1:",J,"]"))
y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
bigLam.nodes <- Rmodel$expandNodeNames("bigLam") #only need this in calcNodes
lam.noID.nodes <- Rmodel$expandNodeNames("lam.noID")
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.noID.nodes,y.ID.nodes)
inds.detected <- which(rowSums(nimbuild$y.ID)>0)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,inds.detected=inds.detected,
                                                 lam.nodes=lam.nodes,lam.noID.nodes=lam.noID.nodes,
                                                 y.ID.nodes=y.ID.nodes,y.noID.nodes=y.noID.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)

#add activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,J=J,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#add block sampler for lam0 and sigma if posteriors correlated (usually are)
conf$addSampler(target = c("lam0","sigma"),
                  type = 'RW_block',
                  control = list(adaptive=TRUE),
                  silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[-c(1:burnin),]))

#posterior correlation
tmp <- cor(mcmc(mvSamples[-c(1:burnin),]))
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)
