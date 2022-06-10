
cat("\nInstall/load packages")
inst.pkgs <- installed.packages()

## required packages:
l.pkgs <- c("expm",
            "Matrix",
            "parallel",
            "gaussquad",
            "splines",
            "scales",
            "mvtnorm",
            "tmvtnorm",
            "MASS",
            "igraph",
            "stringr",
            "devtools")
## check if packages are installed
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})

lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})

if(!("Karen" %in% rownames(inst.pkgs))){
  install_github("delcore-luca/Karen",
                 ref = "master")
}
library(Karen)

rm(list = ls())

source("./source/Karen-sim-MS-fun.R")

##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- substr(currTime, start = 1, stop = 13)

#############
## folders ##
#############

resPath <- "./results"
currResPath <- paste(resPath, "/", currTime, sep = "")
currLogPath <- paste(currResPath, "/processors/", sep = "")

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(currResPath)), dir.create(file.path(currResPath)), FALSE)
ifelse(!dir.exists(file.path(currLogPath)), dir.create(file.path(currLogPath)), FALSE)

## load candidate CDNs:
model.lst <- get.modelList()
nModels <- length(model.lst)
trueCandMod <- commandArgs(trailingOnly = TRUE)
trueMod <- as.numeric(as.vector(unlist(strsplit(trueCandMod, split = "-", fixed = TRUE)))[1])
candMod <- as.numeric(as.vector(unlist(strsplit(trueCandMod, split = "-", fixed = TRUE)))[2])

## CDN parameters of true model:
rcts <- as.vector(unlist(model.lst[[trueMod]]["rct.lst"]))
cnstr <- as.vector(unlist(model.lst[[trueMod]]["constr.lst"]))
latsts <- as.vector(unlist(model.lst[[trueMod]]["latSts.lst"]))
ctps <- unique(setdiff(c(sapply(rcts, function(r){
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

## true dynamic parameters:
if(trueMod == 1){
  th.true <- c(1.3/2, 1.8/2, 1.85/2, 1.95/2, 1.1/2, 2.5 + 1, 2.1 + 1, 3 + 1, 2.7 + 1, 3.1 + 1, .5/2, .45/2, .55/2)
}else{
  th.true <- c(1.3/2, 1.8/2, 1.85/2, 1.95/2, 1.1/2, 1.2/2, 1.4/2, 2.5 + 1, 2.1 + 1, 3 + 1, 2.7 + 1, 3.1 + 1, .5/2, .45/2, .55/2, .48/2)
}
names(th.true) <- tail(rcts, -length(cnstr))

## simulation parameters:
s2.true <- 1e-8
r0.true <- .1
r1.true <- .5
phi.true <- c(th.true, r0.true, r1.true)
names(phi.true) <- c(names(th.true), "r0", "r1")
S <- 1000
nCL <- 3

## initial condition X0 for the true model:
X0 <- rep(0, length(ctps))
names(X0) <- ctps
X0["HSC"] <- 100

nSim <- 100 ## number of independent simulations
ntps <- 15 ## number of time-points
f_NA <- 0 ## fraction of latent data

## initialize list of simulation results:
resAllSim <- vector("list", length = nSim)
AIC_allsim <- rep(NA,nSim)

## simulate all the independent trajectories:
XY_allsim <- lapply(1:nSim, function(sim){
  cat(paste("simulation n. ", sim, "\n", sep = ""))
  XY <- get.sim.trajectories(rct.lst = rcts,
                             constr.lst = cnstr,
                             latSts.lst = latsts,
                             ct.lst = ctps,
                             th = th.true,
                             S = S,
                             nCL = nCL,
                             X0 = X0,
                             s2 = s2.true,
                             r0 = r0.true,
                             r1 = r1.true,
                             f = f_NA,
                             ntps = ntps,
                             trunc = FALSE)
  return(XY)
})

## CDN parameters for the candidate model:
rcts <- as.vector(unlist(model.lst[[candMod]]["rct.lst"]))
cnstr <- as.vector(unlist(model.lst[[candMod]]["constr.lst"]))
latsts <- as.vector(unlist(model.lst[[candMod]]["latSts.lst"]))
ctps <- unique(setdiff(c(sapply(rcts, function(r){
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

for (sim in 1:nSim) {
  cat(paste("simulation n. ", sim, "\t\tTrue model: ", trueMod, "\t Candidate model: ", candMod, "\n", sep = ""))

  ## current simulation:
  XY.sim <- XY_allsim[[sim]]
  Y.sim <- XY.sim$Y
  X.sim <- XY.sim$X

  ## start a cluster:
  nProc <- 1
  cat(paste("\tLoading CPU cluster...\n", sep = ""))
  cat(paste("Cluster type: ", "PSOCK\n", sep = ""))
  cpu <- Sys.getenv("SLURM_CPUS_ON_NODE", nProc)
  hosts <- rep("localhost",cpu)
  cl <- makeCluster(hosts, type = "PSOCK")
  rm(nProc)

  ## initial condition X0 for the candidate model:
  X0 <- rep(0, length(ctps))
  names(X0) <- ctps
  X0["HSC"] <- 100

  ## mean vector and covariance matrix of X0:
  m_0 <- replicate(nCL, X0, simplify = "array")
  colnames(m_0) <- 1:nCL
  P_0 <- Diagonal(length(ctps) * nCL, 1e-5)
  rownames(P_0) <- colnames(P_0) <- rep(1:nCL, each = length(ctps))

  ## fit Karen on data:
  res.fit <- get.fit(rct.lst = rcts,
                     constr.lst = cnstr,
                     latSts.lst = latsts,
                     ct.lst = ctps,
                     Y = Y.sim[,setdiff(ctps, latsts),],
                     m0 = m_0,
                     P0 = P_0,
                     cl = cl,
                     list(nLQR = 3,
                          lmm = 25,
                          pgtol = 0,
                          relErrfct = 1e-5,
                          tol = 1e-9,
                          maxit = 1000,
                          maxitEM = 10,
                          trace = 1,
                          FORCEP = FALSE))

  stopCluster(cl) ## stop the cluster
  phi.curr <- res.fit$fit$par ## MLE parameters

  ## smoothing moments:
  pdf(file = paste(currResPath, "/sim", sim, "_T", ntps, "_trueMod", trueMod, "_candMod", candMod, "_sMoments.pdf", sep = ""), width = 14, height = 5)
  par(mar = c(5,5,2,2), mfrow = c(ceiling(nCL/ceiling(sqrt(nCL))),ceiling(sqrt(nCL))))
  get.sMoments(res.fit = res.fit, X = XY.sim$X)
  dev.off()

  res <- list()
  res$phi.opt <- res.fit$fit$par
  res$bwd.res <- res.fit$bwd.res
  res$Y <- XY.sim$Y
  res$X <- XY.sim$X
  resAllSim[[sim]] <- res
  AIC_allsim[sim] <- res.fit$AIC
}

save.image(paste(currResPath, "/output_trueMod", trueMod, "_candMod", candMod, ".Rdata", sep = ""))

