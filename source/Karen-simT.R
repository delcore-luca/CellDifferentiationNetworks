
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

## CDN parameters:
rcts <- c("HSC->P1",
             "HSC->P2",
             "P1->T",
             "P1->B",
             "P1->NK",
             "P2->G",
             "P2->M",
             "T->0",
             "B->0",
             "NK->0",
             "G->0",
             "M->0"
             ,"HSC->1"
             ,"P1->1"
             ,"P2->1"
)

cnstr <- c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->T\\'\\] + theta\\[\\'P1->B\\'\\] + theta\\[\\'P1->NK\\'\\])",
                "theta\\[\\'HSC->P2\\'\\]=(theta\\[\\'P2->G\\'\\] + theta\\[\\'P2->M\\'\\])")
latsts <- c("HSC", "P1", "P2")

ctps <- unique(setdiff(c(sapply(rcts, function(r){
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

## true dynamic parameters:
th.true <- c(0.65, 0.9, 0.925, 0.975, 0.55, 3.5, 3.1, 4, 3.7, 4.1, 0.25, 0.225, 0.275)
names(th.true) <- tail(rcts, -length(cnstr))
s2.true <- 1e-8
r0.true <- .1
r1.true <- .5
phi.true <- c(th.true, r0.true, r1.true)
names(phi.true) <- c(names(th.true), "r0", "r1")
S <- 1000
nCL <- 3

## initial condition X0:
X0 <- rep(0, length(ctps))
names(X0) <- ctps
X0["HSC"] <- 100

## number of independent simulations:
nSim <- 100
## number of time-points:
ntps <- as.numeric(commandArgs(trailingOnly = TRUE))

## initialize list for simulations results:
resAllSim <- vector("list", length = nSim)

for (sim in 1:nSim) { ## loop over independent simulations
  cat(paste("simulation n. ", sim, "\t\tT = ", ntps, "\n", sep = ""))
  ## simulate trajectories:
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
                             f = 0,
                             ntps = ntps,
                             trunc = FALSE)

  ## start a cluster:
  nProc <- 1
  cat(paste("\tLoading CPU cluster...\n", sep = ""))
  cat(paste("Cluster type: ", "PSOCK\n", sep = ""))
  cpu <- Sys.getenv("SLURM_CPUS_ON_NODE", nProc)
  hosts <- rep("localhost",cpu)
  cl <- makeCluster(hosts, type = "PSOCK")
  rm(nProc)

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
                     Y = XY$Y[,setdiff(ctps, latsts),],
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

  ## scatterplot of estimated parameters against the true ones:
  pdf(file = paste(currResPath, "/sim", sim, "_T", ntps, "_theta.pdf", sep = ""), width = 7, height = 7)
  par(mar = c(5,5,2,2))
  plot(th.true, head(phi.curr,-2),
       xlim = range(th.true, head(phi.curr,-2)),
       ylim = range(th.true, head(phi.curr,-2)),
       xlab = expression(theta[true]),
       ylab = expression(theta[EKF]),
       cex.axis = 2, cex.lab = 2, pch = 21, cex = 5)
  lines(range(th.true, head(phi.curr,-2)),
        range(th.true, head(phi.curr,-2)), lwd = 3, lty = 2, col = "red")
  text(x = th.true,
       y = head(phi.curr,-2),
       labels = names(th.true), font = 2, cex = 2)
  text(x = diff(range(th.true, head(phi.curr,-2)))/10*2, y = max(th.true, head(phi.curr,-2))/10*9, labels = bquote(~ rho[0] == .(round(r0.true,2))), cex = 2)
  text(x = diff(range(th.true, head(phi.curr,-2)))/10*2, y = max(th.true, head(phi.curr,-2))/10*8, labels = bquote(~ widehat(rho)[0] == .(round(phi.curr["r0"],2))), cex = 2)

  text(x = diff(range(th.true, head(phi.curr,-2)))/10*5, y = max(th.true, head(phi.curr,-2))/10*9, labels = bquote(~ rho[1] == .(round(r1.true,2))), cex = 2)
  text(x = diff(range(th.true, head(phi.curr,-2)))/10*5, y = max(th.true, head(phi.curr,-2))/10*8, labels = bquote(~ widehat(rho)[1] == .(round(phi.curr["r1"],2))), cex = 2)
  dev.off()

  ## smoothing moments:
  pdf(file = paste(currResPath, "/sim", sim, "_T", ntps, "_sMoments.pdf", sep = ""), width = 7*ceiling(nCL/ceiling(sqrt(nCL))), height = 3*ceiling(sqrt(nCL)))
  par(mar = c(5,5,2,2), mfrow = c(ceiling(nCL/ceiling(sqrt(nCL))),ceiling(sqrt(nCL))))
  get.sMoments(res.fit = res.fit, X = XY$X)
  dev.off()

  ## clone-average smoothing moments:
  pdf(file = paste(currResPath, "/sim", sim, "_T", ntps, "average_sMoments.pdf", sep = ""), width = 7, height = 5)
  par(mar = c(5,5,2,2), mfrow = c(1,1))
  get.sMoments.avg(res.fit = res.fit, X = XY$X)
  dev.off()

  ## save current results:
  res <- list()
  res$phi.opt <- res.fit$fit$par
  res$bwd.res <- res.fit$bwd.res
  res$Y <- XY$Y
  res$X <- XY$X
  resAllSim[[sim]] <- res
}

## save all the results:
save.image(paste(currResPath, "/output_T", ntps, ".Rdata", sep = ""))


