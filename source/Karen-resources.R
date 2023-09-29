
## required packages:
library(Karen)
library(gaussquad)
library(parallel)
library(Matrix)

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


## initial condition X0:
X0 <- rep(0, length(ctps))
names(X0) <- ctps
X0["HSC"] <- 100

ntps_fNA_nCL <- commandArgs(trailingOnly = TRUE)
# ntps_fNA_nCL <- "5-.3-3-1"

## number of time-points:
ntps <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[1])
## fraction of latent data:
f_NA <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[2])
## n. of clones
nCL <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[3])
## n. of cores:
nCores <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[4])

cat(paste("n. of time points: ", ntps, "\tmissing data: ", paste0(f_NA*100, "%"), "\tn. of clones: ", nCL, "\tn. of cores: ", nCores, "\n", sep = ""))
## simulate trajectories:
set.seed(123)
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

## mean vector and covariance matrix of X0:
m_0 <- replicate(nCL, X0, simplify = "array")
colnames(m_0) <- 1:nCL
P_0 <- Diagonal(length(ctps) * nCL, 1e-5)
rownames(P_0) <- colnames(P_0) <- rep(1:nCL, each = length(ctps))

## start a cluster:
cl <- makeCluster(nCores, type = "FORK")
clusterEvalQ(cl, library("Matrix"))

gc0 <- gc(full=TRUE) # set the zero-level memory
## fit Karen on data:
t.start <- proc.time()
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
                        tol = 1e-4,
                        maxit = 10000,
                        maxitEM = 3,
                        trace = 1,
                        verbose = TRUE,
                        FORCEP = FALSE))
t.end <- proc.time()
stopCluster(cl) ## stop the cluster
gc1 <- gc(full=TRUE) # memory used for inference

########### STORE RESULTS on RESOURCES USAGE ###########
res.usage.mat <- cbind(rbind(gc1 - gc0,
                             total = colSums(gc1 - gc0)),
                       rbind(matrix(0,2,3),
                             head(t.end - t.start,3)))

res.usage.tex <- xtable::xtable(res.usage.mat)
write.table(res.usage.mat, file = paste(currResPath, "/res-usage-f_", f_NA, "-ntps_", ntps, "-nCl_", nCL, "-nCores_", nCores, ".tsv", sep = ""),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
print(res.usage.tex, file = paste(currResPath, "/res-usage-f_", f_NA, "-ntps_", ntps, "-nCl_", nCL, "-nCores_", nCores, ".txt", sep = ""))

phi.curr <- res.fit$fit$par## MLE parameters
phi.curr[phi.curr == 0] <- 1e-8

## scatterplot of estimated parameters against the true ones:
pdf(file = paste(currResPath, "/theta-f_", f_NA, "-ntps_", ntps, "-nCl_", nCL, "-nCores_", nCores, ".pdf", sep = ""), width = 7, height = 7)
par(mar = c(5,5,2,2))
plot(log(th.true), log(head(phi.curr,-2)),
     xlim = log(range(th.true, head(phi.curr,-2))),
     ylim = log(range(th.true, head(phi.curr,-2))),
     xlab = expression(theta[true]),
     ylab = expression(theta[EKF]),
     cex.axis = 2, cex.lab = 2, pch = 21, cex = 5)
lines(log(range(th.true, head(phi.curr,-2))),
      log(range(th.true, head(phi.curr,-2))), lwd = 3, lty = 2, col = "red")
text(x = log(th.true),
     y = log(head(phi.curr,-2)),
     labels = names(th.true), font = 2, cex = 2)
text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.8, y = max(log(th.true), log(head(phi.curr,-2)))*.9, labels = bquote(~ rho[0] == .(round(r0.true,2))), cex = 2)
text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.8, y = max(log(th.true), log(head(phi.curr,-2)))*.7, labels = bquote(~ widehat(rho)[0] == .(round(phi.curr["r0"],2))), cex = 2)

text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.2, y = max(log(th.true), log(head(phi.curr,-2)))*.9, labels = bquote(~ rho[1] == .(round(r1.true,2))), cex = 2)
text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.2, y = max(log(th.true), log(head(phi.curr,-2)))*.7, labels = bquote(~ widehat(rho)[1] == .(round(phi.curr["r1"],2))), cex = 2)
dev.off()

## smoothing moments:
pdf(file = paste(currResPath, "/sMoments-f_", f_NA, "-ntps_", ntps, "-nCl_", nCL, "-nCores_", nCores, ".pdf", sep = ""), width = 7*ceiling(nCL/ceiling(sqrt(nCL))), height = 3*ceiling(sqrt(nCL)))
par(mar = c(5,5,2,2), mfrow = c(ceiling(nCL/ceiling(sqrt(nCL))),ceiling(sqrt(nCL))))
get.sMoments(res.fit = res.fit, X = XY$X)
dev.off()

## clone-average smoothing moments:
pdf(file = paste(currResPath, "/average-sMoments-f_", f_NA, "-ntps_", ntps, "-nCl_", nCL, "-nCores_", nCores, ".pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
get.sMoments.avg(res.fit = res.fit, X = XY$X)
dev.off()

## save all the results:
save.image(paste(currResPath, "/output-f_", f_NA, "-ntps_", ntps, "-nCl_", nCL, "-nCores_", nCores, ".Rdata", sep = ""))


