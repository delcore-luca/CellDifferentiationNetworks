
## required packages:
library(Karen)
library(gaussquad)
library(parallel)
library(Matrix)
library(scales)

get.sMoments <- function(res.fit, X = NULL, cell.cols = NULL){
  V <- res.fit$V # net-effect matrix
  nProc <- length(res.fit$cloneChunks) # nrow(summary(cl)) # number of cores
  Y <- res.fit$Y # simulated measurements
  Y_NA <- Y
  Y_NA[Y_NA == 0] <- NA

  if(!is.null(cell.cols)){
    cols <- cell.cols[rownames(V)]
  }else{
    cols <- palette.colors(nrow(V), palette = "Classic Tableau")
  }

  if(!is.null(X)){
    tps <- as.numeric(rownames(X))
    rownames(X) <- (tps - min(tps))/(max(tps) - min(tps))
  }

  lapply(1:nProc, function(cnk){
    lapply(1:length(res.fit$cloneChunks[[cnk]]), function(cl){
      idx.clones <- matrix(data = 1:(length(res.fit$cloneChunks[[cnk]])*nrow(V)), nrow = nrow(V), ncol = length(res.fit$cloneChunks[[cnk]]))
      mean_smooth <- t(res.fit$bwd.res$m_xt_Yn[[cnk]][,cl,])
      sd_smooth <- sqrt(t(sapply(1:(dim(res.fit$bwd.res$V_xt_Yn[[cnk]][idx.clones[,cl],idx.clones[,cl],])[3]), FUN = function(t){diag(nearestPD(res.fit$bwd.res$V_xt_Yn[[cnk]][idx.clones[,cl],idx.clones[,cl],t]))})))
      rownames(sd_smooth) <- rownames(mean_smooth)

      matplot(as.numeric(rownames(Y_NA)), Y_NA[,,res.fit$cloneChunks[[cnk]][cl]], lty = 1, pch = 20, type = 'p', add = F, col = alpha(cell.cols, alpha = .8), cex = 2,
              cex.axis = 2, cex.lab = 2, xlab = "t", ylab = expression("Y"[t]), main = paste("clone ", cl, sep = ""), cex.main = 2,
              xlim = c(0, max(as.numeric(rownames(Y)))), ylim = range(c(Y_NA[,,res.fit$cloneChunks[[cnk]][cl]], mean_smooth, mean_smooth - 1.96*sd_smooth, mean_smooth + 1.96*sd_smooth), na.rm = T))
      if(!is.null(X)){
        matplot(as.numeric(rownames(X)), X[,,res.fit$cloneChunks[[cnk]][cl]], add = T, pch = 1, cex = 1.5, lwd = 2, col = alpha(cell.cols, alpha = .8))
      }
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth, lwd = 2, lty = 1, type = 'l', add = TRUE, col = cell.cols)
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth - 1.96*sd_smooth, lwd = 2, lty = 3, type = 'l', add = TRUE, col = cell.cols)
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth + 1.96*sd_smooth, lwd = 2, lty = 3, type = 'l', add = TRUE, col = cell.cols)
      # legend(x = "topright", legend = rownames(V), col = cell.cols, pch = 20, lwd = 5, cex = 1.5)
    })
  })
  plot.new()
  legend(x = "center", legend = rownames(V), col = cell.cols, pch = 20, lwd = 5, cex = 1.5)
}

check.cnstrs <- function(rct.lst,
                         constr.lst){}
assignInNamespace("check.cnstrs",check.cnstrs,ns="Karen")

ntps_fNA_nCL <- commandArgs(trailingOnly = TRUE)
# ntps_fNA_nCL <- "5-.8-100-10-20"

## number of time-points:
ntps <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[1])
## fraction of latent data:
f_NA <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[2])
## n. of clones
nCL <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[3])
## n. of cores:
nCores <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[4])
# ## n. of PB cell types:
nSim <- as.numeric(as.vector(unlist(strsplit(ntps_fNA_nCL, split = "-")))[5])

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

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(currResPath)), dir.create(file.path(currResPath)), FALSE)

## CDN parameters:
rcts <- c("HSC->MPP",
          "MPP->CMP",
          "MPP->CLP",
          "CMP->MEP",
          "CMP->GMP",
          "MEP->P",
          "MEP->ERY",
          "GMP->G",
          "GMP->M",
          "CLP->NK",
          "CLP->B",
          "CLP->T",
          "HSC->1",
          "MPP->1",
          "CMP->1",
          "CLP->1",
          "MEP->1",
          "GMP->1",
          "P->0",
          "ERY->0",
          "G->0",
          "M->0",
          "T->0",
          "B->0",
          "NK->0")

cnstr <- c(
  "theta\\[\\'CMP->MEP\\'\\]=(theta\\[\\'MEP->P\\'\\] + theta\\[\\'MEP->ERY\\'\\])",
  "theta\\[\\'CMP->GMP\\'\\]=(theta\\[\\'GMP->G\\'\\] + theta\\[\\'GMP->M\\'\\])",
  "theta\\[\\'MPP->CMP\\'\\]=(theta\\[\\'MEP->P\\'\\] + theta\\[\\'MEP->ERY\\'\\] + theta\\[\\'GMP->G\\'\\] + theta\\[\\'GMP->M\\'\\])",
  "theta\\[\\'MPP->CLP\\'\\]=(theta\\[\\'CLP->T\\'\\] + theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->NK\\'\\])",
  "theta\\[\\'HSC->MPP\\'\\]=(theta\\[\\'MEP->P\\'\\] + theta\\[\\'MEP->ERY\\'\\] + theta\\[\\'GMP->G\\'\\] + theta\\[\\'GMP->M\\'\\] + theta\\[\\'CLP->T\\'\\] + theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->NK\\'\\])"
)
latsts <- c("HSC", "MPP", "CMP", "CLP", "MEP", "GMP")
# latsts <- NULL

ctps <- unique(setdiff(c(sapply(rcts, function(r){
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

## true dynamic parameters:
th.true <- c(c(0.4050, 1.0080/3, 1.0125/3, 1.1475/3, 0.42, .56, .61)*2, c(0.91, .81, .63, .73, .71, .32)*1, c(2.7750, 2.67, 2.81, 2.91, 11.5, 2.5400, 4.9100)/2)*2
names(th.true) <- tail(rcts, -length(cnstr))

s2.true <- 1e-8
r0.true <- .1
r1.true <- .1
phi.true <- c(th.true, r0.true, r1.true)
names(phi.true) <- c(names(th.true), "r0", "r1")
S <- 1000

## initial condition X0:
X0 <- rep(0, length(ctps))
names(X0) <- ctps
X0["HSC"] <- 100

## initialize list for simulations results:
resAllSim <- vector("list", length = nSim)

allParamsKaren <- matrix(data = NA, nrow = length(phi.true), ncol = nSim)
rownames(allParamsKaren) <- names(phi.true)

## define color legend for cell types:
cell.cols <-  c("#DBC9D8", "#7EEB4F", "#726562", "#7CE6CB", "#8FDD8A", "#D556D2", "#8294D9", "#D6DFB6", "#DE687A", "#76908A", "#D6A165", "#DBE261",
                "#89D3E4", "#D492D4", "#7A50E1")
names(cell.cols) <- c("HSC", "MPP", "CMLP", "CMP", "NKP", "CLP", "MEP", "GMP", "P", "ERY", "G", "M", "NK", "T", "B")
cell.cols <- cell.cols[ctps]

set.seed(123)
for (sim in 1:nSim) { ## loop over independent simulations
  cat(paste("simulation n. ", sim, "\n", sep = ""))
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
                             f = f_NA,
                             ntps = ntps,
                             trunc = FALSE)

  cat(paste("\tLoading CPU cluster...\n", sep = ""))
  cat(paste("Cluster type: ", "FORK\n", sep = ""))
  cl <- makeCluster(nCores, type = "FORK")

  ## mean vector and covariance matrix of X0:
  m_0 <- replicate(nCL, X0, simplify = "array")
  colnames(m_0) <- 1:nCL
  P_0 <- Diagonal(length(ctps) * nCL, 100)
  rownames(P_0) <- colnames(P_0) <- rep(1:nCL, each = length(ctps))

  gc0 <- gc(full=TRUE) # set the zero-level memory
  t.start <- proc.time()
  ## fit Karen on data:
  res.fit <- get.fit(rct.lst = rcts,
                     constr.lst = cnstr,
                     latSts.lst = latsts,
                     # latSts.lst = NULL,
                     ct.lst = ctps,
                     Y = XY$Y[,setdiff(ctps, latsts),],
                     # Y = Y0,
                     m0 = m_0,
                     P0 = P_0,
                     cl = cl,
                     list(nLQR = 3,
                          lmm = 25,
                          pgtol = 0,
                          relErrfct = 1e-5,
                          tol = 1e-5,
                          maxit = 1000,
                          maxitEM = 3,
                          trace = 1,
                          verbose = TRUE,
                          FORCEP = FALSE))
  t.end <- proc.time()
  stopCluster(cl) ## stop the cluster
  gc1 <- gc(full=TRUE) # memory used for inference
  phi.curr <- res.fit$fit$par ## MLE parameters
  phi.curr[phi.curr == 0] <- 1e-8

  res.usage.mat <- cbind(rbind(gc1 - gc0,
                               total = colSums(gc1 - gc0)),
                         rbind(matrix(0,2,3),
                               head(t.end - t.start,3)))

  ## scatterplot of estimated parameters against the true ones:
  pdf(file = paste(currResPath, "/theta-sim", sim, ".pdf", sep = ""), width = 7, height = 7)
  par(mar = c(5,5,2,2))
  plot(log(th.true), log(head(phi.curr,-2)),
       xlim = log(range(th.true, head(phi.curr,-2))),
       ylim = log(range(th.true, head(phi.curr,-2))),
       xlab = expression(log~theta[true]),
       ylab = expression(log~theta[EKF]),
       cex.axis = 2, cex.lab = 2, pch = 20, cex = 2)
  lines(log(range(th.true, head(phi.curr,-2))),
        log(range(th.true, head(phi.curr,-2))), lwd = 3, lty = 2, col = "red")
  text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.8, y = max(log(th.true), log(head(phi.curr,-2)))*.9, labels = bquote(~ rho[0] == .(round(r0.true,2))), cex = 2)
  text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.8, y = max(log(th.true), log(head(phi.curr,-2)))*.6, labels = bquote(~ widehat(rho)[0] == .(round(phi.curr["r0"],2))), cex = 2)

  text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.2, y = max(log(th.true), log(head(phi.curr,-2)))*.9, labels = bquote(~ rho[1] == .(round(r1.true,2))), cex = 2)
  text(x = min(c(log(th.true), log(head(phi.curr,-2))))*.2, y = max(log(th.true), log(head(phi.curr,-2)))*.6, labels = bquote(~ widehat(rho)[1] == .(round(phi.curr["r1"],2))), cex = 2)
  dev.off()

  ## smoothing moments:
  pdf(file = paste(currResPath, "/sim", sim, "_T", ntps, "_sMoments.pdf", sep = ""), width = 7*ceiling(nCL/ceiling(sqrt(nCL))), height = 3*ceiling(sqrt(nCL)))
  par(mar = c(5,5,2,2), mfrow = c(ceiling(nCL/ceiling(sqrt(nCL))),ceiling(sqrt(nCL))))
  get.sMoments(res.fit = res.fit, X = XY$X, cell.cols = cell.cols)
  dev.off()

  ## clone-average smoothing moments:
  pdf(file = paste(currResPath, "/sim", sim, "_T", ntps, "average_sMoments.pdf", sep = ""), width = 7, height = 5)
  par(mar = c(5,5,2,2), mfrow = c(1,1))
  get.sMoments.avg(res.fit = res.fit, X = XY$X, cell.cols = cell.cols)
  dev.off()

  ## save current results:
  res <- list()
  res$phi.opt <- res.fit$fit$par
  res$bwd.res <- res.fit$bwd.res
  res$Y <- XY$Y
  res$X <- XY$X
  res$res.usage <- res.usage.mat
  resAllSim[[sim]] <- res

  allParamsKaren[names(res.fit$fit$par),sim] <- phi.curr

  pdf(file = paste(currResPath, "/boxplots-partial.pdf", sep = ""), width = 8, height = 4)
  par(mar = c(4,3,1,1))
  boxplot(log(t(head(allParamsKaren[,1:sim],-2)/th.true)), xaxt = 'n', cex.axis = 1.5, cex.lab = 1.5)
  abline(h = 0, col = "red", lwd = 2)
  axis(side = 1, at = 1:length(th.true),labels = FALSE)
  text(x = 1:length(th.true) + .5,
       y = par("usr")[3] - 0.45,
       labels = names(th.true),
       pos = 2,
       xpd = NA,
       adj = 1,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       cex = 1,
       font = 2)
  dev.off()

  gc()
}

## save all the results:
save.image(paste(currResPath, "/output.Rdata", sep = ""))

pdf(file = paste(currResPath, "/boxplots.pdf", sep = ""), width = 8, height = 4)
par(mar = c(4,3,1,1))
boxplot(t(head(allParamsKaren,-2)/th.true),
        xaxt = 'n',
        cex.axis = 1.5,
        cex.lab = 1.5,
        outline = FALSE)
abline(h = 1, col = "red", lwd = 2)
axis(side = 1, at = 1:length(th.true),labels = FALSE)
text(x = 1:length(th.true) + .5,
     y = par("usr")[3] - 0.45,
     labels = names(th.true),
     pos = 2,
     xpd = NA,
     adj = 1,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1,
     font = 2)
dev.off()

