cat("\nInstall/load packages")
inst.pkgs <- installed.packages()

## required packages:
l.pkgs <- c("Karen")
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

source("./source/Karen-RM-fun.R")

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

macaqueID_mod <- as.vector(unlist(strsplit(commandArgs(trailingOnly = TRUE), split = "-", fixed = TRUE)))
# macaqueID_mod <- as.vector(unlist(strsplit("ZH33-2", split = "-", fixed = TRUE)))
macaqueID <- macaqueID_mod[1]

model.lst <- get.modelList()
nModels <- length(model.lst)

topClones <- 1000
Y <- get.Y(macaqueID, topClones = topClones)
Y <- Y[,,names(head(sort(apply(Y!=0, 3, sum), decreasing = T), topClones)),drop=FALSE]
cat(paste("n. of clones: ", dim(Y)[3], "\n", sep = ""))

nMod <- as.numeric(macaqueID_mod[2])

############### MODEL FITTING #####################
Y0 <- Y
rm(Y)

cat(paste("Macaque ID: ", macaqueID, "\t Model n. ", nMod, "\n", sep = ""))

## cluster parameters:
nProc <- 8
cat(paste("\tLoading CPU cluster...\n", sep = ""))
cat(paste("Cluster type: ", "PSOCK\n", sep = ""))
cpu <- Sys.getenv("SLURM_CPUS_ON_NODE", nProc)
hosts <- rep("localhost",cpu)
cl <- parallel::makeCluster(hosts, type = "PSOCK")
rm(nProc)

## CDN parameters:
rcts <- as.vector(unlist(model.lst[[nMod]]["rct.lst"]))
cnstr <- as.vector(unlist(model.lst[[nMod]]["constr.lst"]))
latsts <- as.vector(unlist(model.lst[[nMod]]["latSts.lst"]))
ctps <- unique(setdiff(c(sapply(rcts, function(r){
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

## initial condition:
X0 <- rep(0, length(ctps))
names(X0) <- ctps
X0["HSC"] <- 1

## mean vector and covariance matrix of X0:
m_0 <- replicate(dim(Y0)[3], X0, simplify = "array")
colnames(m_0) <- dimnames(Y0)[[3]]
P_0 <- Matrix::Diagonal(length(ctps) * dim(Y0)[3], 1e-5)
rownames(P_0) <- colnames(P_0) <- rep(dimnames(Y0)[[3]], each = length(ctps))

## fit Karen on data:
res.fit <- get.fit(rct.lst = rcts,
                   constr.lst = cnstr,
                   latSts.lst = latsts,
                   ct.lst = ctps,
                   Y = Y0,
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
                        verbose = TRUE,
                        FORCEP = TRUE))

## stop cluster and save results:
parallel::stopCluster(cl)
save.image(paste(currResPath, "/rhesusMacaqueID_", macaqueID, "_Model", nMod, "_output.Rdata", sep = ""))

## define color legend for cell types:
cell.cols <-  c("#1F77B4", "#FF7F0E", "#2CA02C", "#E7B928", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F")
names(cell.cols) <- c("HSC", "P1", "P2", "P3", "T", "B", "NK", "G", "M")

## plot smoothing moments:
nCL <- dim(Y0)[3]
pdf(file = paste(currResPath, "/rhesusMacaqueID_", macaqueID, "_Model", nMod, "_sMoments.pdf", sep = ""), width = 7*ceiling(nCL/ceiling(sqrt(nCL))), height = 3*ceiling(sqrt(nCL)))
par(mar = c(5,5,2,2), mfrow = c(ceiling(sqrt(nCL)), ifelse(ceiling(sqrt(nCL)) == 1, ceiling(sqrt(nCL)) + 1, ceiling(sqrt(nCL)))))
get.sMoments(res.fit = res.fit, cell.cols = cell.cols)
dev.off()

## plot clone-average smoothing moments:
pdf(file = paste(currResPath, "/rhesusMacaqueID_", macaqueID, "_Model", nMod, "average_sMoments.pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
get.sMoments.avg(res.fit = res.fit, cell.cols = cell.cols)
dev.off()

# ## load igraphhack functions:
# source_url("https://raw.githubusercontent.com/jevansbio/igraphhack/master/igraphplot2.R")
# environment(plot.igraph2) <- asNamespace('igraph')
# environment(igraph.Arrows2) <- asNamespace('igraph')

## plot cell differentiation network:
legend_image <- grDevices::as.raster(matrix(grDevices::colorRampPalette(c("lightgray", "red", "black"))(99), ncol=1))
pdf(file = paste(currResPath, "/rhesusMacaqueID_",macaqueID, "_Model", nMod, "diffNet.pdf", sep = ""), width = 5, height = 5)
layout(mat = matrix(c(1,1,1,2), ncol = 1))
par(mar = c(1,0,3,0))
get.cdn(res.fit = res.fit,
        edges.lab = F,
        cell.cols = cell.cols)
plot(c(0,1),c(-1,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=seq(0,1,l=5), y = -.2, labels = seq(0,1,l=5), cex = 2, font = 2)
rasterImage(t(legend_image), 0, 0, 1, 1)
dev.off()

## plot data on log-scale:
Y <- res.fit$Y
Y_NA <- Y
Y_NA[Y_NA == 0] <- NA

pdf(file = paste(currResPath, "/rhesusMacaqueID_",macaqueID, "data.pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
matplot(as.numeric(rownames(Y)), log(Y_NA[,,1]),
        lty = 1, pch = 20,
        col = scales::alpha(cell.cols[colnames(Y_NA)], alpha = .2),
        cex = 3, lwd = 3,
        ylim = range(c(log(Y_NA)), na.rm = T), type = 'b',
        cex.axis = 2, cex.lab = 2, ylab = expression(logY[t]), xlab = "time (months)")
lapply(1:dim(Y)[3], function(cl){
  matplot(as.numeric(rownames(Y[,,cl])), log(Y_NA[,,cl]),
          type = "b", lty = 1, pch = 20, add = T,
          col = scales::alpha(cell.cols[colnames(Y_NA)], alpha = .2), cex = 2, lwd = 3)
})
matplot(as.numeric(rownames(Y)), apply(log(Y_NA), c(1,2), mean, na.rm = T), type = 'b', pch = 20, lty = 1,
        col = cell.cols[colnames(Y_NA)], add = T, lwd = 10)
dev.off()
