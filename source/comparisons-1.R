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
            "devtools",
            "Karen",
            "RestoreNet",
            "stringr",
            "ggplot2"
)
## check if packages are installed
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})


lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})

set.seed(123)

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

# rtf <- as.vector(unlist(strsplit(".1-7-.5-3", split = "-", fixed = TRUE)))
rtf <- as.vector(unlist(strsplit(commandArgs(trailingOnly = TRUE), split = "-", fixed = TRUE)))

currDir <- getwd()
setwd(currDir)

resPath <- paste("./results_c1", paste(rtf, collapse = "-"), sep = "_")
currResPath <- paste(resPath, "/", currTime, sep = "")

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(currResPath)), dir.create(file.path(currResPath)), FALSE)

##############################
## update functions of Karen

get.V <- function(ct.lst = ct.lst, rct.lst = rct.lst){
  V <- matrix(data = 0, nrow = length(ct.lst), ncol = length(rct.lst))
  rownames(V) <- ct.lst
  colnames(V) <- rct.lst

  for (r in rct.lst) {
    rgts_prod <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
    V[rgts_prod[1],r] <- -1
    if(rgts_prod[2] != "0"){
      if(rgts_prod[2] == "1"){
        V[rgts_prod[1],r] <- 1
      }else{
        if(rgts_prod[1] != "HSC"){
          V[rgts_prod[1],r] <- 0
        }
        V[rgts_prod[2],r] <- 1
      }
    }
  }
  return(V)
}

assignInNamespace("get.V",get.V,ns="Karen")
Karen:::get.V
rm(get.V)

##################################
## update functions of RestoreNet

get.V <- function (rct.lst)
{
  ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r) {
    as.vector(unlist(str_split(r, "->")))
  }, simplify = "array")))), c("0", "1"))
  V <- matrix(data = 0, nrow = length(ct.lst), ncol = length(rct.lst))
  rownames(V) <- ct.lst
  colnames(V) <- rct.lst
  for (r in rct.lst) {
    rgts_prod <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
    V[rgts_prod[1],r] <- -1
    if(rgts_prod[2] != "0"){
      if(rgts_prod[2] == "1"){
        V[rgts_prod[1],r] <- 1
      }else{
        if(rgts_prod[1] != "HSC"){
          V[rgts_prod[1],r] <- 0
        }
        V[rgts_prod[2],r] <- 1
      }
    }
  }
  return(V)
}

compile.h <- function(rct.lst, envir){
  constr.lst <- c()
  get.h.string <- paste("get.h <- function(x, theta){
                h <- c(",
                        paste(sapply(rct.lst, function(r){
                          rgts <- as.vector(unlist(strsplit(r, split = "->", fixed = T)))[1]
                          prds <- as.vector(unlist(strsplit(r, split = "->", fixed = T)))[2]
                          if(prds == "0"){
                            hx <- paste("x['", rgts, "']*theta['", r, "']", sep = "")
                          }else{
                            hx <- paste("x['", rgts, "']*theta['", r, "']", sep = "")
                          }

                        }, simplify = "array"), collapse = ",\n"),
                        ")
return(h)
}", sep = "")

  for (constr in constr.lst) {
    constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
    pattern <- constr[1]
    replacement <- constr[2]

    get.h.string <- gsub(pattern = pattern,
                         replacement = replacement,
                         x = get.h.string)
  }

  eval(parse(text=get.h.string), envir = envir)
}

assignInNamespace("get.V",get.V,ns="RestoreNet")
assignInNamespace("compile.h",compile.h,ns="RestoreNet")

rm(get.V)
rm(compile.h)

###################
## CDN parameters:
rcts <- c("HSC->1",
          "HSC->P1",
          "HSC->P2",
          "P1->0",
          "P2->0",
          "P1->A",
          "P1->B",
          "P1->C",
          "P2->D",
          "P2->E",
          "A->0",
          "B->0",
          "C->0",
          "D->0",
          "E->0")
cnstr <- c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->A\\'\\] + theta\\[\\'P1->B\\'\\] + theta\\[\\'P1->C\\'\\])",
           "theta\\[\\'HSC->P2\\'\\]=(theta\\[\\'P2->D\\'\\] + theta\\[\\'P2->E\\'\\])")
latsts <- NULL
ctps <- unique(setdiff(c(sapply(rcts, function(r){
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

th.true <- c(NA, NA, 2.2, 2.24, 2.16, 1.8, 1.68, 2.52, 2.232, 2.304, 1.224, 1.512, 6.24, 3.4912, 3.2728)
names(th.true) <- c("HSC->P1",
                    "HSC->P2",
                    "P1->A",
                    "P1->B",
                    "P1->C",
                    "P2->D",
                    "P2->E",
                    "A->0",
                    "B->0",
                    "C->0",
                    "D->0",
                    "E->0"
                    ,"HSC->1"
                    ,"P1->0"
                    ,"P2->0"
)
th.true["HSC->P1"] <- th.true["P1->A"] + th.true["P1->B"] + th.true["P1->C"]
th.true["HSC->P2"] <- th.true["P2->D"] + th.true["P2->E"]

th.true <- th.true[rcts]

s2.true <- 1e-8
r0.true <- r1.true <- as.numeric(rtf[1])

phi.true <- c(th.true, r0.true, r1.true)
names(phi.true) <- c(names(th.true), "r0", "r1")
S <- 1000
nCL <- as.numeric(rtf[4])

## initial condition X0:
X0 <- rep(0, length(ctps))
names(X0) <- ctps
X0["HSC"] <- 100

## number of independent simulations:
nSim <- 100
## number of time-points:
ntps <- as.numeric(rtf[2])
## fraction of latent data:
f_NA <- as.numeric(rtf[3])

allParamsKaren <- matrix(data = NA, nrow = length(phi.true), ncol = nSim)
rownames(allParamsKaren) <- names(phi.true)

allParamsRestoreNet <- matrix(data = NA, nrow = length(phi.true) - 1, ncol = nSim)
rownames(allParamsRestoreNet) <- head(names(phi.true),-1)
rownames(allParamsRestoreNet)[nrow(allParamsRestoreNet)] <- "s2"

allParamsKaren_noCnstr <- matrix(data = NA, nrow = length(phi.true), ncol = nSim)
rownames(allParamsKaren_noCnstr) <- names(phi.true)

for (sim in 1:nSim) { ## loop over independent simulations
  cat(paste("simulation n. ", sim, "\t\tf = ", f_NA, "\n", sep = ""))

  ## simulate trajectories:
  XY <- get.sim.trajectories(rct.lst = rcts,
                             constr.lst = NULL,
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
  Y0 <- XY$Y

  ## Save Data:
  obsTimes <- c(0, as.numeric(rownames(Y0)))
  obsTimes <- (obsTimes - min(obsTimes))/(max(obsTimes) - min(obsTimes))
  rownames(Y0) <- obsTimes <- tail(obsTimes,-1)
  Y0 <- sapply(1:nCL, function(cl){rbind(X0, Y0[,,cl])}, simplify = "array")
  rownames(Y0)[1] <- "0"
  dimnames(Y0)[[3]] <- 1:nCL
  saveRDS(Y0, paste(currResPath, "/sim", sim, "_f", f_NA, "_r", r0.true, "_t", ntps, "-Y0.rds", sep = ""))


  Y0 <- XY$Y

  ## Apply Karen:

  ## mean vector and covariance matrix of X0:
  m_0 <- replicate(dim(Y0)[3], X0, simplify = "array")
  colnames(m_0) <- dimnames(Y0)[[3]]
  P_0 <- Matrix::Diagonal(length(ctps) * dim(Y0)[3], 1) # 1e-5
  rownames(P_0) <- colnames(P_0) <- rep(dimnames(Y0)[[3]], each = length(ctps))

  ## start a cluster:
  nProc <- 1
  cat(paste("\tLoading CPU cluster...\n", sep = ""))
  cat(paste("Cluster type: ", "PSOCK\n", sep = ""))
  cpu <- Sys.getenv("SLURM_CPUS_ON_NODE", nProc)
  hosts <- rep("localhost",cpu)
  cl <- makeCluster(hosts, type = "PSOCK")
  rm(nProc)

  ## fit Karen on data:
  res.fit <- get.fit(rct.lst = rcts,
                     constr.lst = cnstr,
                     latSts.lst = latsts,
                     ct.lst = ctps,
                     Y = Y0,
                     m0 = m_0,
                     P0 = P_0,
                     cl = cl,
                     list(nLQR = 3, ## 3
                          lmm = 25, ## 25
                          pgtol = 0,
                          relErrfct = 1e-4, ## 1e-5
                          tol = 1e-9, ## 1e-9
                          maxit = 1000, ## 1000
                          maxitEM = 3, ## 10
                          trace = 1,
                          verbose = TRUE,
                          FORCEP = FALSE))
  parallel::stopCluster(cl)

  phi.curr <- res.fit$fit$par ## MLE parameters
  allParamsKaren[names(res.fit$fit$par),sim] <- res.fit$fit$par

  ## start a cluster:
  nProc <- 1
  cat(paste("\tLoading CPU cluster...\n", sep = ""))
  cat(paste("Cluster type: ", "PSOCK\n", sep = ""))
  cpu <- Sys.getenv("SLURM_CPUS_ON_NODE", nProc)
  hosts <- rep("localhost",cpu)
  cl <- makeCluster(hosts, type = "PSOCK")
  rm(nProc)

  ## fit Karen without constraints:
  cat("Fit Karen with no constraints...")
  res.fit <- get.fit(rct.lst = rcts,
                     constr.lst = NULL,
                     latSts.lst = latsts,
                     ct.lst = ctps,
                     Y =  XY$Y[,setdiff(ctps, latsts),],
                     m0 = m_0,
                     P0 = P_0,
                     cl = cl,
                     list(nLQR = 3, ## 3
                          lmm = 25, ## 25
                          pgtol = 0,
                          relErrfct = 1e-4, ## 1e-5
                          tol = 1e-9, ## 1e-9
                          maxit = 1000, ## 1000
                          maxitEM = 3, ## 10
                          trace = 1,
                          verbose = TRUE,
                          FORCEP = FALSE))
  parallel::stopCluster(cl)
  cat("done\n")

  phi.curr <- res.fit$fit$par ## MLE parameters
  allParamsKaren_noCnstr[names(res.fit$fit$par),sim] <- res.fit$fit$par

  ####################
  ## Apply RestoreNet:

  obsTimes <- c(0, as.numeric(rownames(Y0)))
  obsTimes <- (obsTimes - min(obsTimes))/(max(obsTimes) - min(obsTimes))
  rownames(Y0) <- obsTimes <- tail(obsTimes,-1)

  Y0 <- Y0[,,as.numeric(which(apply(Y0, 3, function(Ycl){sum(Ycl) > 0})))]

  res.null <- fit.null(Y = Y0,
                       lmm = 25,
                       factr = (1e-5/.Machine$double.eps),
                       pgtol = 0,
                       rct.lst = rcts)

  phi.curr <- res.null$fit$par

  allParamsRestoreNet[,sim] <- res.null$fit$par

}

allParamsKaren[c("HSC->P1", "HSC->P2"),] <- apply(allParamsKaren, 2, FUN = function(x){
  c(x["P1->A"] + x["P1->B"] + x["P1->C"],
    x["P2->D"] + x["P2->E"])
})

## save all the results:
save.image(paste(currResPath, "/output_f", f_NA, "_r", r0.true, "_t", ntps, ".Rdata", sep = ""))

allParamsGLS <- t(read.table(paste0("./GLS-results/c1_", paste(rtf, collapse = "-"), "_GLS"), sep = "\t", header = F)[,1:length(th.true)])

all.res <- data.frame(variable = rep(rownames(head(allParamsKaren,-2)), times = 4*nSim),
                      value = as.numeric(c(c(head(allParamsKaren,-2)/th.true),
                                           c(head(allParamsKaren_noCnstr,-2)/th.true),
                                           c(head(allParamsRestoreNet,-1)/th.true),
                                           c(allParamsGLS/th.true))),
                      method = c(rep("Karen", each = length(th.true)*nSim),
                                 rep("Karen_noCnstr", each = length(th.true)*nSim),
                                 rep("RestoreNet", each = length(th.true)*nSim),
                                 rep("GLS", each = length(th.true)*nSim)))

# methods color
methodsColors <- c(branchCorr = "#FFCC00",
                   GLS = "#F8766D",
                   RestoreNet = "#00BA38",
                   Karen = "#00BFC4",
                   Karen_noCnstr = "#C77CFF")

# New sorting order
methods_order <- c("GLS", 
                   "RestoreNet", 
                   "Karen", 
                   "Karen_noCnstr")
# Re-order the levels
all.res$method <- factor(as.character(all.res$method), levels=methods_order)
# Re-order the data.frame
all.res <- all.res[order(all.res$method),]

yLim <- c(+Inf,-Inf)
for (m in unique(all.res$method)) {
  for (v in unique(all.res$variable)) {
    curr.whiskers <- boxplot.stats(all.res[which(all.res$method == m &
                                                   all.res$variable == v), "value"])$stats[c(1, 5)]
    yLim[1] <- ifelse(curr.whiskers[1] <= yLim[1], curr.whiskers[1], yLim[1])
    yLim[2] <- ifelse(curr.whiskers[2] >= yLim[2], curr.whiskers[2], yLim[2])
  }
}

pdf(file = paste(currResPath, "/allBoxplots_f", f_NA, "_r", r0.true, "_t", ntps, "_theta.pdf", sep = ""), width = 15, height = 5)
(ggplot(all.res, aes(x = variable, y = value, fill = method)) +  # ggplot function
    geom_boxplot(position = position_dodge(preserve = "single"),outlier.shape = NA) +
    scale_fill_manual(values=methodsColors[unique(as.character(all.res$method))]) +
    labs(x = "parameter", y = "estimates") +
    theme(axis.title.y = element_text(size = rel(2), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(2), angle = 00)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20,color="black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size=20,color="black")) +
    theme(legend.text = element_text(size=30)) +
    theme(legend.title = element_text(size=30)) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept=1, color = "red", size=1.5) +
    coord_cartesian(ylim = yLim)
)
dev.off()



