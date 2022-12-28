get.modelList <- function(){
  model.lst <- list()

  model.lst[[1]] <- list(rct.lst = c("HSC->P1",
                                     "P1->T",
                                     "P1->B",
                                     "P1->NK",
                                     "P1->G",
                                     "P1->M",
                                     "T->0",
                                     "B->0",
                                     "NK->0",
                                     "G->0",
                                     "M->0",
                                     "HSC->1",
                                     "P1->1"),
                         constr.lst = c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->T\\'\\] + theta\\[\\'P1->B\\'\\] + theta\\[\\'P1->NK\\'\\] + theta\\[\\'P1->G\\'\\] + theta\\[\\'P1->M\\'\\])"),
                         latSts.lst = c("HSC", "P1"))

  model.lst[[2]] <- list(rct.lst = c("HSC->P1",
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
                                     ,"P2->1"),
                         constr.lst = c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->T\\'\\] + theta\\[\\'P1->B\\'\\] + theta\\[\\'P1->NK\\'\\])",
                                        "theta\\[\\'HSC->P2\\'\\]=(theta\\[\\'P2->G\\'\\] + theta\\[\\'P2->M\\'\\])"),
                         latSts.lst = c("HSC", "P1", "P2"))

  model.lst[[3]] <- list(rct.lst = c("HSC->P1",
                                     "HSC->P2",
                                     "P1->T",
                                     "P1->B",
                                     "P1->NK",
                                     "P2->NK",
                                     "P2->G",
                                     "P2->M",
                                     "T->0",
                                     "B->0",
                                     "NK->0",
                                     "G->0",
                                     "M->0"
                                     ,"HSC->1"
                                     ,"P1->1"
                                     ,"P2->1"),
                         constr.lst = c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->T\\'\\] + theta\\[\\'P1->B\\'\\] + theta\\[\\'P1->NK\\'\\])",
                                        "theta\\[\\'HSC->P2\\'\\]=(theta\\[\\'P2->NK\\'\\] + theta\\[\\'P2->G\\'\\] + theta\\[\\'P2->M\\'\\])"),
                         latSts.lst = c("HSC", "P1", "P2"))

  model.lst[[4]] <- list(rct.lst = c("HSC->P1",
                                     "HSC->P2",
                                     "HSC->P3",
                                     "P1->T",
                                     "P1->B",
                                     "P3->NK",
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
                                     ,"P3->1"),
                         constr.lst = c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->T\\'\\] + theta\\[\\'P1->B\\'\\])",
                                        "theta\\[\\'HSC->P2\\'\\]=(theta\\[\\'P2->G\\'\\] + theta\\[\\'P2->M\\'\\])",
                                        "theta\\[\\'HSC->P3\\'\\]=theta\\[\\'P3->NK\\'\\]"),
                         latSts.lst = c("HSC", "P1", "P2", "P3"))

  return(model.lst)
}

get.Y <- function(macaqueID, topClones){
  if(macaqueID == "all"){
    Y_all <- lapply(c("ZH33", "ZH17", "ZG66"), function(mID){
      Y.ID <- get.Y.ID(mID, topClones = ifelse(mID == "ZH33", round(topClones/3) + 1, round(topClones/3)))
      dimnames(Y.ID)[[3]] <- paste(dimnames(Y.ID)[[3]], mID, sep = "-")
      return(Y.ID)
    })
    tps_all <- sort(unique(as.vector(unlist(lapply(Y_all, function(Y){as.numeric(rownames(Y))})))))
    Y <- array(data = 0, dim = c(length(tps_all),
                                 ncol(Y_all[[1]]),
                                 sum(as.vector(unlist(lapply(Y_all, function(Y){dim(Y)[[3]]}))))),
               dimnames = list(tps_all,
                               colnames(Y_all[[1]]),
                               as.vector(unlist(lapply(Y_all, function(Y){dimnames(Y)[[3]]})))))
    for (nID in 1:length(Y_all)) {
      Y[dimnames(Y_all[[nID]])[[1]],
        dimnames(Y_all[[nID]])[[2]],
        dimnames(Y_all[[nID]])[[3]]] <- Y_all[[nID]]
    }
  }else{
    Y <- get.Y.ID(macaqueID, topClones = topClones)
  }
  return(Y)
}

get.Y.ID <- function(macaqueID, topClones){

  data(Y_RM)
  Y <- Y_RM[[macaqueID]]
  min_ab <- min(sapply(c("ZH33", "ZH17", "ZG66"), function(macaqueID){
    Y <- Y_RM[[macaqueID]]
    min_ab <- min(setdiff(c(apply(Y, c(1,2), sum)), 0))
    return(min_ab)
  }))

  f <- min_ab/apply(Y, c(1,2), sum)
  f[is.infinite(f)] <- 0
  Y_rescaled <- Y
  Y_rescaled <- sapply(1:dim(Y)[3], function(cl){round(Y[,,cl]*f)}, simplify = "array")
  dimnames(Y_rescaled)[3] <- dimnames(Y)[3]
  Y <- Y_rescaled
  Y <- Y[,,names(which(apply(Y!=0, 3, sum) != 0))]
  # topClones <- dim(Y)[3] # 1000
  Y <- Y[,,names(head(sort(apply(Y!=0, 3, sum), decreasing = T), topClones)),drop=FALSE]

  return(Y)
}

