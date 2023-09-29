get.modelList <- function(){
  model.lst <- list()

  model.lst[[1]] <- list(rct.lst = c("HSC->MPP",
                                     "MPP->CMLP",
                                     "CMLP->MEP",
                                     "CMLP->GMP",
                                     "MEP->P",
                                     "MEP->ERY",
                                     "GMP->G",
                                     "GMP->M",
                                     "CMLP->NK",
                                     "CMLP->B",
                                     "CMLP->T",
                                     "HSC->1",
                                     "MPP->1",
                                     "CMLP->1",
                                     "MEP->1",
                                     "GMP->1",
                                     "P->0",
                                     "ERY->0",
                                     "G->0",
                                     "M->0",
                                     "NK->0",
                                     "B->0",
                                     "T->0"),
                         constr.lst = c(
                           "theta\\[\\'MEP->P\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'MEP->ERY\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'GMP->G\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'CMLP->MEP\\'\\]=2*(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'CMLP->GMP\\'\\]=2*(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'CMLP->NK\\'\\]=(theta\\[\\'CMLP->B\\'\\] + theta\\[\\'CMLP->T\\'\\])/2",
                           "theta\\[\\'MPP->CMLP\\'\\]=(4*(theta\\[\\'GMP->M\\'\\]) + 3*(theta\\[\\'CMLP->B\\'\\] + theta\\[\\'CMLP->T\\'\\])/2)",
                           "theta\\[\\'HSC->MPP\\'\\]=(4*(theta\\[\\'GMP->M\\'\\]) + 3*(theta\\[\\'CMLP->B\\'\\] + theta\\[\\'CMLP->T\\'\\])/2)",
                           ##
                           "theta\\[\\'P->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'ERY->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'G->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'NK->0\\'\\]=(theta\\[\\'B->0\\'\\] + theta\\[\\'T->0\\'\\])/2",
                           "theta\\[\\'MEP->1\\'\\]=theta\\[\\'GMP->1\\'\\]"
                         ),
                         latSts.lst = c("HSC", "MPP", "CMLP", "MEP", "GMP", "P", "ERY", "G", "NK"))

  model.lst[[2]] <- list(rct.lst = c("HSC->MPP",
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
                                     "NK->0"),
                         constr.lst = c(
                           "theta\\[\\'MEP->P\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'MEP->ERY\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'GMP->G\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           ##
                           "theta\\[\\'CMP->MEP\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'CMP->GMP\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'MPP->CMP\\'\\]=(4*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'CLP->NK\\'\\]=(theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->T\\'\\])/2",
                           "theta\\[\\'MPP->CLP\\'\\]=3*(theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->T\\'\\])/2",
                           "theta\\[\\'HSC->MPP\\'\\]=(4*(theta\\[\\'GMP->M\\'\\]) + 3*(theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->T\\'\\])/2)",
                           ##
                           "theta\\[\\'P->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'ERY->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'G->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'NK->0\\'\\]=(theta\\[\\'B->0\\'\\] + theta\\[\\'T->0\\'\\])/2",
                           "theta\\[\\'MEP->1\\'\\]=(theta\\[\\'GMP->1\\'\\])"
                         ),
                         latSts.lst = c("HSC", "MPP", "CMP", "CLP", "MEP", "GMP", "P", "ERY", "G", "NK"))

  model.lst[[3]] <- list(rct.lst = c("HSC->MPP",
                                     "MPP->CMP",
                                     "MPP->CLP",
                                     "CMP->MEP",
                                     "CMP->GMP",
                                     "MEP->P",
                                     "MEP->ERY",
                                     "GMP->G",
                                     "GMP->M",
                                     "CMP->NK",
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
                                     "NK->0"),
                         constr.lst = c(
                           "theta\\[\\'MEP->P\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'MEP->ERY\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'GMP->G\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           ##
                           "theta\\[\\'CMP->MEP\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'CMP->GMP\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'CMP->NK\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'MPP->CMP\\'\\]=(6*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'CLP->NK\\'\\]=(theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->T\\'\\])/2",
                           "theta\\[\\'MPP->CLP\\'\\]=(3*(theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->T\\'\\])/2)",
                           "theta\\[\\'HSC->MPP\\'\\]=(6*(theta\\[\\'GMP->M\\'\\]) + 3*(theta\\[\\'CLP->B\\'\\] + theta\\[\\'CLP->T\\'\\])/2)",
                           ##
                           "theta\\[\\'P->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'ERY->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'G->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'NK->0\\'\\]=(theta\\[\\'B->0\\'\\] + theta\\[\\'T->0\\'\\])/2",
                           "theta\\[\\'MEP->1\\'\\]=(theta\\[\\'GMP->1\\'\\])"
                         ),
                         latSts.lst = c("HSC", "MPP", "CMP", "CLP", "MEP", "GMP", "P", "ERY", "G", "NK"))

  model.lst[[4]] <- list(rct.lst = c("HSC->MPP",
                                     "MPP->CMP",
                                     "MPP->NKP",
                                     "MPP->CLP",
                                     "CMP->MEP",
                                     "CMP->GMP",
                                     "MEP->P",
                                     "MEP->ERY",
                                     "GMP->G",
                                     "GMP->M",
                                     "NKP->NK",
                                     "CLP->B",
                                     "CLP->T",
                                     "HSC->1",
                                     "MPP->1",
                                     "CMP->1",
                                     "NKP->1",
                                     "CLP->1",
                                     "MEP->1",
                                     "GMP->1",
                                     "P->0",
                                     "ERY->0",
                                     "G->0",
                                     "M->0",
                                     "T->0",
                                     "B->0",
                                     "NK->0"),
                         constr.lst = c(
                           "theta\\[\\'MEP->P\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'MEP->ERY\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           "theta\\[\\'GMP->G\\'\\]=(theta\\[\\'GMP->M\\'\\])",
                           ##
                           "theta\\[\\'CMP->MEP\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'CMP->GMP\\'\\]=(2*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'MPP->CMP\\'\\]=(4*(theta\\[\\'GMP->M\\'\\]))",
                           "theta\\[\\'NKP->NK\\'\\]=(theta\\[\\'CLP->T\\'\\] + theta\\[\\'CLP->B\\'\\])/2",
                           "theta\\[\\'MPP->NKP\\'\\]=(theta\\[\\'CLP->T\\'\\] + theta\\[\\'CLP->B\\'\\])/2",
                           "theta\\[\\'MPP->CLP\\'\\]=(theta\\[\\'CLP->T\\'\\] + theta\\[\\'CLP->B\\'\\])",
                           "theta\\[\\'HSC->MPP\\'\\]=(4*(theta\\[\\'GMP->M\\'\\]) + 3*(theta\\[\\'CLP->T\\'\\] + theta\\[\\'CLP->B\\'\\])/2)",
                           ##
                           "theta\\[\\'P->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'ERY->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'G->0\\'\\]=(theta\\[\\'M->0\\'\\])",
                           "theta\\[\\'NK->0\\'\\]=(theta\\[\\'B->0\\'\\] + theta\\[\\'T->0\\'\\])/2",
                           "theta\\[\\'MEP->1\\'\\]=(theta\\[\\'GMP->1\\'\\])",
                           "theta\\[\\'NKP->1\\'\\]=(theta\\[\\'CLP->1\\'\\])"
                         ),
                         latSts.lst = c("HSC", "MPP", "CMP", "NKP", "CLP", "MEP", "GMP", "P", "ERY", "G", "NK"))

  return(model.lst)
}

check.cnstrs <- function(rct.lst,
                         constr.lst){}
assignInNamespace("check.cnstrs",check.cnstrs,ns="Karen")

get.cdn2 <- function(res.fit, edges.lab = FALSE, AIC = FALSE, cell.cols = NULL, modelStructure){

  if(!is.null(cell.cols)){
    cols <- cell.cols[rownames(res.fit$V)]
  }else{
    cols <- palette.colors(nrow(res.fit$V), palette = "Classic Tableau")
    names(cols) <- rownames(res.fit$V)
  }

  mycircle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size  <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
    }

    mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
           vertex.size, vertex.frame.width,
           FUN=function(x, y, bg, fg, size, lwd) {
             symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                     circles=size, add=TRUE, inches=FALSE)
           })
  }

  add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                   plot=mycircle, parameters=list(vertex.frame.color=1,
                                                  vertex.frame.width=1))

  phi.curr <- res.fit$fit$par
  phi.curr <- c(Karen:::eval.constraints(theta = head(phi.curr,-2), constr.lst = res.fit$constr.lst), phi.curr)

  adjNet <- matrix(data = NA, nrow = nrow(res.fit$V), ncol = nrow(res.fit$V))
  rownames(adjNet) <- colnames(adjNet) <- rownames(res.fit$V)
  for (j in 1:length(head(phi.curr,-2))) {
    r <- as.vector(unlist(strsplit(names(phi.curr)[j], split = "->", fixed = TRUE)))[1]
    p <- as.vector(unlist(strsplit(names(phi.curr)[j], split = "->", fixed = TRUE)))[2]
    if(p != "0" & p != "1"){
      adjNet[p,r] <- phi.curr[j]
    }else{
      if(p == "1"){
        adjNet[r,r] <- phi.curr[j]
      }
    }
  }

  adjNet[adjNet == 0] <- 1e-8

  adjNet <- t(t(adjNet)/colSums(adjNet, na.rm = T))
  adjNet[is.nan(adjNet)] <- 0
  adjNet[is.na(adjNet)] <- 0
  net <- graph_from_adjacency_matrix(t(adjNet), mode = "directed", weighted = TRUE)

  E(net)$width <- 15
  V(net)[names(cols)]$color <- cols

  edge.start <- ends(net, es=E(net), names=F)[,1]
  edge.loops <- which(apply(ends(net, es=E(net), names=F), 1, function(e){e[1]==e[2]}))
  V(net)$label.cex <- c(1.8, rep(1.8, nrow(res.fit$V) - 1))
  arrow.width <- rep(2.5, length(E(net)))
  arrow.width[edge.loops] <- 1
  E(net)[edge.loops]$width <- 7

  col.pal <- colorRampPalette(c("lightgray", "red", "black"))(100)
  rnks <- round(head(rank(c(E(net)$weight, seq(0,1,length.out = 100 - length(E(net)$weight)))), length(E(net))))
  E(net)[order(E(net)$weight, decreasing = FALSE)]$color <- col.pal[sort(rnks)]

  latSts.lst <- setdiff(res.fit$latSts.lst, "HSC")

  yCords <- seq(-1,1,length.out = 5)

  if(modelStructure == 1){
    coords <- rbind(c(0, yCords[5]), # HSC
                    c(0,  yCords[4]), # MPP
                    c(0,  yCords[3]), # CMLP
                    cbind(seq(-17, -7, length.out = 2),  yCords[2]), # MEP, GMP
                    cbind(seq(-20, 0, length.out = 4),  yCords[1]), # P, ERY, G, M
                    cbind(seq(3, 17, length.out = 3),  yCords[2])) # NK, B, T
  }
  if(modelStructure == 2 | modelStructure == 3){
    coords <- rbind(c(0,.3), # HSC
                    c(0, .25), # MPP
                    cbind(seq(-2, 2, length.out = 2), .2), # CMP and CLP
                    cbind(seq(-4, -1, length.out = 2), .15), # MEP and GMP
                    cbind(seq(-5, 1, length.out = 4), .1), # P, ERY, G, M
                    cbind(seq(2, 6, length.out = 3), .15)) # NK, B, T
  }
  if(modelStructure == 4){
    coords <- rbind(c(0,.3), # HSC
                    c(0, .25), # MPP
                    cbind(seq(-7, 7, length.out = 3), .2), # CMP, NKP and CLP
                    cbind(seq(-11, -4, length.out = 2), .15), # MEP and GMP
                    cbind(seq(-13, -1, length.out = 4), .1), # P, ERY, G, M
                    cbind(c(3, seq(7, 11, length.out = 2)), .15)) # NK, B, T
  }
  coords <- tail(coords, length(rownames(res.fit$V)))

  Karen:::plot.network(net,
                       edge.color = E(net)$color,
                       vertex.size = c(30, rep(30, nrow(res.fit$V) - 1)),
                       edge.label = ifelse(rep(edges.lab, length(E(net))), round(E(net)$weight,2), ""),
                       edge.label.cex = 2,
                       edge.label.font = 2,
                       layout = coords,
                       vertex.label.family = "Helvetica",
                       edge.label.color = "black",
                       vertex.frame.color = "black",
                       vertex.frame.width = 4,
                       vertex.frame.cex = 2,
                       vertex.label.font = 2,
                       vertex.label.color = "white",
                       vertex.shape="fcircle",
                       edge.arrow.size = arrow.width
  )
  if(AIC){
    mtext(side = 3, text = paste("AIC = ", round(2*res.fit$fit$value + 2*length(res.fit$fit$par),3), sep = ""), cex = 2, font = 2)
  }
}

get.sMoments <- function(res.fit, X = NULL, cell.cols = NULL){
  V <- res.fit$V # net-effect matrix
  nProc <- length(res.fit$cloneChunks) # number of cores
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

      matplot(as.numeric(rownames(Y_NA)), Y_NA[,,res.fit$cloneChunks[[cnk]][cl]], lty = 1, pch = 20, type = 'p', add = F, col = alpha(cols, alpha = .8), cex = 2,
              cex.axis = 2, cex.lab = 2, xlab = "t", ylab = expression("Y"[t]), main = paste("clone ", cl, sep = ""), cex.main = 2,
              xlim = c(0, max(as.numeric(rownames(Y)))), ylim = range(c(Y_NA[,,res.fit$cloneChunks[[cnk]][cl]], mean_smooth, mean_smooth - 1.96*sd_smooth, mean_smooth + 1.96*sd_smooth), na.rm = T))
      if(!is.null(X)){
        matplot(as.numeric(rownames(X)), X[,,res.fit$cloneChunks[[cnk]][cl]], add = T, pch = 1, cex = 1.5, lwd = 2, col = alpha(cols, alpha = .8))
      }
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth, lwd = 2, lty = 1, type = 'l', add = TRUE, col = cols)
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth - 1.96*sd_smooth, lwd = 2, lty = 3, type = 'l', add = TRUE, col = cols)
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth + 1.96*sd_smooth, lwd = 2, lty = 3, type = 'l', add = TRUE, col = cols)
    })
  })
  plot.new()
  legend(x = "center", legend = rownames(V), col = cols, pch = 20, lwd = 5, cex = 1.5)
}
assignInNamespace("get.sMoments",get.sMoments,ns="Karen")
