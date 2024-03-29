
##########################
## SUPPORTING FUNCTIONS ##
##########################

get.modelList <- function(){
  model.lst <- list()

  model.lst[[1]] <- list(rct.lst = c("HSC->P1",
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

  model.lst[[2]] <- list(rct.lst = c("HSC->P1",
                                     "HSC->P2",
                                     "HSC->P3",
                                     "P1->T",
                                     "P1->B",
                                     "P2->NK",
                                     "P2->G",
                                     "P2->M",
                                     "P3->B",
                                     "P3->NK",
                                     "T->0",
                                     "B->0",
                                     "NK->0",
                                     "G->0",
                                     "M->0",
                                     "HSC->1",
                                     "P1->1",
                                     "P2->1",
                                     "P3->1"),
                         constr.lst = c("theta\\[\\'HSC->P1\\'\\]=(theta\\[\\'P1->T\\'\\] + theta\\[\\'P1->B\\'\\])",
                                        "theta\\[\\'HSC->P2\\'\\]=(theta\\[\\'P2->G\\'\\] + theta\\[\\'P2->M\\'\\] + theta\\[\\'P2->NK\\'\\])",
                                        "theta\\[\\'HSC->P3\\'\\]=(theta\\[\\'P3->B\\'\\] + theta\\[\\'P3->NK\\'\\])"),
                         latSts.lst = c("HSC", "P1", "P2", "P3"))

  return(model.lst)
}
