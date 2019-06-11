#' Test contrast from PIMSeq results
#' 
#' @description Tests contrasts using results obtained from PIMSeq function. This function
#' can be called independently after fitting PIM for gene expression data.
#' 
#' @param pim.res PIMSeq output 
#' @param contrasts a vector of coefficients for the linear combination of the factors 
#' @param only.conditions a logical value If only.conditions=TRUE (default), contrast will
#' be tested among the coefficients of condition variables only. If only.conditions=FALSE,
#' contrast will be tested considering all coefiecients in the model (conditions and
#' nuisance variables). The dimension of the contrast matrix/vector depends on this argument. 
#' @param gene.warning.return a logical value to show warning message if the contrast 
#' test for a particular gene encounters an error. Default is FALSE.
#' @param ... additional arguiments
#' 
#' 
#' @return a data frame containg the estimated linear factor, test statistic, standard 
#' error, p value, estimated PI, and FDR adjusted p value.
#' 
#' @examples
#'  # examples
#' @export


testPIMcontrast<- function(pim.res,  contrasts, only.conditions=TRUE, 
                           gene.warning.return=FALSE, ...){
  if(is.null(pim.res)) stop("PIM results not found.") 
  if(is.null(contrasts)) stop("Must specify contrasts.")
  if(anyNA(contrasts))  stop("NAs not allowed in contrasts")
  
  contrasts <- as.matrix(contrasts) 
  if(only.conditions){
    name.X <- pim.res$add.PIM.results$PIMSeq.inputs$condition
  }
  else{
    name.X <- c(pim.res$add.PIM.results$PIMSeq.inputs$condition,
                pim.res$add.PIM.results$PIMSeq.inputs$nuisance_vars)
  }
  
  per.gene.test <- lapply(pim.res$add.PIM.results$fit.model, function(pim.res.gene){
    sub.coef <- sort(unique(do.call("c", lapply(name.X, 
                                           function(nm) grep(nm, names(pim.res.gene$b)))))) 
    b    <- as.matrix(pim.res.gene$b[sub.coef])
    v    <- pim.res.gene$v[sub.coef,sub.coef]
    if(anyNA(b) | anyNA(v) | anyNA(v<0)){
      if(gene.warning.return) warning("Either the coefficients or the var-covar matrix contain NA.")
      
      lf   <- NA
      v.lf <- NA
      z.lf <- NA
      
      p.val <- NA
      PI.DE <- NA 
    }
    else{
      if(nrow(contrasts) != length(b)) stop("Number of columns of contrast matrix must match number of coefficients in fit")
      
      lf   <- t(contrasts)%*%b
      v.lf <- as.numeric(t(as.matrix(contrasts)) %*% v %*% as.matrix(contrasts))
      z.lf <- as.numeric(lf/sqrt(v.lf))
      
      p.val <- as.numeric(2*(1-pnorm(abs(z.lf))))
      PI.DE <- as.numeric(exp(lf)/(1+exp(lf)))
    }
     
    list(ID = as.character(pim.res.gene$tag.name), contrast=as.numeric(lf), 
         Test.Stat=z.lf,  std.error=v.lf, p.value=p.val, PI=PI.DE)
  })
  
  test.contrasts <- as.data.frame(do.call('rbind', per.gene.test))
  test.contrasts$ID         <- as.character(test.contrasts$ID) 
  test.contrasts$contrast   <- as.numeric(test.contrasts$contrast) 
  test.contrasts$Test.Stat  <- as.numeric(test.contrasts$Test.Stat)
  test.contrasts$std.error  <- as.numeric(test.contrasts$std.error)
  test.contrasts$p.value    <- as.numeric(test.contrasts$p.value) 
  test.contrasts$p.adjusted <- p.adjust(test.contrasts$p.value, method="BH")
  test.contrasts$PI         <- as.numeric(test.contrasts$PI)
  
  test.contrasts
} 




























