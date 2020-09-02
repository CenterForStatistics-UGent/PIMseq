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
#'  # see examples of PIMSeq function.
#' @export


testPIMcontrast<- function(pim.res,  contrasts, only.conditions=TRUE, 
                           gene.warning.return=FALSE, ...){
  if(is.null(pim.res)) stop("PIM results not found.") 
  if(is.null(contrasts)) stop("Must specify contrasts.")
  if(anyNA(contrasts))  stop("NAs not allowed in contrasts")
  if(!is.matrix(contrasts)) stop("Contrast should be a matrix class. Rows showd be contrasts and the number of columns should be equal to the number of coefficients in the model.")
  
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
    }else{
      if(nrow(contrasts)>1){
        if(ncol(contrasts) != length(b)) stop("Number of columns of contrast matrix must match number of coefficients in fit")
        lf   <- contrasts %*% b
        v.lf <- contrasts %*% v %*% t(contrasts)
        v.lf.inv <- try(solve(v.lf), silent = TRUE)
        PI.DE.all <- as.numeric(exp(lf)/(1+exp(lf))) 
        PI.DE <- PI.DE.all[which.max(abs(PI.DE.all-0.5))]
        if(all(class(v.lf.inv) != "try-error")){
          z.lf  <- as.numeric((t(lf) %*% v.lf.inv %*% lf))
          p.val <- as.numeric(1-pchisq(q=z.lf, df = nrow(contrasts)))
        }else{
          z.lf  <- NA
          p.val <- NA 
        } 
        out.list <- data.frame(ID = as.character(pim.res.gene$tag.name),
                               PI=PI.DE, Chis.quare=z.lf,  p.value=p.val)
       
        
      }else{
        if(ncol(contrasts) != length(b)) stop("Number of columns of contrast matrix must match number of coefficients in fit")
        lf   <- as.numeric(contrasts%*%b)
        v.lf <- as.numeric(contrasts %*% v %*% t(contrasts))
        PI.DE.all <- as.numeric(exp(lf)/(1+exp(lf))) 
        PI.DE <- PI.DE.all[which.max(abs(PI.DE.all-0.5))]
        if(v.lf>0){
          z.lf <- as.numeric(lf/sqrt(v.lf))
          p.val <- as.numeric(2*(1-pnorm(abs(z.lf))))
        }else{
          z.lf  <- NA
          p.val <- NA  
        }
        
        out.list <- data.frame(ID = as.character(pim.res.gene$tag.name), 
                               PI=PI.DE,
                               contrast=as.numeric(lf),
                               std.error=as.numeric(v.lf),
                               Z=z.lf,  p.value=p.val)
      }
    } 
    out.list
  })
  
  test.contrasts <- as.data.frame(do.call('rbind', per.gene.test))
  test.contrasts$p.adjusted <- p.adjust(test.contrasts$p.value, method="BH") 
  
  test.contrasts
} 




























