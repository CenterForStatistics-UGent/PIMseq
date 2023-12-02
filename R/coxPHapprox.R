# Fit CoxPH approximation of PIM
# 
# @param expres.mat expression matrix
# @param model.formula model formula
# @param X.mat a matrix of primary factors
# @param U.mat a matrix of nuisance factors
# @param ties a character string specifying the method for tie handling. see coxph() function
# under survivl pckage. We set the default to "exact", because the “exact partial likelihood” is 
# equivalent to a conditional logistic model, and is appropriate when the times are a 
# small set of discrete values, which is the case for RNA-seq data.
# @param ... additional arguiments to be passed to coxph function
# 
# @return a list object
# 
# @examples
# #example
# 
# @export
# @import survival
# 
coxphApprox <- function(expres.mat, model.formula, X.mat, U.mat, 
                        ties="efron", BPPARAM, ...){
  message("Using CoxPH as an approximation to the PIM for faster computation.")
  
  pim.models.list <- BiocParallel::bplapply(seq_len(nrow(expres.mat)), function(tag){
    #pim.models.list <- lapply(seq_len(nrow(expres.mat)), function(tag){ 
    y <- log2(as.vector(expres.mat[tag, ])+1) 
    if(!is.null(U.mat)){
      d <- as.data.frame(cbind.data.frame(y,X.mat,U.mat))
      d <- as.data.frame(do.call(data.frame, d)) 
    }
    else if(is.null(U.mat)){
      d <- as.data.frame(cbind.data.frame(y,X.mat))  
    }
    else{
      stop("Please specifiy at least one covariate!")
    } 
    
    # Fit models 
    Cns <- rep(1, nrow(d))
    fml <- as.formula(paste("survival::Surv(y, Cns)", model.formula))
    # if(tie.adjust){
    #   model.tag <- try(coxph(fml, data = d, ties=ties, ...), silent=TRUE)
    # }else{
    #   model.tag <- try(coxph(fml, data = d, ...), silent=TRUE)
    # }
    model.tag <- try(coxph(fml, data = d, ties=ties, ...), silent=TRUE)
    
    
    if(class(model.tag)[[1]] != "try-error"){
      
      b=-1*coef(model.tag)
      names(b)      <- colnames(pim::vcov(model.tag))
      augmented.MPI <- c(Augmented.MPI=NA,  Std.Error=NA,  Wald.statistic=NA, p.value=NA)
      tag.pim.res   <- list(b=b, v=pim::vcov(model.tag), tag.name = rownames(expres.mat)[tag], 
                            augmented.MPI=augmented.MPI)
      tag.pim.res
    }
    else{ 
      if(!is.null(U.mat)){
        r <- sum(sapply(1:ncol(X.mat), function(x){
          if(is.factor(X.mat[,x])){length(levels(X.mat[,x]))-1}
          else{1}
        }))+ sum(sapply(1:ncol(U.mat), function(x){
          if(is.factor(U.mat[,x])){length(levels(U.mat[,x]))-1}
          else{1}
        })) 
        
        names <- c(as.character(sapply(1:ncol(X.mat), function(x){
          if(is.factor(X.mat[,x])){paste0(colnames(X.mat)[x],levels(X.mat[,x])[-1])}
          else{colnames(X.mat)[x]}
        })),
        as.character(sapply(1:ncol(U.mat), function(x){
          if(is.factor(U.mat[,x])){paste0(colnames(U.mat)[x],levels(U.mat[,x])[-1])}
          else{colnames(U.mat)[x]}
        })))
      }
      else{
        r <- sum(sapply(1:ncol(X.mat), function(x){
          if(is.factor(X.mat[,x])){length(levels(X.mat[,x]))-1}
          else{1}
        }))
        
        names <- c(as.character(sapply(1:ncol(X.mat), function(x){
          if(is.factor(X.mat[,x])){paste0(colnames(X.mat)[x],levels(X.mat[,x])[-1])}
          else{colnames(X.mat)[x]}
        })))
      }
      
      # main.test <- gdata::unmatrix(main.test,byrow=FALSE)
      # main.test
      
      b=rep(NA, r) ; names(b) <- names
      v=matrix(rep(NA, r*r), r, r) ; rownames(v) <- colnames(v) <- names
      augmented.MPI=c(Augmented.MPI=NA,  Std.Error=NA,  Wald.statistic=NA, p.value=NA)
      tag.pim.res <- list(b=b, v=v, tag.name = rownames(expres.mat)[tag], augmented.MPI=augmented.MPI)
      tag.pim.res
    }
    
  },
  BPPARAM=BPPARAM) 
  pim.models.list  
}