# Fit PIM
# 
# @param PIMlist input from createPIMList() function 
# @param link a character vector with a single value that determines the used link function. 
# Possible values are "logit", "probit" and "identity". The default is "logit".
# @param n.cores number of cores for parallel computing
# @param verbose a logical value whether to show the progress status 
# @param coxph.aprox a logical value. If TRUE, coxph algorithm will be used to approximate the PIM
# estimates for faster computation
# @param ... additional arguiments to be passed
# 
# @return a list object
# 
# @examples
# #example
# 
# 
# @import pim
# @importFrom parallel detectCores parLapplyLB makeCluster stopCluster 
# 
fit.PIM <- function(PIMlist, link, BPPARAM, coxph.aprox, verbose,  ...){
  if(verbose) message("fitting PIM ....")
  
  expres.mat    <- PIMlist$expression.mat
  model.formula <- PIMlist$model.formula
  X.mat         <- PIMlist$X.mat  
  U.mat         <- PIMlist$U.mat
  
  #LS <- colSums(expres.mat)
  
  if(coxph.aprox){
    coxphApprox(expres.mat=expres.mat, model.formula=model.formula,
                X.mat=X.mat, U.mat=U.mat, BPPARAM=BPPARAM, ...)
  }
  else{
    # if(n.cores == "available"){
    #   cl <- parallel::detectCores()-1
    #   message(".... parallel computing on all available cores (minus 1)")
    #   cl <- parallel::makeCluster(cl)
    # }
    # else if(n.cores>1 & n.cores != "available"){
    #   cl <- n.cores
    #   message(".... parallel computing on ", n.cores, " cores")
    #   cl <- parallel::makeCluster(cl)
    # }
    # else if(n.cores==1){
    #   if(ncol(expres.mat)>200)
    #   {message(".... running PIM on a single core: parallel computing can be used for faster computation")}
    #   cl <- n.cores
    #   cl <- parallel::makeCluster(cl)
    
    
    #parallel::clusterExport(cl, c("expres.mat", "X.mat", "U.mat"), envir = environment())
    #parallel::clusterEvalQ(cl, c(library('pim')))
    #pim.models.list <- parallel::parLapplyLB(cl, seq_len(nrow(expres.mat)), function(tag){
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
      design      <- pim::new.pim.formula(as.formula(paste("y", model.formula)), data = d) # design
      tag.pim.env <- pim::new.pim.env(d, ...) # Creat PIM enviroment 
      Y.pobs      <- pim::response(design) # Calculate pseudo observations  
      
      # model.tag   <- tryCatch(pim.fit(x=pim::model.matrix(design), y=Y.pobs, 
      #                          penv = tag.pim.env, ...), silent=TRUE)
      model.tag   <- try(suppressWarnings(pim.fit(x=pim::model.matrix(design), y=Y.pobs, 
                                      penv = tag.pim.env, link = link)), 
                              #error = function(e){"error"},
                        silent = TRUE )
      
      if(class(model.tag)[[1]] != "try-error"){
        
        b=coef(model.tag)
        names(b)      <- colnames(pim::vcov(model.tag))
        augmented.MPI <- augmentedPI(pim.fit=model.tag, y=y, a=X.mat[,1], x=U.mat, 
                                     link=link, y.pobs=NULL, ...)
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
        }else{
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
    BPPARAM = BPPARAM)
    #parallel::stopCluster(cl)
    pim.models.list
  } 
}
