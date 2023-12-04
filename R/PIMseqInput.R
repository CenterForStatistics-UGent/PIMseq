# Prepare inputs for PIM 
# 
# @param SCExp a SingleCellExperiment object contating gene expression data, gene and 
# cell level annotations
# @param name.X a character name(s) that refer to column(s) in colData(SCExp) 
# indicating the main factor(s) across which we test for DE
# @param name.U a vector of characters that refer to column(s) in colData(SCExp) 
# indicating variables to be used for normalization and/or batch correction. If norma_vars=NULL,
# PIMSeq assumes the data in SCExp is a normalized data. 
# @param assay.name a character value indicating the name of the assay (gene expression matrix)
# in SCExp. 
# @param verbose a logical value whether to show the progress status
# @param ... additional arguiments
# @return a list object
# 
# @examples 
#  #example
#  
# @import SingleCellExperiment
# @importFrom SummarizedExperiment assays
createPIMList <- function(SCExp, name.X,  name.U,  assay.name, verbose, ...){
 
  #Validate inputs
  if(verbose) message("input validation ....")
  valid.input <- Validate.input(SCExp, name.X, name.U, assay.name, ...)
  if(!valid.input$valid){stop(valid.input$message)}
  
  if(verbose) message("creating PIMlist ....")
  
  #Creat matrix of factors
  if(!is.null(name.X)){ 
    unique.name.X <- unique(do.call('c', lapply(name.X, function(x) unlist(strsplit(x, "*", fixed = TRUE)))))
    check.X <- sapply(unique.name.X, function(x){
      if(x %in% names(colData(SCExp))) {TRUE}
      else {FALSE}
    })
    if(!all(check.X)) {stop(paste("The condition name(s)", unique.name.X[!check.X], "are not found in colData()."))}
    
    
    X.mat <- as.data.frame(colData(SCExp)[ ,unique.name.X])
    colnames(X.mat) <- unique.name.X
  }
  else{
    stop("Please specify the condition name.")
  }
  
  if(!is.null(name.U)){
    unique.name.U <- unique(do.call('c', lapply(name.U, function(x) unlist(strsplit(x, "*", fixed = TRUE)))))
    check.U <- sapply(unique.name.U, function(u){
      if(u %in% names(colData(SCExp))) {TRUE}
      else {FALSE}
    })
    if(!all(check.U)) {stop(paste("The additional covariate name(s)", unique.name.U[!check.U], "are not found in colData()."))}
    
    U.mat  <- as.data.frame(colData(SCExp)[ ,unique.name.U]) 
    colnames(U.mat) <- unique.name.U
  }
  else{
    U.mat <- NULL
  }
  
  
  if(!is.null(name.X) & !is.null(name.U)){
    model.formula <- paste("~", paste(c(name.X, name.U), collapse = "+"))
  }
  else if(!is.null(name.X) & is.null(name.U)){
    model.formula <- paste("~", paste(name.X, collapse = "+"))
  } 
  
    
  expression.mat <- assays(SCExp)[[assay.name]]
   
  #Creat PIMList object
  list(expression.mat = expression.mat,
       model.formula=model.formula,
       X.mat = X.mat,
       U.mat = U.mat)
} 
