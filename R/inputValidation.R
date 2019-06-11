# A function that checks the validity of input data
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
# 
# 
# @return a list object
# 
# @examples 
# # example 
# @import SingleCellExperiment
# @importFrom SummarizedExperiment assays
# 
Validate.input <- function(SCExp, name.X, name.U, assay.name){  
  if(is.null(SCExp)) { 
    return(list(valid  = FALSE, 
                message= "Error: SCExp can not be NULL."))
  }
  else if(is.null(name.X)) { 
    return(list(valid  = FALSE, 
                message= "Error: Condition name can not NULL."))
  }
  else if(class(SCExp) != "SingleCellExperiment"){ 
    return(list(valid  = FALSE, 
                message= "Error: Input is not a 'SingleCellExperiment' class")) 
  }
  else if(is.null(assays(SCExp)[[assay.name]])){ 
    return(list(valid  = FALSE, 
                message= "Error: Assay is not found in SCExp")) 
  }
  else if(!is.null(assays(SCExp)[[assay.name]])){ 
    if(anyNA(assays(SCExp)[[assay.name]])){
      return(list(valid  = FALSE, 
                  message= "Error: Assay contains NA.")) 
    }
    else{
      return(list(valid  = TRUE, 
                  message= NULL))
    }
  }
  # else if(!all(name.X %in% names(colData(SCExp)))){ 
  #   return(list(valid  = FALSE, 
  #               message= paste("Error: The following factor(s) are not found in the colData(SCExp).", 
  #                              name.X[which(!(name.X %in% names(colData(SCExp))))])))
  # }
  # else if(!all(name.U %in% names(colData(SCExp)))){ 
  #   return(list(valid  = FALSE, 
  #               message= paste("Error: The following factor(s) are not found in the colData(SCExp)..", 
  #                              name.U[which(!(name.U %in% names(colData(SCExp))))]))) 
  # }
  else{ 
    return(list(valid = TRUE,
                message = NULL))
  }
}