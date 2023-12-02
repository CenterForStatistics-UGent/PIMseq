# internal function for testing global DE
golbalDEtest <- function(fit.model, SCExp, nuisance.vars, condition, link, ...){
  test.contrasts <- as.data.frame(do.call('rbind', 
        lapply(fit.model, function(mod){ 
    sub.coef <- sort(unique(do.call("c", 
              lapply(condition,  function(nm) grep(nm, names(mod$b)))))) 
    
    b    <- as.matrix(mod$b[sub.coef])
    v    <- mod$v[sub.coef, sub.coef]
    if(length(sub.coef)==1){
      z.stat  <- try(b/sqrt(v), silent = TRUE)
      if(class(z.stat)[[1]] == "try-error"){
        z.stat <- 0
      }
      pval <- as.numeric(2*pnorm(abs(z.stat), lower.tail = FALSE)) 
      PI.DE <- exp(b)/(1+exp(b))
      data.frame(ID = as.character(mod$tag.name), 
                 PI=PI.DE, Z=z.stat, p.value=pval)
    }
    else if(length(sub.coef)>1){
      lrt  <- try(as.numeric(t(b) %*% solve(v) %*% b), silent = TRUE)
      if(class(lrt)[[1]] == "try-error"){
        lrt <- 0
      }
      pval <- as.numeric(pchisq(lrt, nrow(b), lower.tail = FALSE))
      max.abs.beta <- b[which.max(abs(b))]
      PI.DE <- exp(max.abs.beta)/(1+exp(max.abs.beta))
      data.frame(ID = as.character(mod$tag.name), 
                 PI=PI.DE,
                 Chi.square=lrt, p.value=pval)
    } 
  }))) 
  test.contrasts$p.adjusted <- p.adjust(test.contrasts$p.value, method="BH") 
  
  all.coefficients <- as.data.frame(t(sapply(fit.model, function(mod){
    b         <- mod$b 
    names(b)  <- paste0("beta:", names(b))
    std.error <- sqrt(diag(mod$v))
    names(std.error) <- paste0("SE:", rownames(mod$v))
    tag <- mod$tag.name
    names(tag) <- "ID"
    
    c(tag, b, std.error)
  })))
  all.coefficients[, -1] <- apply(all.coefficients[,-1], 2, function(x) as.numeric(as.character(x)))
  
  augmented.MP.df <- as.data.frame(t(sapply(fit.model, function(mod){
    mod$augmented.MPI
  })))
  augmented.MP.df$p.adjusted <- p.adjust(augmented.MP.df$p.value, method = "BH")
  augmented.MP.df$ID  <- test.contrasts$ID
  
  list(test.contrasts=test.contrasts, all.coefficients=all.coefficients, 
       augmented.MP=augmented.MP.df, 
       add.PIM.results = list(fit.model=fit.model, 
                              PIMSeq.inputs=list(SCExp=SCExp, condition=condition, 
                                                 nuisance.vars=nuisance.vars, link=link)
                              )
       )
}
