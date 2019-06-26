# Augmented probablistic index
# 
# @description If the model includes nuisance factors that are not of primary interest (for example 
# library size, batch, cell cycle ets), then the estimated probablistic index will be marginalized 
# over these factors to increase efficiency and statistical power. At the current version, it is 
# implemeted only for a two group comparison study. See Vermeulen et al 2015 for further.
# 
# @param pim.fit PIM output for a particular gene
# @param y a vector of the response variable 
# @param a grouping factor 
# @param x a matrix of covariates. A variable is presented in a column
# @param link a character vector with a single value that determines the used link function. 
# Possible values are "logit", "probit" and "identity". The default is "logit".
# @param y.pobs pseudo observations from a pim result
# @param alpha siginificance level
# 
# @return a list object
# 
# @references Vermeulen, K., Thas, O., & Vansteelandt, S. (2015). Increasing the power of the Mann‐Whitney test in randomized experiments through flexible covariate adjustment. \emph{Statistics in medicine}, 34(6), 1012-1030.
# @examples
# #example
#  
# @keywords internal
augmentedPI <- function(pim.fit, y, a, x,link="logit", y.pobs=NULL, alpha=0.05){
  
  a <- as.factor(a)
  if(length(unique(a))!=2){
    #message("Augmented marginal PI is not applicable for a grouping factor with level morethan 2.")
    return(c(Augmented.MPI=NA, 
             Std.Error=NA, 
             Wald.statistic=NA, 
             p.value=NA))
  }
  else if(is.null(x)){
   # message("Augmented marginal PI is not applicable onl if there is additional covariate.")
    return(c(Augmented.MPI=NA, 
             Std.Error=NA, 
             Wald.statistic=NA, 
             p.value=NA))
  }
  else{
    a <- ifelse(a==levels(a)[1], 0, 1)
    n <- length(a)
    p <- sum(a)/n
    
    
    x.dat    <- x
    if(is.null(colnames(x.dat))) {colnames(x.dat) <- paste0("U", 1:ncol(x))}
    
    x <- as.matrix(do.call(cbind, lapply(x, function(xx){
      if(is.factor(xx)) {
        as.numeric(xx)}
      ##as.numeric(ifelse(xx==levels(xx)[1], 0, 1))}
      else if(is.character(xx)){
        xx=as.factor(xx)
        as.numeric(xx)
        #as.numeric(ifelse(xx==levels(xx)[1], 0, 1))
        #message("One of the character variables in the U matrix is converted to a factor and a dummy variable is created")
      }
      else if(is.numeric(xx)){as.numeric(xx)}
      else{stop("Unable to identify the type of one of the variables in the U matrix")}
    })))
    #x <- as.matrix(x)
    
    if(is.null(colnames(x))) {colnames(x) <- colnames(x.dat)}
    data.sub <- as.data.frame(cbind(y,a,x)) 
    #data.sub <- as.data.frame(do.call(data.frame, data.sub))
    
    # Point estimate:
    # pim.fit <- pim::pim(formula=as.formula(paste("y~a+",paste(colnames(x.dat), collapse="+"))),
    #                     link=link, data=data.sub)
    coef <- pim::coef(pim.fit) 
    MW   <- 1-wilcox.test(y~a,exact=FALSE)$statistic/(sum(a)*sum(1-a))
    
    if(link=="probit"){
      augMW<-MW
      for(i in 1:n){
        augMW<-augMW+sum((1/(n*(n-1))-(1-a[i])*a[-i]/
                            (sum(a)*(n-sum(a))))*pnorm(coef[1]+t(-x[i,]+t(x[-i,]))
                                                       %*%coef[2:(dim(x)[2]+1)]))}
    }
    else if(link=="logit"){
      augMW<-0
      for(i in 1:n){
        augMW<-augMW+sum(expit(coef[1]+t(-x[i,]+t(x[-i,])) %*%coef[2:(dim(x)[2]+1)]))/(n*(n-1))
      }
    }
    else{stop("not a valid link function")}
    
    # Standard error:
    pseudo.y <- pseudo(y) 
    a1.hat<-sapply(1:length(a), A=a, pseudo.Y=pseudo.y, FUN=vec.a1hat)
    a2.hat<-sapply(1:length(a), A=a, pseudo.Y=pseudo.y, FUN=vec.a2hat)
    phi0<-vec.phi0(a,a1.hat,a2.hat,augMW)
    
    if(link=="probit"){
      pred<-t(sapply(coef=coef,x=x,1:dim(x)[1],vec.pred.probit))
      alphahat<-sapply(A=a,pred=pred,1:length(a),vec.alphahat)
      phiest<-phi0+(a-p)*(alphahat-mean(alphahat)+mean((1-a)*a1.hat/(1-p)^2-a*a2.hat/p^2))
    }
    else if(link=="logit"){
      pred<-t(sapply(1:dim(x)[1], coef=coef,x=x,vec.pred.logit))
      alphahat<-sapply(1:length(a),A=a,pred=pred,vec.alphahat)
      phiest<-phi0+(a-p)*alphahat
    }
    else{print("not a valid link function")}
    
    
    se.augMW<-sqrt(mean(phiest^2)/n) 
    W <-(augMW-0.5)/se.augMW # Wald test statistic:
    p.value<-2*pnorm(abs(W),lower.tail=FALSE) # p-value Wald test:
    # CI<-augMW+c(-1,1)*se.augMW*qnorm(1-alpha/2)
    
    
    return(c(Augmented.MPI=augMW, 
             Std.Error=se.augMW, 
             Wald.statistic=W, 
             p.value=p.value))
  }  
  
}


