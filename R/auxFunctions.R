# Auxiliary functions 1: for augmented PI
# 
# @param x a numerical value
# 
# @return a numerical value
#  
#  
# @keywords internal
expit<-function(x){exp(x)/(1+exp(x))}



# Auxiliary functions 2: for augmented PI
# 
# @param Y a numerical value
# @keywords internal
pseudo<-function(Y){
  I<-matrix(rep(Y,length(Y)),ncol=length(Y),byrow=TRUE)
  I<-ifelse(I<Y,1,ifelse(I==Y,0.5,0))
  return(t(I))
}



# Auxiliary functions 3: for augmented PI
# 
# @param A a vector of covariates
# @param pseudo.Y a matrix of pseudo obs
# @param i indices
# @keywords internal

vec.a1hat<-function(A,pseudo.Y,i){
  sum(A[-i]*pseudo.Y[i,-i])/((sum(A)/length(A))*(length(A)-1))
}



# Auxiliary functions 4: for augmented PI
# 
# @param A a vector of covariates
# @param pseudo.Y a matrix of pseudo obs
# @param i indices
# @keywords internal
#
vec.a2hat<-function(A, pseudo.Y, i){
  sum((1-A[-i])*pseudo.Y[-i,i])/((sum(1-A)/length(A))*(length(A)-1))
}



# Auxiliary functions 5: for augmented PI
# 
# @param a1.hat a1.hat
# @param a2.hat a2.hat
# @param est est
# @keywords internal

vec.phi0<-function(A,a1.hat,a2.hat,est){
  p<-sum(A)/length(A)
  (1-A)/(1-p)*a1.hat+A/p*a2.hat-2*est
}


# Auxiliary functions 6: for augmented PI
# 
# @param coef coef
# @keywords internal

vec.pred.probit<-function(coef,x,i){
  pnorm(coef[1]+t(-x[i,]+t(x[,]))%*%coef[2:(dim(x)[2]+1)])
}


# Auxiliary functions 7: for augmented PI
# @keywords internal

vec.pred.logit<-function(coef,x,i){
  expit(coef[1]+t(-x[i,]+t(x[,]))%*%coef[2:(dim(x)[2]+1)])
}


# Auxiliary functions 8: for augmented PI
# 

# @param pred ????
# @keywords internal

vec.alphahat<-function(A,pred,i){
  p<-sum(A)/length(A)
  sum(pred[i,-i]/(1-p)-pred[-i,i]/p)/(length(A)-1)
}


# Auxiliary functions 5: Added function to creat dummay variable for factor variables in X
# 
# @param ref reference group
# @importFrom dummies dummy
# @keywords internal
creatDummay <- function(x, ref=1){
  x.d <- dummy(x)[,-ref]
  x.d
}








