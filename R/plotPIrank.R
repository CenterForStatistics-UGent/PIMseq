#' PI rank plot
#' 
#' @description Plots gene ranking based on the estimated PI. Low-regulated genes (PI <0.5) appear
#' in the left edge and up-regulated genes appear in the right edge. Error bars indicate the (1-alpha)*100%
#' confidence interval estimate of the PI.
#' 
#' @param pim.res outputs from \emph{PIMSeq()} function
#' @param contrast a vector of coefficients for the linear combination of the factors 
#' @param add.PI.interval a logical value to add error bar
#' @param PI.thrld a vector of numerical values for the down and up regulated genes PI threshold
#' @param FDR.thrld FDR threshold, a numerical value between 0 and 1
#' @param xlab xlab parameter for \emph{plot()}
#' @param ylab ylab parameter for \emph{plot()}
#' @param add.legend a logical value to add a legend
#' @param cols.seg a vector of colors for the error bars. Size must be equal to the number of genes
#' @param las las parameter for \emph{plot()}
#' @param ylim ylim parameter for \emph{plot()}
#' @param col col parameter for \emph{plot()}
#' @param pch pch parameter for \emph{plot()}
#' @param cex cex parameter for \emph{plot()}
#' @param cex.lab cex.lab parameter for \emph{plot()}
#' @param lty.PI lty parameter for \emph{plot()}, for PI horizontal line(s)
#' @param lwd.PI lwd parameter for \emph{plot()}, for PI horizontal line(s)
#' @param col.PI col parameter for \emph{plot()}, for PI horizontal line(s)
#' @param legend.pos legend position
#' @param  pt.cex.legend pch size for legend 
#' @param  cex.legend size of legend
#' @param ... other graphical parameters for \emph{plot()}
#' 
#' 
#' @return plot
#' 
#' @export
#' @examples
#' # see examples in PIMSeq() function
#' 
#' @importFrom graphics plot abline segments points legend
#' 
plotPIrank <- function(pim.res, contrast, add.PI.interval=TRUE, FDR.thrld=0.05,  
                       PI.thrld=c(0.4, 0.6), xlab="gene rank", ylab="probabilstic index", 
                       add.legend=TRUE, cols.seg=NULL, las=1, cex.lab=1.25, ylim=c(0,1), 
                       col=NULL, pch=20, cex=0.5, lty.PI=3, lwd.PI=1.25, col.PI="blue",
                       legend.pos= "bottomright", pt.cex.legend = 1.25, 
                       cex.legend=1, ...){
  
  expit  <- function(x){exp(x)/(1+exp(x))}
  
  df <- testPIMcontrast(pim.res, contrasts = contrast)
  
  #df <- merge(pim.contrast, pim.res$all.coefficients, by="ID") 
  
  if(add.PI.interval){
    df$lPI <- expit(df$contrast-1.96*df$std.error)
    df$uPI <- expit(df$contrast+1.96*df$std.error)
  }
  df  <- df[order(df$PI), ] 
  
  if(is.null(cols.seg)) cols.seg <- ifelse(df$p.value>=FDR.thrld, "gray", "rosybrown2")
  
  plot(1:nrow(df), df$PI, type="n", xlab=xlab,
       ylab = ylab, las=las, cex.lab=cex.lab, ylim=ylim, ...)
  if(add.PI.interval){
    segments(1:nrow(df), df$lPI, 1:nrow(df), df$uPI, col=cols.seg, ...) 
  }
  
  if(is.null(col)){
    col <- factor(df$p.adjusted<FDR.thrld)
  }
  points(1:nrow(df), df$PI, col=col, pch=pch, cex=cex, ...)
  
  abline(h=PI.thrld, col=col.PI, lty=lty.PI, lwd=lwd.PI) 
  
  if(add.legend){
    legend("bottomright", c(paste("FDR>=", FDR.thrld), paste("FDR<", FDR.thrld)), 
           col = unique(col), pch = pch, cex=cex.legend, pt.cex = pt.cex.legend, ...)
  } 
}