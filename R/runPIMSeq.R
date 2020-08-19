#' PIM for testing differential gene expression in single cell RNA-seq data.
#' 
#' @description Uses the \emph{pim} package for testing differential gene expression in bulk (with moderate to higher 
#' number of samples) and single cell RNA-seq data.
#' 
#' 
#' @param SCExp a SingleCellExperiment object containing gene expression data, gene and 
#' cell level annotations
#' @param condition a character name(s) that refer to column(s) in colData(SCExp) 
#' indicating the primary factor(s). (see details)
#' @param nuisance.vars a vector of characters that refer to column(s) in colData(SCExp) 
#' indicating nuisance/technical factors such as normalization variables.
#' If nuisance_vars=NULL, PIMSeq assumes the data in SCExp is a normalized data. (see details)
#' @param assay.name a character value indicating the name of the assay (gene expression matrix)
#' in SCExp. Default is counts. See Example 2 for its usage.
#' @param link a character value for the type of link function. Possible values are 
#' "logit", "probit" and "identity". The default is "logit".
#' @param coxph.aprox a logical value. If TRUE, coxph algorithm will be used to approximate the PIM
#' parameters for faster computation
#' @param BPPARAM BiocParallelParam for parallel computing (see BiocParallel). Default is
#' BiocParallel::SerialParam() -- no parallelization
#' @param verbose a logical value whether to show the progress status
#' @param ... additional arguments to be passed to pim such as compare, model, ...
#' 
#' @details This method tests for differential gene expression (DGE) between pre-defined groups
#' such as treatment and cell types. The name of the primary factror(s) can be passed using 
#' the \emph{condition} argument. The primary factor can be a grouping factor or continuous measurement. 
#' In addition, more than one primary factors can be used for testing their main effect and/or 
#' interaction effect. For testing a particular set of contrasts, you can use the \emph{contrast}
#' argument, which is a vector of the contrast/linear factor coefficients (see example). 
#' If \emph{contrast=NULL}, the any difference among the groups will be tested.
#' 
#' PIMSeq uses a semi-parametric regression framework called probabilistic index models (PIM) 
#' which allows it to test DGE by accounting for other sources of unwanted variations, for 
#' example library size, batch, cell cycle, etc. The name of such factor(s) can be passed 
#' using the \emph{nuisance_vars} argument.
#' 
#' PIMSeq normalizes library size difference between samples/cells by incorporating the log
#' library size into the PIM as a nuisance factor. However, it is also possible to use pre-normalized
#' data (normalized with different approaches).
#' 
#' Although PIM is suggested for single cell RNA sequencing data, it can also be used for bulk 
#' RNA-seq data if the number of samples is at least 20.
#' 
#' For large datasets (for example more than 500 cells/samples), PIMSeq can be quite slow. To 
#' speed up the computation, it is possible to use parallel computation by setting the number
#' of cores to be used using \emph{n.cores} argument. If \emph{n.cores}="available", then PIMseq
#' will detect the number of cores in the machine and use all of them except 1. 
#' 
#' PIMSeq calls the \emph{pim.fit()} function in \emph{pim} pckage to fit PIM, and it is possible
#' to use the arguments in the \emph{pim.fit()} by declaring them directly in \emph{PIMSeq()} function.
#' 
#' 
#' @return a list object containing four items: 
#' \describe{
#' \item{test.contrasts}{(data frame) result from testing DGE for each gene using the standard
#' PIM. It includes the test statistic, raw p values, FDR adjusted p values, and the estimated
#' probabilistic index (PI).}
#' \item{all.coefficients}{(data frame) contains the estimated coefficients of the PIM and their standard
#' error. This result is optional.}
#' \item{augmented.MP}{(data frame) contains results from testing DGE using augmented wilcoxon-rank sume
#' test. It contains the marginal estimate of PI and its standard error along with raw p values and
#' FDR adjusted p values.}
#' \item{add.PIM.results}{(list) additional results to be used for further 
#' analysis such as testing contrasts.}
#' }
#' 
#' @references 
#' \itemize{
#'   \item Thas, O., Neve, J. D., Clement, L., & Ottoy, J. P. (2012). Probabilistic index models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 74(4), 623-671.
#'   \item Vermeulen, K., Thas, O., & Vansteelandt, S. (2015). Increasing the power of the Mann‐Whitney test in randomized experiments through flexible covariate adjustment. \emph{Statistics in medicine}, 34(6), 1012-1030.
#' }
#' 
#' @examples  
#' \donttest{
#' library(PIMseq)
#' library(SingleCellExperiment) 
#' 
#' # -----------------------------------------------------------------------------
#' # Example 1
#' # Analysis of Neuroblastoma cell line single cell RNA-seq data for testing DE between
#' # nutlin treated and control (vehicle) cells
#' 
#' data("scNGP.data")
#' dim(scNGP.data)
#' scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data)>0)>10, ]
#' dim(scNGP.data2)
#' 
#' colData(scNGP.data2)$treatment <- 
#' ifelse(colData(scNGP.data2)$characteristics..treatment=="nutlin",1,0)
#' table(colData(scNGP.data2)$treatment)
#' colData(scNGP.data2)$logLS <- log(colSums(counts(scNGP.data2)))
#' 
#' res.PIM <- PIMSeq(SCExp = scNGP.data2, condition = "treatment", 
#'                               nuisance.vars = "logLS", BPPARAM=BiocParallel::SerialParam())
#' 
#' head(res.PIM$test.contrasts)  # result based on the standard PIM
#' head(res.PIM$all.coefficients) # (optional) estimated coefficients and their standard error
#' head(res.PIM$augmented.MP) # (optional) DGE test using augmented PI (if comparing only 2 groups)
#' 
#' #number of detected DE genes at 5% FDR
#' table(res.PIM$test.contrasts$p.adjusted <0.05)
#' 
#' # number of detected DE genes at 5% FDR and PI <0.3 or PI >0.7
#' table(res.PIM$test.contrasts$p.adjusted <0.05 & 
#'         (res.PIM$test.contrasts$PI<0.3 | res.PIM$test.contrasts$PI>0.7))
#' 
#' # number of down-regulated DE genes at 5% FDR and PI <0.3
#' table(res.PIM$test.contrasts$p.adjusted <0.05 & res.PIM$test.contrasts$PI<0.3)
#' 
#' # number of up-regulated DE genes at 5% FDR and PI > 0.7
#' table(res.PIM$test.contrasts$p.adjusted <0.05 & res.PIM$test.contrasts$PI>0.7)
#' 
#' # rank genes based on their estimated PI  
#' plotPIrank(res.PIM, contrast = 1, PI.thrld = c(0.3, 0.7))
#' title(main="Gene ranking based on PI")
#' 
#' # alternative plots to volcano plot
#' plot(res.PIM$test.contrasts$PI, -log10(res.PIM$test.contrasts$p.value), 
#'      col=factor(res.PIM$test.contrasts$p.adjusted<0.05), las=1,
#'      xlab="PI", ylab="-log10(p value)")
#' legend("bottomright", c("FDR<0.05", "FDR>=0.05"), col=c(2,1), pch=19)
#' 
#' 
#' 
#' 
#' 
#' # ------------------------------------------------------------------------------------
#' # Example 2
#' # Analysis of NORMALIZED Neuroblastoma cell line single cell RNA-seq data for 
#' # testing DE between nutlin treated and control (vehicle) cells
#' 
#' # We paricularly use the counts per millions of reads (CPM) to normalize the data.
#' 
#' data("scNGP.data")
#' dim(scNGP.data)
#' scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data)>0)>10, ] 
#' 
#' log2.CPM <- log2(counts(scNGP.data2) %*% diag(1e6/colSums(counts(scNGP.data2))) + 1)
#' names(assays(scNGP.data2))
#' assays(scNGP.data2)$log2CPM <- log2.CPM
#' names(assays(scNGP.data2))
#' 
#' colData(scNGP.data2)$treatment <- 
#'  ifelse(colData(scNGP.data2)$characteristics..treatment=="nutlin",1,0)
#' table(colData(scNGP.data2)$treatment)
#' 
#' res.PIM <- PIMSeq(SCExp = scNGP.data2, condition = "treatment",
#'                               assay.name = "log2CPM", nuisance.vars = NULL)
#' 
#' 
#' head(res.PIM$test.contrasts)  # result based on the standard PIM
#' 
#' # number of detected DE genes at 5% FDR
#' table(res.PIM$test.contrasts$p.adjusted <0.05) 
#' 
#' 
#' 
#' 
#' 
#' # -----------------------------------------------------------------------------
#' # Example 3
#' # Applying PIM for complex experimental design and testing contrasts.
#' # For this purpose, let us simulate artificial data for a 2x2 factorial design.
#' # i.e. two factors each with two levels
#' example.data <- make.example.data(n.gene = 1000, n.sample = 50, n.group = 2, 
#'                                   n.batch = 2, p.DEgenes = 0.2)
#' table(example.data$Group, example.data$Batch)
#' 
#' colData(example.data)$X1 <- as.factor(example.data$Group)
#' colData(example.data)$X2 <- as.factor(example.data$Batch)
#' colData(example.data)$logLS <- log(example.data$E.LS)
#' 
#' keep <- rowSums(counts(example.data)>0)>10
#' table(keep)
#' example.data <- example.data[keep, ]  # filter genes
#' dim(example.data)
#' 
#' 
#' res.PIM2 <- PIMSeq(SCExp = example.data, condition = c("X1", "X2", "X1:X2"), 
#'                    nuisance.vars = "logLS")
#' head(res.PIM2$test.contrasts)
#' table(res.PIM2$test.contrasts$p.adjusted<0.05)  # number of DE genes, but b/n which group
#' 
#' #DE between levels of factor 1 when the level of factor 2 is 1
#' contrast1 <- testPIMcontrast(res.PIM2, contrasts = c(1, 0, 0)) 
#' 
#' #DE between levels of factor 1 when the level of factor 2 is 2
#' contrast2 <- testPIMcontrast(res.PIM2, contrasts = c(1, 0, 1))
#' 
#' #DE between levels of factor 2 when the level of factor 1 is 1 
#' contrast3 <- testPIMcontrast(res.PIM2, contrasts = c(0, 1, 0)) 
#' 
#' #DE between levels of factor 2 when the level of factor 1 is 2
#' contrast4 <- testPIMcontrast(res.PIM2, contrasts = c(0, 1, 1))
#'  
#' table(contrast1$p.adjusted<0.05)
#' table(contrast2$p.adjusted<0.05)
#' table(contrast3$p.adjusted<0.05)
#' table(contrast4$p.adjusted<0.05)
#' }
#' 
#' 
#' @export  
#' 
#' @import  pim 
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' @importFrom survival coxph Surv 
#' @importFrom BiocParallel bplapply SerialParam 
#' @importFrom stats pnorm pchisq p.adjust as.formula rgamma rlnorm rnbinom wilcox.test model.frame
#' @importFrom SummarizedExperiment assays
PIMSeq <- function(SCExp,
                   condition, 
                   nuisance.vars=NULL, 
                   assay.name ="counts", 
                   link="logit", 
                   verbose=TRUE, 
                   coxph.aprox=FALSE, 
                   BPPARAM=BiocParallel::SerialParam(), ...)
  {
  
  #Creat PIMlist object
  PIMlist <- createPIMList(SCExp=SCExp, 
                           name.X=condition, 
                           name.U=nuisance.vars,  
                           assay.name=assay.name,
                           verbose=verbose, ...)
  
  #Fit PIM  
  fit.model <- fit.PIM(PIMlist, link=link, coxph.aprox=coxph.aprox,
                       BPPARAM=BPPARAM, verbose=verbose, ...)   
  
  #Global DE test 
  glob.DE.test <- golbalDEtest(fit.model, SCExp=SCExp,  condition=condition, 
                               nuisance.vars=nuisance.vars, ...)
  return(glob.DE.test)
}
