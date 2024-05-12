 ### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: May 12 2024 (17:05) 
##           By: Brice Ozenne
##     Update #: 1408
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.lmm (documentation)
##' @title Summary Output for a Linear Mixed Model
##' @description Summary output for a linear mixed model fitted with \code{lmm}.
##'
##' @param object [lmm] output of the \code{lmm} function.
##' @param level [numeric,0-1] confidence level for the confidence intervals.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param print [logical] should the output be printed in the console.
##' @param columns [character vector] Columns to be output for the fixed effects.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param digits [interger, >0] number of digits used to display estimates.
##' @param digits.df [interger, >0] number of digits used to display degrees of freedom.
##' @param digits.p.value [interger, >0] number of digits used to display p-values.
##' @param hide.data [logical] should information about the dataset not be printed.
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param hide.cor [logical] should information about the correlation structure not be printed.
##' @param type.cor [character] should the correlation matrix be display (\code{"matrix"}) or the parameter values (\code{"param"}).
##' @param hide.sd [logical] should information about the standard deviation not be printed.
##' @param hide.var [logical] should information about the variance not be printed.
##' @param hide.re [logical] should information about the random effect not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.
##'
##' @return A list containing elements displayed in the summary: \itemize{
##' \item \code{correlation}: the correlation structure.
##' \item \code{variance}: the variance structure.
##' \item \code{sd}: the variance structure expressed in term of standard deviations.
##' \item \code{mean}: the mean structure.
##' }
##'
##' @keywords methods

## * summary.lmm (code)
##' @export
summary.lmm <- function(object, level = 0.95, robust = FALSE,
                        print = TRUE, columns = NULL, digits = 3, digits.df = 1, digits.p.value = 3, 
                        hide.data = FALSE, hide.fit = FALSE, hide.cor = NULL, type.cor = NULL, hide.var = NULL, hide.sd = NULL, hide.re = NULL, hide.mean = FALSE, ...){

    ## ** extract from object
    param.value <- object$param
    param.type <- stats::setNames(object$design$param$type, object$design$param$name)
    
    param.mu <- param.value[names(which(param.type=="mu"))]
    param.sigma <- param.value[names(which(param.type=="sigma"))]
    param.k <- param.value[names(which(param.type=="k"))]
    param.rho <- param.value[names(which(param.type=="rho"))]
    data <- object$data
    call <- object$call
    structure <- object$design$vcov
    structure.ranef <- structure$ranef

    logLik <- stats::logLik(object)
    nobs <- stats::nobs(object)
    method.fit <- object$args$method.fit
    type.information <- object$args$type.information
    nobsByCluster <- lengths(object$design$index.cluster)
    formula <- object$formula
    df <- !is.null(object$df)
    options <- LMMstar.options()

    n.strata <- object$strata$n
    U.strata <- object$strata$levels
    
    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    valid.columns <- c("estimate","se","df","lower","upper","null","statistic","p.value","")
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
        }
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            columns <- union(options$columns.summary, unname(columns))
        }
        if(!is.null(names(columns)) && all(names(columns)=="remove")){
            columns <- setdiff(options$columns.summary, unname(columns))
        }
    }else{
        columns <- options$columns.summary
    }    

    if(!is.null(type.cor)){
        type.cor <- match.arg(type.cor, c("matrix","param"))
    }

    if(inherits(structure,"RE")){
        if(is.null(hide.cor)){hide.cor <- TRUE}
        if(is.null(hide.sd)){hide.sd <- TRUE}
        if(is.null(hide.var)){hide.var <- TRUE}
        if(is.null(hide.re)){hide.re <- FALSE}
    }else{
        if(is.null(hide.cor)){hide.cor <- is.null(object$formula$cor)}
        if(is.null(hide.sd)){hide.sd <- FALSE}
        if(is.null(hide.var)){hide.var <- TRUE}
        if(is.null(hide.re)){hide.re <- TRUE}
    }

    ## ** welcome message
    if(print){
        if(length(param.rho) == 0){
            cat("\t\tLinear regression \n")
        }else{
            cat("\t\tLinear Mixed Model \n")
        }
        cat(" \n")
    }

    
    ## ** data message    
    if(print && !hide.data){
        if(inherits(call$data,"call") || inherits(call$data,"name")){
            cat("Dataset:", deparse(call$data), "\n\n")
        }else{
            cat("Dataset:\n\n")
        }
        if(nobs["missing.obs"]>0){
            if(nobs["missing.cluster"]>0){
                cat("  - ", nobs["cluster"], " clusters were analyzed, ",nobs["missing.cluster"]," were excluded because of missing values \n" , sep = "")
            }else{
                cat("  - ", nobs["cluster"], " clusters \n" , sep = "")
            }
            cat("  - ", nobs["obs"], " observations were analyzed, ",nobs["missing.obs"]," were excluded because of missing values \n",  sep = "")
        }else{
            cat("  - ", nobs["cluster"], " clusters \n" , sep = "")
            cat("  - ", nobs["obs"], " observations \n",  sep = "")
        }
        if(length(unique(nobsByCluster))==1){
            cat("  - ", nobsByCluster[1], " observations per cluster \n", sep = "")
        }else{
            cat("  - between ", min(nobsByCluster), " and ",max(nobsByCluster)," observations per cluster \n", sep = "")
        }

        cat("\nSummary of the outcome and covariates: \n\n")
        data.XY <- data[all.vars(stats::terms(formula$mean))]
        str.XY <- utils::capture.output(utils::str(data.XY))[-1]
        str.XY[1] <- paste0(" ",str.XY[1])
        cat(paste0("  ",str.XY,"\n"))

        reference.level <- levels(object)$reference
        if(!is.null(reference.level)){
            cat("    reference level: ",paste(paste(names(reference.level),reference.level,sep="="),collapse=";")," \n", sep = "")
        }
        cat("\n")
    }
    
    ## ** optim message    
    if(print && !hide.fit){
        cat("Estimation procedure \n\n")
        if(method.fit == "REML"){
            cat("  - Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - Maximum Likelihood (ML) \n")
        }
        cat("  - log-likelihood :", as.double(logLik), "\n",sep="")
        cat("  - parameters: mean = ",length(param.mu),", variance = ",length(c(param.sigma,param.k)),", correlation = ",length(param.rho),"\n", sep = "")
            abs.score <- abs(object$score)
            abs.diff <- abs(object$opt$previous.estimate-object$param)
            name.score <- names(which.max(abs.score))[1]
            name.diff <- names(which.max(abs.diff))[1]
            
        cat("  - convergence: ",object$opt$cv>0," (",object$opt$n.iter," iterations) \n",
            "    largest |score| = ",max(abs.score)," for ",name.score,"\n",
            if(!is.null(name.diff)){paste0("            |change|= ",max(abs.diff)," for ",name.diff,"\n")},
            sep = "")

        cat(" \n")
    }

    ## ** vcov structure
    if(print && (!hide.cor || !hide.var || !hide.sd || !hide.re)){
        cat("Residual variance-covariance: ")
        if(is.na(structure$name$strata)){
            txt.strata <- ""
        }else{
            txt.strata <- "stratified "
        }
        
        if(length(c(param.sigma,param.k))==0){
            hide.sd <- TRUE
            hide.var <- TRUE
        }
        if(inherits(structure,"RE")){
            if(structure.ranef$crossed==FALSE && structure.ranef$nested==FALSE){
                cat("random intercept \n", sep = "")
            }else if(structure.ranef$crossed==FALSE && structure.ranef$nested==TRUE){
                cat("nested random intercepts \n", sep = "")
            }else if(structure.ranef$crossed==TRUE && structure.ranef$nested==FALSE){
                cat("cross random intercepts \n", sep = "")
            }else{
                cat("random effects \n", sep = "")
            }        
        }else if(inherits(structure,"ID")){
            cat(txt.strata,"identity \n\n",sep="")         
        }else if(inherits(structure,"IND")){
            cat(txt.strata,"diagonal \n\n",sep="")
        }else if(inherits(structure,"UN")){
            cat(txt.strata,"unstructured \n\n",sep="")
        }else if(inherits(structure,"CS")){
            if(all(is.na(structure$name$cor))){
                cat(txt.strata,"compound symmetry \n\n",sep="")
            }else if(structure$type == "heterogeneous"){
                cat(txt.strata,"block unstructured \n\n",sep="")
            }else if(structure$type == "homogeneous"){
                cat(txt.strata,"block compound symmetry \n\n",sep="")
            }else if(structure$type == "heterogeneous0"){
                cat(txt.strata,"crossed unstructured \n\n",sep="")
            }else if(structure$type == "homogeneous0"){
                cat(txt.strata,"crossed compound symmetry \n\n",sep="")
            }
        }else if(inherits(structure, "TOEPLITZ")){            
            if(all(is.na(structure$name$cor))){
                cat(txt.strata,"Toeplitz \n\n",sep="")
            }else if(structure$type == "heterogeneous"){
                n.block <- length(unique(structure$X$cor[,3]))-1
                cat(txt.strata,paste0("unstructured with ",n.block," constant subdiagonal",if(n.block>1){"s"}," \n\n"),sep="")
            }else if(tolower(structure$type) == "lag"){
                cat(txt.strata,"block Toeplitz \n\n",sep="")
            }else if(structure$type == "homogeneous"){
                n.block <- length(unique(structure$X$cor[,3]))-1
                cat(txt.strata,paste0("block compound symmetry with ",n.block," specific subdiagonal",if(n.block>1){"s"}," \n\n"),sep="")
            }
        }else if(inherits(structure, "CUSTOM")){
            cat("user-defined structure \n\n")
        }
    }
    ## *** correlation
    if(object$time$n>1 && !hide.cor){
        if(print){
            cat("  - correlation structure:",deparse(formula$cor),"\n")
        }

        ## find unique correlation patterns
        if(identical(type.cor,"param") || (is.null(type.cor) && object$time$n>10)){
            table.cor <- rbind(coef(object,effect="correlation"))
        }else{
            table.cor <- lapply(stats::sigma(object, simplify = FALSE), stats::cov2cor)
        }
        if(print){
            if(identical(type.cor,"param") || (is.null(type.cor) && object$time$n>10)){
                table.cor.print <- table.cor
                rownames(table.cor.print) <- "    "
                print(table.cor.print)
                cat("\n")
            }else{
                table.cor.print  <- lapply(table.cor, function(iCor){ ## iCor <- table.cor[[1]]
                    if(is.matrix(iCor)){
                        if(!is.null(rownames(iCor))){
                            rownames(iCor) <- paste0("    ",rownames(iCor))
                        }
                    }else{
                        iCor <- lapply(iCor, function(iiCor){
                            if(!is.null(rownames(iiCor))){
                                rownames(iiCor) <- paste0("    ",rownames(iiCor))
                            }
                            return(iiCor)
                        })
                    }
                    return(iCor)
                })
                if(length(table.cor)==1){ ## only one (unique) pattern
                    print(table.cor.print[[1]], digits = digits)
                    if(!hide.var || !hide.sd){cat("\n")}
                }else{
                    print(table.cor.print, digits = digits)
                }                        
            }
        }
    }else{
        table.cor <- NULL
    }
    
    ## *** variance
    if(!hide.var || !hide.sd){
        name.sigma <- names(coef(object, transform.k = "sd", effects = "variance"))
        index.ref <- which(names(coef(object, effects = "variance", transform.names = FALSE)) %in% names(param.sigma))
        if(print){
            cat("  - variance structure:",deparse(formula$var),"\n")
        }
    }

    if(!hide.var){
        table.var <- cbind(estimate = coef(object, transform.k = "var", effects = "variance"),
                           estimate.ratio = coef(object, transform.sigma = "none", transform.k = "square", effects = "variance", transform.names = FALSE))
        table.var[index.ref,"estimate.ratio"] <- 1
        test.k <- NROW(table.var) > length(index.ref)
        rownames(table.var) <- name.sigma
    }else{
        table.var <- NULL
    }
    if(!hide.sd){
        table.sd <- cbind(estimate = coef(object, transform.k = "sd", effects = "variance"),
                          estimate.ratio = coef(object, transform.sigma = "none", transform.k = "none", effects = "variance", transform.names = FALSE))
        table.sd[index.ref,"estimate.ratio"] <- 1
        test.k <- NROW(table.sd) > length(index.ref)
        rownames(table.sd) <- name.sigma
    }else{
        table.sd <- NULL
    }
    if(print && (!hide.var || !hide.sd)){
        printtable <- matrix(NA, ncol = 0, nrow = length(name.sigma))
        if(!hide.var){
            printtable <- cbind(printtable, data.frame(variance = unname(table.var[,"estimate"]),stringsAsFactors = FALSE))
            if(test.k){
                printtable <- cbind(printtable, data.frame(ratio = unname(table.var[,"estimate.ratio"]),stringsAsFactors = FALSE))
            }
        }
        if(!hide.sd){
            printtable <- cbind(printtable, data.frame("standard deviation" = unname(table.sd[,"estimate"]),stringsAsFactors = FALSE))
            if(test.k){
                printtable <- cbind(printtable, data.frame(ratio = unname(table.sd[,"estimate.ratio"]),stringsAsFactors = FALSE))
            }
        }
        rownames(printtable) <- paste0("    ",name.sigma)
        print(printtable, digits = digits)
    }
    if(print && !hide.re){
        if(n.strata==1){
            cat("  - variance decomposition: ",deparse(structure.ranef$formula),"\n",sep="")
        }else{
            cat("  - variance decomposition: ",paste0(object$strata$var," ",deparse(structure.ranef$formula)),"\n",sep="")
        }
        table.reA <- nlme::ranef(object, effects = "variance", format = "long", scale = "absolute", simplify = FALSE)
        table.reR <- nlme::ranef(object, effects = "variance", format = "long", scale = "relative", simplify = FALSE)
        if(n.strata==1){
            printtable <- cbind(variance = table.reA$estimate, "%" = 100*table.reR$estimate, sd = sqrt(table.reA$estimate))
            rownames(printtable) <- paste0("    ",table.reA$type)
        }else{
            printtable <- stats::setNames(lapply(U.strata, function(iU){ ## iU <- "Male"
                iPrintTable <- cbind(variance = table.reA[table.reA$strata==iU,"estimate"],
                                     "%" = 100*table.reR[table.reA$strata==iU,"estimate"],
                                     sd = sqrt(table.reA[table.reA$strata==iU,"estimate"]))
                rownames(iPrintTable) <- paste0("    ",table.reA[table.reA$strata==iU,"type"])
                return(iPrintTable)
            }), U.strata)
        }
        print(printtable, digits = digits, na.print = "" , quote = FALSE)
    }else{
        table.reA <- NULL
        table.reR <- NULL
    }

    
    if(print && (!hide.cor || !hide.var || !hide.sd || !hide.re)){
        cat("\n")
    }

    ## ** mean structure
    if(!hide.mean && length(param.mu)>0){
        table.mean <- confint(object,
                              level = level,
                              robust = robust,
                              effects = "mean",
                              columns = c("estimate","se","df","lower","upper","statistic","null","p.value"))

        if(print){
            cat("Fixed effects:",deparse(formula$mean),"\n\n")
            .printStatTable(table = table.mean, df = df, level = level, robust = robust,
                            method.p.adjust = NULL,
                            backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                            columns = columns, 
                            col.df = "df", name.statistic = c("z-statistic","t-statistic"),
                            type.information = type.information,
                            digits = digits,
                            digits.df = digits.df,
                            digits.p.value = digits.p.value,
                            decoration = TRUE, legend = TRUE)
            cat("\n")
        }
    }else{
        table.mean <- NULL
    }
    
    ## ** export
    return(invisible(list(correlation = table.cor,
                          variance = table.var,
                          sd = table.sd,
                          reA = table.reA,
                          reR = table.reR,
                          mean = table.mean)))
}

## * summary.Wald_lmm (documentation)
##' @title Summary of Testing for a Linear Mixed Models
##' @description Estimates, p-values, and confidence intevals for linear hypothesis testing, possibly adjusted for multiple comparisons.
##' 
##' @param object an \code{Wald_lmm} object, output of \code{anova}.
##' @param print [logical] should the output be printed in the console.
##' Can be a vector of length 2 where the first element refer to the global tests and the second to the individual tests.
##' @param seed [integer] value that will be set before adjustment for multiple comparisons to ensure reproducible results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param columns [character vector] Columns to be displayed for each null hypothesis.
##' Can be any of \code{"type"}, \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.##' 
##' @param legend [logical] should explanations about the content of the table be displayed.
##' @param digits [interger, >0] number of digits used to display estimates.
##' @param digits.df [interger, >0] number of digits used to display degrees of freedom.
##' @param digits.p.value [interger, >0] number of digits used to display p-values.
##' @param sep [character] character string used to separate the type of test (e.g. mean, variance) and the name of the test.
##' @param ... arguments \code{method}, \code{level}, and \code{backtransform} passed to \code{\link{confint.Wald_lmm}}
##'
##'
##' @details By default adjustment for multiple comparisons via a single step max-test adjustment,
##'  either using the multcomp package (equal degrees of freedom, \code{method="single-step"}) or the copula package (unequal degrees of freedom, \code{method="single-step2"}).
##' See the argument \code{method} of \code{\link{confint.Wald_lmm}} for other adjustments for multiple comparisons. \cr
##' When multiple multivariate Wald tests are performed, adjustment for multiple comparisons for the univariate Wald tests is performed within each multivariate Wald test.
##' The number of tests ajusted for equal the first degree of freedom of the multivariate Wald statistic. \cr
##'
##' Adding the value \code{"type"} in argument \code{"columns"} ensures that the type of parameter that is being test (mean, variance, correlation) is output.
##'
##' @return \code{NULL}
##' 
##' @keywords methods
 
## * summary.Wald_lmm (code)
##' @export
summary.Wald_lmm <- function(object, print = TRUE, seed = NULL, columns = NULL, legend = TRUE,
                             digits = 3, digits.df = 1, digits.p.value = 3, sep = ": ",
                             ...){

    ## ** normalize input
    if(length(print)==1){
        print.univariate <- print
        print.multivariate <- print
    }else if(length(print)>2){
        stop("Argument \'print\' should have length at most 2. \n",
             "The first element refering to global test and the second to individual hypotheses. \n")
    }else{
        print.multivariate <- print[1]
        print.univariate <- print[2]
    }
    options <- LMMstar.options()
    valid.columns <- c("null","estimate","se","statistic","df","lower","upper","p.value","","type")
    if(identical(columns,"all")){
        columns.multivariate <- setdiff(valid.columns, c("estimate", "se", "lower", "upper"))
        columns.univariate <- valid.columns
    }else  if(is.null(columns)){
        columns.univariate <- options$columns.anova
        if(any(object$multivariate$type!="all")){
            columns.multivariate <- union(c("type","statistic"), setdiff(options$columns.anova, c("estimate", "se", "lower", "upper")))
        }else{
            columns.multivariate <- union(c("statistic"), setdiff(options$columns.anova, c("estimate", "se", "lower", "upper")))
        }
    }else{
        columns.univariate <- tolower(columns)
        if(any(columns.univariate %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns.univariate[columns.univariate %in% valid.columns == FALSE], collapse ="\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns.univariate), collapse ="\" \""),"\".\n")
        }
        if(!is.null(names(columns.univariate)) && all(names(columns.univariate)=="add")){
            columns.univariate <- union(options$columns.anova, unname(columns.univariate))
            columns.multivariate <- setdiff(union("statistic",columns.univariate), c("estimate", "se", "lower", "upper"))
        }else if(!is.null(names(columns.univariate)) && all(names(columns.univariate)=="remove")){
            columns.univariate <- setdiff(options$columns.anova, unname(columns.univariate))
            columns.multivariate <- setdiff(setdiff(union(options$columns.anova,"statistic"), unname(columns.univariate)), c("estimate", "se", "lower", "upper"))
        }else{
            columns.multivariate <- setdiff(columns.univariate, c("estimate", "se", "lower", "upper"))
        }
    }
    if(length(columns.univariate)==0){
        print.univariate <- FALSE
    }
    if(length(columns.multivariate)==0){
        print.multivariate <- FALSE
    }

    if("df" %in% columns.multivariate){
        index.df <- which(columns.multivariate == "df")
        if(index.df == 1){
            columns.multivariate <- c("df.num", "df.denom", columns.multivariate[(index.df+1):length(columns.multivariate)])
        }else if(index.df == length(columns.multivariate)){
            columns.multivariate <- c(columns.multivariate[1:(index.df-1)], "df.num", "df.denom")
        }else{
            columns.multivariate <- c(columns.multivariate[1:(index.df-1)], "df.num", "df.denom", columns.multivariate[(index.df+1):length(columns.multivariate)])
        }
    }
    type.information <- object$object$type.information
    df <- object$args$df
    robust <- object$args$robust
    ci <- object$args$ci

    transform.sigma <- object$args$transform.sigma
    transform.k <- object$args$transform.k
    transform.rho <- object$args$transform.rho

    ## ** ensure reproducibility
    if(!is.null(seed)){
        old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit( assign(".Random.seed", old.seed, envir = .GlobalEnv, inherits = FALSE) )
        set.seed(seed)
    }

    ## ** extract information
    ## *** multivariate tests
    if(print.multivariate>0 && any(!is.na(object$multivariate$df.num))){
        table.multivariate <- object$multivariate[,setdiff(columns.multivariate,c("type","")),drop=FALSE]
        nchar.type <- nchar(object$args$type[[1]])
        maxchar.type <- max(nchar.type)
        if("type" %in% columns.multivariate){
            vec.white <- sapply(maxchar.type-nchar.type, function(iN){paste(rep(" ", iN), collapse = "")})
            sparsetype <- object$args$type[[1]]            
            sparsetype[duplicated(sparsetype)] <- sapply(sparsetype[duplicated(sparsetype)], function(iWord){paste(rep(" ", times = nchar(iWord)), collapse="")})
            rownames(table.multivariate) <- paste(vec.white,sparsetype,sep,object$multivariate$test,sep="")
        }else if(any(duplicated(object$multivariate$test))){
            test.nduplicated <- !duplicated(object$args$type[[1]])
            rownames(table.multivariate) <- sapply(1:NROW(table.multivariate), function(iIndex){
                if(test.nduplicated[iIndex]){
                    paste(paste(rep(" ", maxchar.type-nchar.type[iIndex]), collapse = ""),object$args$type[[1]][iIndex],sep,object$multivariate$test[iIndex],sep="")
                }else{
                    paste(paste(rep(" ", maxchar.type+nchar(sep)), collapse = ""),object$multivariate$test[iIndex],sep="")
                }
            })
        }else if(any(!identical(object$args$type[[1]],"all"))){
            rownames(table.multivariate) <- object$multivariate$test
        }
        if(print.multivariate>0.5){
            cat("\t\tMultivariate Wald test \n\n")
        }

        .printStatTable(table = table.multivariate, robust = robust, df = df, level = NULL, type.information = type.information,
                        method.p.adjust = NULL,
                        backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                        columns = setdiff(columns.multivariate,"type"), col.df = c("df.num","df.denom"), name.statistic = c("Chi2-statistic","F-statistic"),
                        digits = digits, digits.df = c(df.num = 0, df.denom = digits.df), digits.p.value = digits.p.value,
                        decoration = legend, legend = legend)
        
        cat("\n")
    }

    ## *** local tests
    if(print.univariate>0 && ci){

        table.univariate <- confint(object, columns = union(setdiff(columns.univariate,""),"type"), ...)
        if(is.null(columns) && all(is.na(table.univariate$lower)) && all(is.na(table.univariate$upper))){
            columns.univariate <- setdiff(columns.univariate, c("lower","upper"))
        }

        error <- attr(table.univariate,"error")
        n.sample <- attr(table.univariate,"n.sample")
        method.p.adjust <- attr(table.univariate,"method")
        if(length(method.p.adjust)>1){
            method.p.adjust <- method.p.adjust[1]
            message("Different adjustment for multiple comparisons were used but only the first is described. \n")
        }
        level <- attr(table.univariate,"level")
        backtransform <- attr(table.univariate,"backtransform")

        ## type of adjustment
        if(length(object$args$type[[1]])>1){
            if(length(unique(object$args$type[[1]]))==1){
                factor.p.adjust <- "covariate name"
            }else if(all(duplicated(object$args$type[[1]])==FALSE)){
                factor.p.adjust <- "type of parameter"
            }else{
                factor.p.adjust <- "covariate name and type of parameter"
            }
        }else{
            factor.p.adjust <- NULL
        }

        ## incorporate type
        if(method.p.adjust %in% c("average","pool.fixse","pool.se","pool.gls","pool.gls1","pool.rubin") == FALSE){
            table.univariate$type <-  c("mu" = "mean", "sigma" = "variance", "k" = "variance", "rho" = "correlation")[object$univariate$type]
        }
        nchar.type <- nchar(table.univariate$type)
        maxchar.type <- max(nchar.type)
        if("type" %in% columns.univariate){
            vec.white <- sapply(maxchar.type-nchar.type, function(iN){paste(rep(" ", iN), collapse = "")})
            rownames(table.univariate) <- paste(vec.white,table.univariate$type,sep,rownames(table.univariate),sep="")
        }else if(length(unique(table.univariate$type))>1){
            test.nduplicated <- !duplicated(table.univariate$type)            
            rownames(table.univariate) <- sapply(1:NROW(table.univariate), function(iIndex){ ## iIndex <- 1
                if(test.nduplicated[iIndex]){
                    paste(paste(rep(" ", maxchar.type-nchar.type[iIndex]), collapse = ""),table.univariate$type[iIndex],sep,rownames(table.univariate)[iIndex],sep="")
                }else{
                    paste(paste(rep(" ", maxchar.type+nchar(sep)), collapse = ""),rownames(table.univariate)[iIndex],sep="")
                }
            })
        }
        table.univariate$type <- NULL
        if(inherits(object,"rbindWald_lmm") && length(unique(setdiff(object$univariate$method,"none"))>1)){
            warning("Different methods have been used to adjust for multiple comparisons - text describing the adjustment will not be accurate.")
        }
        if(print.univariate>0.5){
            cat("\t\tUnivariate Wald test \n\n")
        }
        if(attr(object$object,"independence")==FALSE && method.p.adjust == "pool.se"){
            attr(method.p.adjust,"warning") <- "WARNING: uncertainty about the weights assumes independence between parameters from different models.\n"
        }else if(method.p.adjust == "pool.fixse"){
            attr(method.p.adjust,"warning") <- "WARNING: uncertainty about the weights has been ignored.\n"
        }else if(method.p.adjust %in% c("pool.gls","pool.gls1")){
            attr(method.p.adjust,"warning") <- "WARNING: uncertainty about the weights has been ignored.\n"
            if(!is.null(error) && any(!is.na(error))){
                attr(method.p.adjust,"warning") <- paste0(attr(method.p.adjust,"warning"),
                                                          "           ",error," principal components have been ignored when pooling (singular vcov).\n")
            }
        }

        if("single-step2" %in% method.p.adjust){
            digits.p.value2 <- c(digits.p.value,1/n.sample)
        }else{
            digits.p.value2 <- digits.p.value
        }
        .printStatTable(table = table.univariate, robust = robust, df = df, level = level, type.information = type.information,
                        method.p.adjust = method.p.adjust, factor.p.adjust = factor.p.adjust, error.p.adjust = error, seed = seed, n.sample = n.sample,
                        backtransform = backtransform, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                        columns = setdiff(columns.univariate,"type"), col.df = c("df"), name.statistic = c("t-statistic","z-statistic"),
                        digits = digits, digits.df = digits.df, digits.p.value = digits.p.value2,
                        decoration = legend, legend = legend)
        
        if(print.univariate>0.5){
            cat("\n")
        }


    }

    ## ** export
    return(invisible(NULL))
}


## * summary.LRT_lmm
##' @export
summary.LRT_lmm <- function(object, digits = 3, digits.df = 1, digits.p.value = 3, columns = NULL, ...){

    ## ** normalize input
    valid.columns <- c("null","logLikH0","logLikH1","statistic","df","p.value","")

    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(is.null(columns) || identical(columns,"all")){
        columns <- valid.columns
    }else{
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse ="\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse ="\" \""),"\".\n")
        }
    }

    ## ** display
    cat("\t\tLikelihood ratio test \n")

    name1 <- deparse(attr(object,"call")$object)
    name2 <- deparse(attr(object,"call")$effects)
    if(attr(object,"type")=="1-2"){
        cat("\t\t(",name1," vs. ",name2,")\n\n",sep="")
    }else if(attr(object,"type")=="2-1"){
        cat("\t\t(",name2," vs. ",name1,")\n\n",sep="")
    }
    table <- as.data.frame(object)
    if("null" %in% columns){
        table$null
        cat("  ","Null hypothesis: ",table$null,".\n\n",sep="")
    }

    table[names(table)[names(table) %in% setdiff(columns,"null") == FALSE]] <- NULL
    rownames(table) <- ""
    
    .printStatTable(table = table, robust = NULL, df = FALSE, level = NULL, type.information = NULL,
                    method.p.adjust = NULL,
                    backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                    columns = setdiff(columns,"null"), col.df = c("df"), name.statistic = c("statistic"),
                    digits = digits, digits.df = digits.df, digits.p.value = digits.p.value,
                    decoration = TRUE, legend = TRUE)

    cat("\n")

    ## ** export
    return(invisible(NULL))
}

## * summary.effect_lmm (code)
##' @export
summary.effect_lmm <- function(object, columns = NULL, print = TRUE, ...){

    ## ** normalize user input
    if("print" %in% names(match.call())==FALSE && all(is.na(object$multivariate$df.num))){
        print <- c(0,0.5)        
    }
    if(object$args$effect[[1]][1]=="identity" && "columns" %in% names(match.call())==FALSE){
        object$univariate$p.value <- NULL
        columns <- c("estimate","se","df","lower","upper")
    }

    ## ** prepare
    outcome.txt <- switch(object$args$effect[[1]][2],
                          "none" = "outcome",
                          "change" = "change in outcome",
                          "auc" = "area under the outcome curve",
                          "auc-b" = "area under the outcome curve above baseline")
    contrast.txt <- switch(object$args$effect[[1]][1],
                           "identity" = "Average",
                           "difference" = "Difference in average")

    ## ** display
    if(is.null(object$args$variable)){
        cat("\t\t",contrast.txt," counterfactual ",outcome.txt,"\n\n", sep = "")
    }else{
        cat("\t\t",contrast.txt," counterfactual ",outcome.txt,"\n\t\t w.r.t \'",object$args$variable,"\' values \n\n", sep = "")
    }
    summary.Wald_lmm(object, print = print, columns = columns, ...)
}

## * summary.mlmm (documentation)
##' @title Summary of Multiple Linear Mixed Models
##' @description Estimates, p-values, and confidence intevals for multiple linear mixed models.
##' 
##' @param object an \code{mlmm} object, output of \code{mlmm}.
##' @param digits [integer,>0] number of digits used to display numeric values.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}.
##' @param print [logical] should the output be printed in the console.
##' Can be a vector of length 2 where the first element refer to the global tests and the second to the individual tests.
##' @param hide.data [logical] should information about the dataset not be printed.
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param ... other arguments are passed to \code{\link{summary.Wald_lmm}}.
##'
##' @keywords methods

## * summary.mlmm (code)
##' @export
summary.mlmm <- function(object, digits = 3, method = NULL, print = NULL, hide.data = FALSE, hide.fit = FALSE, ...){

    options <- LMMstar.options()

    if(is.null(method)){
        method <- "none"
    }
    if(is.null(print)){
        print <- c(0,1/2)
    }

    ## extract models
    ls.model <- object$model
    method.fit <- object$object$method.fit
    optimizer <- ls.model[[1]]$args$control$optimizer
    logLik <- sapply(ls.model, logLik)
    cv <- sapply(ls.model, function(iM){iM$opt$cv})
    n.iter <- sapply(ls.model, function(iM){iM$opt$n.iter})
    param.value <- coef(object, effects = "all")

    nparam.mu  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="mu")})
    nparam.sigma  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="sigma")})
    nparam.k  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="k")})
    nparam.rho  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="rho")})

    M.nobs <- stats::nobs(object)
    call <- attr(object,"call")

    ## ** welcome message
    if(any(print>0)){
        cat("	Linear Mixed Models stratified according to \"",eval(call$by),"\" \n\n",sep="")
    }

    ## ** data message    
    if(!hide.data){
        if(inherits(call$data,"call")){
            cat("Dataset:", deparse(call$data), "\n")
        }
        cat("Strata : \"",paste(names(ls.model),collapse = "\", \""),"\"\n\n",sep="")
        if(any(M.nobs[,"missing.obs"]>0)){
            if(any(M.nobs[,"missing.cluster"])){
                cat("  - ", paste(M.nobs[,"cluster"], collapse=", "), " clusters were analyzed \n",
                    "    ", paste(M.nobs[,"missing.cluster"], collapse=", "), " were excluded because of missing values \n" , sep = "")
            }else{
                cat("  - ", paste(M.nobs[,"cluster"], collapse=", "), " clusters \n" , sep = "")
            }
            cat("  - ", paste(M.nobs[,"obs"], collapse = ", "), " observations were analyzed \n",
                "    ", paste(M.nobs[,"missing.obs"],collapse=", "), " were excluded because of missing values \n",  sep = "")
        }else{
            cat("  - ", paste(M.nobs[,"cluster"], collapse=", "), " clusters \n" , sep = "")
            cat("  - ", paste(M.nobs[,"obs"], collapse = ", "), " observations \n", sep = "")
        }
        cat("\n")
    }
    
    ## ** optim message    
    if(!hide.fit){
        cat("Estimation procedure \n\n")
        if(method.fit == "REML"){
            cat("  - Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - Maximum Likelihood (ML) \n")
        }
        cat("  - log-likelihood :", paste(round(as.double(logLik), digits = digits),collapse = ", "), "\n",sep="")
        cat("  - parameters: mean = ",paste(nparam.mu,collapse = ", "),", variance = ",paste(nparam.sigma+nparam.k,collapse = ", "),", correlation = ",paste(nparam.rho,collapse = ", "),"\n", sep = "")
        cat("  - convergence: ",paste(cv>0,collapse = ", ")," (",paste(n.iter,collapse = ", ")," iterations) \n", sep = "")
        cat(" \n")
    }

    ## ** extract test
    if(any(print>0)){
        name.param <- unique(object$univariate$parameter)
        if(length(name.param)==1){
            cat("Statistical inference for ",name.param," \n\n",sep="")
        }else{
            cat("Statistical inference \n\n")
        }
        out <- summary.Wald_lmm(object, method = method, print = print, ...)
    }

    ## ** export
    return(invisible(out))
}

## * summary.partialCor
##' @title Summary for partial correlation
##' @description Display estimated partial correlation and associated p-values and confidence intevals.
##' 
##' @param object a \code{partialCor} object, output of \code{partialCor}.
##' @param digits [integer,>0] number of digits used to display numeric values.
##' @param detail [integer,>0] passed to \code{print.confint_lmm}. If above 0.5 also display when a back-transformation has been used.
##' @param ... other arguments are passed to \code{print.confint_lmm}.
##'
##' @keywords methods
##' 
##' @export
summary.partialCor <- function(object, digits = 3, detail = TRUE, ...){

    args <- attr(object,"args")

    cat("\n\t\tPartial correlation \n\n")

    message.backtransform <- attr(object,"backtransform")
    attr(object,"backtransform") <- NULL
    ## display estimates
    if(!is.null(attr(object,"parameter")) && length(attr(object,"parameter"))==1){
        cat("\tParameter: ", attr(object,"parameter"),"\n\n",sep="")
        rownames(object) <- NULL
    }
    attr(detail,"summary") <- TRUE

    .printStatTable(table = object, df = args$df, level = args$level, robust = FALSE,
                    method.p.adjust = NULL,
                    backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                    columns = names(object),
                    col.df = "df", name.statistic = c("z-statistic","t-statistic"),
                    type.information = NULL,
                    digits = digits,
                    digits.df = 1,
                    digits.p.value = digits,
                    decoration = TRUE, legend = TRUE)

    ## legend (transformation)
    test.backtransform <- !is.null(message.backtransform) && any(!is.na(message.backtransform$FUN))
    if(test.backtransform){
        message.backtransform <- message.backtransform[!is.na(message.backtransform$FUN),,drop=FALSE]

            if(any(message.backtransform[,setdiff(names(message.backtransform), "FUN")] == FALSE)){
                warning("Could not back-transform everything.\n")
            }

        if(NROW(object)==1){
            short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
            txt <- unique(short2text[intersect(names(short2text),intersect(names(object),names(message.backtransform)))])
        }else{
            short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
            txt <- unique(short2text[intersect(names(short2text),intersect(names(object),names(message.backtransform)))])
        }
        cat("  ",paste(txt,collapse = ", ")," have been back-transformed",sep="")
        if(detail>=0.5){
            cat(" (",paste(message.backtransform$FUN,collapse="/"),"). \n", sep ="")
        }
    }
    cat("\n")

    if(!is.null(attr(object,"R2"))){
        cat("\t\tCoefficient of determination (R2)\n\n")

        table.R2 <- attr(object,"R2")
        .printStatTable(table = table.R2, df = args$df, level = args$level, robust = FALSE,
                        method.p.adjust = NULL,
                        backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                        columns = names(object),
                        col.df = "df", name.statistic = c("z-statistic","t-statistic"),
                        type.information = NULL,
                        digits = digits,
                        digits.df = 1,
                        digits.p.value = digits,
                        decoration = TRUE, legend = TRUE)

        cat("\n")
    }

    return(invisible(NULL))
}



## * print.resample
##' @export
summary.resample <- function(object, digits = 3, ...){
    args <- attr(object,"args")
    n.sample <- attr(object,"n.sample")

    if(args$type %in% c("perm-var","perm-res")){
        if(args$studentized){
            cat("\tStudentized permutation test\n\n", sep = "")
        }else{
            cat("\tPermutation test\n\n", sep = "")
        }
    }else if(args$type == "boot"){
        if(args$studentized){
            cat("\tNon-parametric studentized bootstrap\n\n", sep = "")
        }else{
            cat("\tNon-parametric bootstrap\n\n", sep = "")
        }
    }

    
    base::print.data.frame(object, digits = digits)

    cat(rep("-",ncharTable(object, digits = digits)),"\n",sep="")
    cat(paste0("(based on ",n.sample," samples - ",round((1-n.sample/args$n.sample)*100, digits = digits),"% failed) \n"))
    cat("\n")
    return(invisible(NULL))
}

## * .printStatTable (documentation)
##' @description Display a table containing the model coefficients and their uncertainty, as well as a legendn.
##' Inspired from \code{\link{stats::printCoefmat}}.
##'
##' @param table [data.frame] table containing the coefficients to be displayed.
##' @param robust [logical or NULL] are robust standard error used?
##' @param df [logical or NULL] are degrees of freedom calculated by Satterthwaite approximation?
##' @param level [numeric or NULL] confidence level.
##' @param type.information [character] type of information matrix.
##' @param method.p.adjust [character or NULL] adjustment method for multiple comparisons.
##' \code{"none"} corresponds to no adjustment.
##' @param factor.p.adjust [character or NULL] Are p-values adjusted within a certain factor?
##' @param error.p.adjust [numeric or NULL] Numeric error performed when adjusting for multiple comparisons.
##' @param seed [integer, >0] Random number generator (RNG) state used when adjusting for multiple comparisons.
##' @param n.sample [integer, >0] Number of samples used  when adjusting for multiple comparisons.
##' @param backtransform [data.frame or NULL] Should estimates and their uncertainty be back-transformed?
##' @param transform.sigma,transform.k,transform.rho [character or NULL] Transformation used on certain type of parameters.
##' @param decoration [logical] should a decoration be displayed between the table and the legend?
##' @param columns [character vector] columns from argument \code{table} to be displayed.
##' @param col.df [character vector] columns containing the degrees of freedom. If several, they will be merged.
##' @param name.statistic [character vector] how to rename the statistic depending on whether degrees of freedom have been computed.
##' @param digits [interger, >0] number of digits used to display estimates.
##' @param digits.df [interger, >0] number of digits used to display degrees of freedom.
##' @param digits.p.value [interger, >0] number of digits used to display p-values.
##' @param decoration [logical] should an horizontal bar be displayed just after the table.
##' @param legend [logical] should explanations about the content of the table be displayed.
##' @param space [character] horizontal space.
##' 

## * .printStatTable (code)
##' @noRd
.printStatTable <- function(table, robust, df, level, type.information,
                            method.p.adjust = NULL, factor.p.adjust, error.p.adjust, seed, n.sample,
                            backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                            columns, col.df, name.statistic,
                            digits, digits.df, digits.p.value,
                            decoration, legend, space = "  "){

    ## ** check input
    if(any(setdiff(columns,"") %in% names(table) == FALSE)){
        missing.col <- setdiff(columns,"")[setdiff(columns,"") %in% names(table) == FALSE]
        stop("Inconsistency between argument \'columns\' and \'table\'. \n",
             "Could not find column(s) \"",paste(missing.col, collapse="\" \""),"\" in \'table\'. \n")
    }
    if(is.null(transform.sigma)){
        transform.sigma <- NA
    }
    if(is.null(transform.k)){
        transform.k <- NA
    }
    if(is.null(transform.rho)){
        transform.rho <- NA
    }

    if(length(digits.p.value)==1){
        digits.p.value <- c(digits.p.value,1e-4)
    }

    ## ** add stars    
    if("p.value" %in% names(table)){
        if(all(is.na(table$p.value))){
            table.print <- table[,setdiff(names(table),"p.value"),drop=FALSE]
        }else{
            table.print <- cbind(table,
                                 stats::symnum(table$p.value, corr = FALSE, na = FALSE,
                                               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                               symbols = c("***", "**", "*", ".", " "))
                                 )
            names(table.print)[NCOL(table.print)] <- ""
        }
    }else{
        table.print <- table
    }
    

    ## ** round
    columns.num <- intersect(setdiff(columns,c(col.df,"null","p.value")), names(table))
    for(iCol in columns.num){ 
        table.print[[iCol]] <- as.character(round(table[[iCol]], digits = digits))
        if(any(!is.na(table[[iCol]])) && any(table[[iCol]]!=0 & table.print[[iCol]] == "0")){
            table.print[[iCol]][table[[iCol]]!=0 & table.print[[iCol]] == "0"] <- paste0("<0.",paste0(rep(0,digits-1),collapse=""),"1")
        }
    }
    if("p.value" %in% names(table.print)){
        table.print$p.value <- as.character(format.pval(table.print$p.value, digits = digits.p.value[1], eps = digits.p.value[2]))
    }
    if(!is.null(col.df) && all(col.df %in% columns)){
        if(identical(col.df,"df")){
            table.print$df <- as.character(round(table.print$df, digits = digits.df))
        }else{
            table.print[[col.df[1]]] <- paste0("(",apply(do.call(cbind,lapply(col.df, function(iCol){
                iDF <- formatC(table.print[[iCol]], digits = digits.df[[iCol]], format = "f")
                nchar.iDF <- nchar(iDF)
                maxchar.iDF <- max(nchar.iDF)
                if(all(nchar.iDF==maxchar.iDF)){
                    return(iDF)
                }else{
                    iWS <- sapply(maxchar.iDF-nchar.iDF, function(iN){paste(rep(" ",iN),collapse="")})
                    return(paste0(iWS,iDF))
                }
            })), 1, paste, collapse = ","),")")
            table.print[col.df[-1]] <- NULL
            names(table.print)[names(table.print)== col.df[1]] <- "df"
            columns[columns==col.df[1]] <- "df"
            columns <- setdiff(columns,col.df[-1])
        }
    }

    ## ** rename statistic
    if(!is.null(name.statistic) && ("df" %in% names(table.print)) && ("statistic" %in% names(table.print))){
        if(df){
            columns[columns=="statistic"] <- name.statistic[2]
            names(table.print)[names(table.print)=="statistic"] <- name.statistic[2]
        }else{
            columns[columns=="statistic"] <- name.statistic[1]
            names(table.print)[names(table.print)=="statistic"] <- name.statistic[1]
        }
    }

    ## ** add space to rownames
    if(identical(rownames(table.print),"1")){
        rownames(table.print) <- space
    }else{
        rownames(table.print) <- paste(space,rownames(table.print))
    }
    ## ** subset columns
    table.print <- table.print[,names(table.print) %in% columns,drop=FALSE]

    ## ** print table
    print(table.print)

    ## ** add decoration below table
    if(decoration){
        cat(space,paste(rep("-", ncharTable(table.print, digits = digits)),collapse=""),"\n")
    }
    
    ## ** legend
    if(legend){

        ## *** significance level
        if("" %in% columns){
            cat(space,"Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1.\n",sep="")
        }

        if(!is.null(level) && is.null(method.p.adjust) && any(c("lower", "upper") %in% columns)){
            ## *** ci 
            if(NROW(table)==1){
                txt.ci <- paste0("the ",100*level,"% confidence intervals of the coefficient")
            }else{
                txt.ci <- paste0(100*level,"% pointwise confidence intervals for each coefficient")
            }
            if("lower" %in% columns && "upper" %in% columns){
                cat(space,"Columns lower and upper contain ",txt.ci,".\n", sep = "")
            }else if("lower" %in% columns){
                cat(space,"Column lower contains ",txt.ci,".\n", sep = "")
            }else if("upper" %in% columns){
                cat(space,"Column upper contains ",100*level,"% ",txt.ci,".\n", sep = "")
            }
        }else if(!is.null(method.p.adjust) && method.p.adjust %in% c("average","pool.se","pool.fixse","pool.gls","pool.gls1","pool.rubin")){

            if(method.p.adjust == "average"){
                if(NROW(table)>1){
                    cat(space,"Estimates have been averaged within ",factor.p.adjust,".\n", sep="")
                }else{
                    cat(space,"Estimates have been averaged.\n", sep="")
                }
            }else if(method.p.adjust %in% c("pool.se","pool.fixse")){
                if(NROW(table)>1){
                    cat(space,"Estimates have been averaged, weighted by the inverse of their variance, within ",factor.p.adjust,".\n", sep="")
                }else{
                    cat(space,"Estimates have been averaged, weighted by the inverse of their variance.\n", sep="")
                }
                if(any(c("se","lower","upper","p.value") %in% columns) && !is.null(attr(method.p.adjust,"warning"))){
                    cat(space,attr(method.p.adjust,"warning"),sep="")
                }
            }else if(method.p.adjust %in% c("pool.gls","pool.gls1")){
                if(NROW(table)>1){
                    cat(space,"Estimates have been averaged, weighted via GLS ",factor.p.adjust,".\n", sep="")
                }else{
                    cat(space,"Estimates have been averaged, weighted via GLS.\n", sep="")
                }
                if(any(c("se","lower","upper","p.value") %in% columns) && !is.null(attr(method.p.adjust,"warning"))){
                    cat(space,attr(method.p.adjust,"warning"),sep="")
                }
            }else if(method.p.adjust == "pool.rubin"){
                if(NROW(table)>1){
                    cat(space,"Estimates have been pooled within ",factor.p.adjust," using Rubin's rule.\n", sep="")
                }else{
                    cat(space,"Estimates have been pooled using Rubin's rule.\n", sep="")
                }
            }
        }else if(!is.null(level) && !is.null(method.p.adjust) && any(c("p.value","lower", "upper") %in% columns) && NROW(table)>1){

            ## *** adjustment for multiple comparisons 
            txt.cip <- paste("Columns",paste(intersect(c("lower","upper","p.value"),columns),collapse="/"))

            if(method.p.adjust == "none"){ ## 
                ## cat(space,txt.cip," not adjusted for multiple comparisons.\n", sep="")
            }else{
                if(!is.null(attr(table,"Madjust-within")) && attr(table,"Madjust-within")){
                    if(method.p.adjust %in% c("single-step", "single-step2")){
                        cat(space,txt.cip," adjusted for multiple comparisons (within coefficient) -- max-test.\n", sep="")
                    }else{
                        cat(paste0(space,txt.cip," adjusted for multiple comparisons (within coefficient) -- ",method.p.adjust,".\n", sep=""),sep="")
                    }
                }else{
                    if(method.p.adjust %in% c("single-step", "single-step2")){
                        cat(space,txt.cip," adjusted for multiple comparisons -- max-test.\n", sep="")
                    }else{
                        cat(paste0(space,txt.cip," adjusted for multiple comparisons -- ",method.p.adjust,".\n", sep=""),sep="")
                    }
                }
            }

            if(method.p.adjust != "none" && !is.null(factor.p.adjust)){
                cat(space,"(adjustment within ",factor.p.adjust,"). \n",sep="")
            }
            if(method.p.adjust == "single-step"){
                if(!is.null(error.p.adjust) && any(!is.na(error.p.adjust)) && any(abs(stats::na.omit(error.p.adjust))>1e-12)){
                    if(!is.null(seed)){
                        cat(space,"(error when computing the adjusted ",tolower(txt.cip)," by numerical integration: ", signif(max(error.p.adjust, na.rm=TRUE), digits = digits.p.value[1])," with seed for the RNG: ",seed,")\n",sep="")
                    }else{
                        cat(space,"(error when computing the adjusted ",tolower(txt.cip)," by numerical integration: ", signif(max(error.p.adjust, na.rm=TRUE), digits = digits.p.value[1]),")\n",sep="")
                    }
                }else if(!is.null(seed)){
                    cat(space,"(seed for the RNG: ",seed,")\n",sep="")
                }
            }else if(method.p.adjust == "single-step2"){
                if(!is.null(seed)){
                    cat(space,"(",n.sample," samples have been used with seed for the RNG",seed,")\n",sep="")
                }else{
                    cat(space,"(",n.sample," samples have been used)\n",sep="")
                }
            }
        }

        ## *** type of standard error 
        if(!is.null(robust) && !is.null(type.information)){
            if(robust && "se" %in% columns){
                cat(space,"Robust standard errors are derived from the ",type.information," information (column se). \n", sep = "")
            }else if("se" %in% columns){
                cat(space,"Model-based standard errors are derived from the ",type.information," information (column se). \n", sep = "")
            }
        }

        ## *** type of degree of freeedom 
        if(identical(df,TRUE) && "df" %in% columns){
            cat(space,"Degrees of freedom were computed using a Satterthwaite approximation (column df). \n", sep = "")
        }

        ## *** backtransformation
        if(!is.null(backtransform)){
            message.backtransform <- backtransform[!is.na(backtransform$FUN),,drop=FALSE]

            if(any(message.backtransform[,setdiff(names(message.backtransform), "FUN")] == FALSE)){
                warning("Could not back-transform everything.\n")
            }
            if(NROW(table)==1){
                short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),intersect(columns,names(message.backtransform)))])
            }else{
                short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),intersect(columns,names(message.backtransform)))])
            }
            substr(txt[1], 1, 1) <- toupper(substr(txt[1], 1, 1))
            cat("  ",paste(txt,collapse = ", ")," have been back-transformed",sep="")
            cat(" (",paste0(paste(rownames(message.backtransform),collapse = "/")," parameters with ",paste(message.backtransform$FUN,collapse="/")),"). \n", sep ="")
        }else if(any(!is.na(c(transform.sigma,transform.k,transform.rho)))){
            vec.transform <- stats::na.omit(c(sigma=transform.sigma,k=transform.k,rho=transform.rho))
            if(NROW(table)==1){
                short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),columns)])
            }else{
                short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),columns)])
            }
            substr(txt[1], 1, 1) <- toupper(substr(txt[1], 1, 1))
            cat("  ",paste(txt,collapse = ", ")," have been transformed (",paste0(paste(names(vec.transform),collapse = "/")," parameters with ",paste(vec.transform,collapse="/")),"). \n", sep ="")
        }
    }

    ## ** export
    return(invisible(NULL))

}


######################################################################
### summary.R ends here
 
   
