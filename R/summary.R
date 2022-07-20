### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: jul 20 2022 (17:05) 
##           By: Brice Ozenne
##     Update #: 571
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
##' This is a modified version of the \code{nlme::summary.gls} function.
##' @name summary
##'
##' @param object [lmm] output of the \code{lmm} function.
##' @param digits [integer,>0] number of digits used to display numeric values.
##' When a vector of length 2, the second values is used for the number of digits for p-values.
##' @param level [numeric,0-1] confidence level for the confidence intervals.
##' @param type.cor [character] should the correlation matrix be display (\code{"matrix"}) or the parameter values (\code{"param"}).
##' @param print [logical] should the output be printed in the console.
##' @param columns [character vector] Columns to be output for the fixed effects. Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param hide.data [logical] should information about the dataset not be printed.
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param hide.cor [logical] should information about the correlation structure not be printed.
##' @param hide.sd [logical] should information about the standard deviation not be printed.
##' @param hide.var [logical] should information about the variance not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.
##'
##' @return A list containing elements displayed in the summary: \itemize{
##' \item \code{correlation}: the correlation structure.
##' \item \code{variance}: the variance structure.
##' \item \code{sd}: the variance structure expressed in term of standard deviations.
##' \item \code{mean}: the mean structure.
##' }

## * summary.lmm (code)
##' @rdname summary
##' @export
summary.lmm <- function(object, digits = 3, level = 0.95, type.cor = NULL, robust = FALSE, print = TRUE, columns = NULL,
                        hide.data = FALSE, hide.fit = FALSE, hide.cor = is.null(object$formula$cor), hide.var = TRUE, hide.sd = FALSE, hide.mean = FALSE, ...){

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

    logLik <- stats::logLik(object)
    nobs <- stats::nobs(object)
    method.fit <- object$method
    nobsByCluster <- object$design$cluster$nobs
    formula <- object$formula
    df <- !is.null(object$df)
    options <- LMMstar.options()

    n.cluster.original <- object$cluster$n
    n.cluster.design <- object$design$cluster$n
    
    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    valid.columns <- c("estimate","se","statistic","df","lower","upper","null","p.value","")
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
        }
    }else{
        columns <- options$columns.summary
    }    

    if(!is.null(type.cor)){
        type.cor <- match.arg(type.cor, c("matrix","param"))
    }
    if(length(digits)==1){
        digits <- rep(digits,2)
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
        cat("Dataset:", deparse(call$data), "\n\n")

        if(nobs["missing"]>0){
            if(n.cluster.original-n.cluster.design>0){
                cat("  - ", nobs["cluster"], " clusters were analyzed, ",n.cluster.original-n.cluster.design," were excluded because of missing values \n" , sep = "")
            }else{
                cat("  - ", nobs["cluster"], " clusters \n" , sep = "")
            }
            cat("  - ", sum(nobsByCluster), " observations were analyzed, ",nobs["missing"]," were excluded because of missing values \n",  sep = "")
        }else{
            cat("  - ", nobs["cluster"], " clusters \n" , sep = "")
            cat("  - ", sum(nobsByCluster), " observations \n",  sep = "")
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
        if(object$opt$name!="gls"){
            abs.score <- abs(object$score)
            abs.diff <- abs(object$opt$previous.estimate-object$param)
            name.score <- names(which.max(abs.score))[1]
            name.diff <- names(which.max(abs.diff))[1]
            
            cat("  - convergence: ",object$opt$cv>0," (",object$opt$n.iter," iterations) \n",
                "    largest |score| = ",max(abs.score)," for ",name.score,"\n",
                if(!is.null(name.diff)){paste0("            |change|= ",max(abs.diff)," for ",name.diff,"\n")},
                sep = "")
        }
        cat(" \n")
    }

    ## ** vcov structure
    if(print && (!hide.cor || !hide.var || !hide.sd)){
        cat("Residual variance-covariance: ")
        if(is.na(structure$name$strata)){
            txt.strata <- ""
        }else{
            txt.strata <- "stratified "
        }
        
        if(length(param.rho)==0){
            if(length(c(param.sigma,param.k))==1){
                cat(txt.strata,"identity \n\n",sep="")
            }else{
                cat(txt.strata,"diagonal \n\n",sep="")
            }
        }else if(structure$type == "UN"){
            cat(txt.strata,"unstructured \n\n",sep="")
        }else if(structure$type == "CS"){
            if(all(is.na(structure$name$cor[[1]]))){
                cat(txt.strata,"compound symmetry \n\n",sep="")
            }else if(structure$heterogeneous){
                cat(txt.strata,"block unstructured \n\n",sep="")
            }else{
                cat(txt.strata,"block compound symmetry \n\n",sep="")
            }
        }else if(structure$type == "TOEPLITZ"){
            cat(txt.strata,"Toeplitz \n\n",sep="")
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
            table.cor <- lapply(stats::sigma(object, simplifies = FALSE), stats::cov2cor)
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
                    print(table.cor.print[[1]], digits = digits[1])
                    cat("\n")
                }else{
                    print(table.cor.print, digits = digits[1])
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
            cat("  - variance structure:",deparse(formula$var.design),"\n")
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
        print(printtable, digits = digits[1])
    }
    if(print && (!hide.cor || !hide.var || !hide.sd)){
        cat("\n")
    }
    
    ## ** mean structure
    if(!hide.mean){
        if(print){
            cat("Fixed effects:",deparse(call$formula),"\n\n")
        }
        table.mean <- confint(object,
                              level = level,
                              robust = robust,
                              effects = "mean",
                              columns = c("estimate","se","df","lower","upper","statistic","p.value"))

        printtable.mean <- table.mean[,-which(names(table.mean) %in% columns == FALSE),drop=FALSE]
        if(length(object$df)>0){
            names(printtable.mean) <- gsub("^statistic","t-statistic",names(printtable.mean))
        }else{
            names(printtable.mean) <- gsub("^statistic","z-statistic",names(printtable.mean))
        }
   
        if(print){
            toPrint <- capture.output(stats::printCoefmat(printtable.mean, digits = digits[1],
                                                          has.Pvalue = "p.value" %in% columns,
                                                          P.values = "p.value" %in% columns,
                                                          eps.Pvalue = 10^{-digits[2]},
                                                          signif.legend = FALSE))
            sapply(toPrint, function(iTxt){cat("  ",iTxt,"\n",sep="")})
            cat(" ",paste(rep("-", max(nchar(toPrint))),collapse=""),"\n")
            if("" %in% columns){
                cat("  Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n")
            }
            if(robust && "se" %in% columns){
                cat("  Uncertainty was quantified using robust standard errors (column se). \n", sep = "")
            }else if("se" %in% columns){
                cat("  Uncertainty was quantified using model-based standard errors (column se). \n", sep = "")
            }
            if(!is.null(object$df)){
                cat("  Degrees of freedom were computed using a Satterthwaite approximation (column df). \n", sep = "")
            }
            if("lower" %in% columns && "upper" %in% columns){
                cat("  The columns lower and upper indicate a ",100*level,"% confidence interval for each coefficient.\n", sep = "")
            }else if("lower" %in% columns){
                cat("  The column lower indicates a ",100*level,"% confidence interval for each coefficient.\n", sep = "")
            }else if("upper" %in% columns){
                cat("  The column upper indicate a ",100*level,"% confidence interval for each coefficient.\n", sep = "")
            }

            if(("lower" %in% columns) || ("upper" %in% columns) || (!is.null(object$df)) || (robust && "se" %in% columns)){
                cat("\n")
            }
        }
    }else{
        table.mean <- NULL
    }
    
    ## ** export
    return(invisible(list(correlation = table.cor,
                          variance = table.var,
                          sd = table.sd,
                          mean = table.mean)))
}

## * summary.Wald_lmm (documentation)
##' @title Summary of Testing for a Linear Mixed Models
##' @description Estimates, p-values, and confidence intevals for linear hypothesis testing, possibly adjusted for multiple comparisons.
##' 
##' @param object an \code{anova_lmm} object, output of \code{anova}.
##' @param print [logical] should the output be printed in the console.
##' Can be a vector of length 2 where the first element refer to the global tests and the second to the individual tests.
##' @param seed [integer] value that will be set before adjustment for multiple comparisons to ensure reproducible results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param columns [character vector] Columns to be displayed for each null hypothesis.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param digits [interger] number of digits used to display estimates.
##' @param digits.p.value [interger] number of digits used to display p-values.
##' @param ... arguments \code{method}, \code{level}, and \code{backtransform} passed to \code{\link{confint_anova_lmm}}
##'
##'
##' @details By default adjustment for multiple comparisons via a single step max-test adjustment,
##' either using the multcomp package (equal degrees of freedom) or the copula package (unequal degrees of freedom).
 
## * summary.Wald_lmm (code)
##' @export
summary.Wald_lmm <- function(object, print = TRUE, seed = NULL, columns = NULL,
                             digits = max(3L, getOption("digits") - 2L),
                             digits.p.value = max(3L, getOption("digits") - 2L),
                             ...){

    if(!is.null(seed)){
        old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit( assign(".Random.seed", old.seed, envir = .GlobalEnv, inherits = FALSE) )
        set.seed(seed)
    }

    if(length(print)==1){
        print.indiv <- print
        print.global <- print
    }else if(length(print)>2){
        stop("Argument \'print\' should have length at most 2. \n",
             "The first element refering to global test and the second to individual hypotheses. \n")
    }else{
        print.global <- print[1]
        print.indiv <- print[2]
    }
    options <- LMMstar.options()
    valid.columns <- c("null","estimate","se","statistic","df","lower","upper","p.value","partial.r","")
    if(identical(columns,"all")){
        columns.global <- setdiff(valid.columns, c("estimate", "se", "lower", "upper"))
        columns.indiv <- valid.columns
    }else  if(is.null(columns)){
        columns.indiv <- options$columns.anova
        columns.global <- union("statistic", setdiff(options$columns.anova, c("estimate", "se", "lower", "upper")))
    }else{
        columns.indiv <- tolower(columns)
        if(any(columns.indiv %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns.indiv[columns.indiv %in% valid.columns == FALSE], collapse ="\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns.indiv), collapse ="\" \""),"\".\n")
        }
        columns.global <- setdiff(columns.indiv, c("estimate", "se", "lower", "upper"))
    }
    if("df" %in% columns.global){
        index.df <- which(columns.global == "df")
        if(index.df == 1){
            columns.global <- c("df.num", "df.denom", columns.global[(index.df+1):length(columns.global)])
        }else if(index.df == length(columns.global)){
            columns.global <- c(columns.global[1:(index.df-1)], "df.num", "df.denom")
        }else{
            columns.global <- c(columns.global[1:(index.df-1)], "df.num", "df.denom", columns.global[(index.df+1):length(columns.global)])
        }
    }
    columns.global <- gsub("^partial.r$","partial.r2", columns.global)
    object.df <- object$args$df
    object.robust <- object$args$robust
    object.ci <- object$args$ci

    table.multivariate <- object$multivariate
    table.multivariate$type.original <- object$args$type[[1]]
    type <- unique(table.multivariate$type.original)

    if(object.ci){
        table.univariate <- confint(object, columns = union(c("type","test","method"),setdiff(columns.indiv,"")), ...)
        if(NROW(table.multivariate)==1){
            table.univariate$type.original <- table.multivariate$type.original
        }else{
            typetest2type.original <- stats::setNames(table.multivariate$type.original,paste(table.multivariate$type,table.multivariate$test,sep="|"))
            table.univariate$type.original <- typetest2type.original[paste(table.univariate$type,table.univariate$test,sep="|")]
        }
        univariate.method <- attr(table.univariate,"method")
    }

    for(iType in type){ ## iType <- type[1]

        ## ** Type of test
        if(iType == "all" && (print.global>0.5 || print.indiv>0.5)){
            cat("\n\t", "|| User-specified linear hypotheses || \n", sep="")
            ## print(object$call)
        }else if(print.global>0.5 || print.indiv>0.5){
            cat("\n\t","|| ",iType," coefficients || \n", sep="")
        }

        ## ** Global tests
        object.print <- table.multivariate[table.multivariate$type.original==iType,,drop=FALSE]
        object.print <- cbind(object.print,
                              " " = stats::symnum(object.print$p.value, corr = FALSE, na = FALSE, 
                                                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                  symbols = c("***", "**", "*", ".", " "))
                              )
        names(object.print)[NCOL(object.print)] <- ""
        object.print$p.value <- as.character(signif(object.print$p.value, digits = digits.p.value))
        rownames(object.print) <- object.print$test

            if(print.global){

                txt.test <- "Multivariate Wald test (global null hypothesis)"
                if(all(is.na(object.print$statistic)) && !is.null(attr(object.print$statistic,"error"))){
                    cat("\n - ",txt.test,": ",attr(object.print$statistic,"error"),"\n", sep="")
                }else{
                    cat("\n - ",txt.test,"\n", sep="")
                    object.print <-  object.print[,names(object.print) %in% columns.global,drop=FALSE]
                    if(object.df){
                        names(object.print) <- gsub("^statistic","F-statistic",names(object.print))
                    }else{
                        names(object.print) <- gsub("^statistic","chi2-statistic",names(object.print))
                    }
                    print(object.print, digits = digits, row.names = any(rownames(object.print)!="1"))
                }
            }

            ## ** Individual specific tests
            if(object.ci){
                                
                if(print.indiv){

                    cat("\n - Univariate Wald test (individual null hypotheses) \n", sep="")
                    object.print <- table.univariate[table.univariate$type.original==iType,,drop=FALSE]
                    n.hypoPerTest <- table(paste(object.print$type,object.print$test,sep="|"))
                    if(length(univariate.method)>1){
                        warning("Different methods have been used to adjust for multiple comparisons - text describing the adjustment will not be accurate.")
                    }
                    object.print <- object.print[,names(object.print) %in% columns.indiv,drop=FALSE]
                    if(object.df){
                        names(object.print) <- gsub("^statistic","t-statistic",names(object.print))
                    }else{
                        names(object.print) <- gsub("^statistic","z-statistic",names(object.print))
                    }

                    stats::printCoefmat(object.print, digits = digits,
                                        has.Pvalue = "p.value" %in% columns.indiv,
                                        P.values = "p.value" %in% columns.indiv,
                                        eps.Pvalue = 10^{-digits.p.value},
                                        signif.legend = FALSE)
                }
            }
        }


        ## ** Legend
        if(print.global || print.indiv){

            cat("---\n")
            if("" %in% columns.global || "" %in% columns.indiv){
                cat("Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n")
            }
            if(object.robust && "se" %in% columns.indiv){
                cat("  Standard errors: robust\n")
            }else if("se" %in% columns.indiv){
                cat("  Standard errors: model-based\n")
            }

            if("back-transform" %in% names(attributes(table.univariate))){
                message.backtransform <- attr(table.univariate,"back-transform")[!is.na(attr(table.univariate,"back-transform")$FUN),,drop=FALSE]

                if(any(message.backtransform[,setdiff(names(message.backtransform), "FUN")] == FALSE)){
                    warning("Could not back-transform everything.\n")
                }
                if(NROW(object.print)==1){
                    short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
                    txt <- unique(short2text[intersect(names(short2text),intersect(names(object.print),names(message.backtransform)))])
                }else{
                    short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
                    txt <- unique(short2text[intersect(names(short2text),intersect(names(object.print),names(message.backtransform)))])
                }
                substr(txt[1], 1, 1) <- toupper(substr(txt[1], 1, 1))
                cat("  ",paste(txt,collapse = ", ")," have been back-transformed",sep="")
                if(print>0.5){
                    cat(" (",paste0(paste(rownames(message.backtransform),collapse = "/")," parameters with ",paste(message.backtransform$FUN,collapse="/")),"). \n", sep ="")
                }
            }else if(any(!is.na(c(object$args$transform.sigma,object$args$transform.k,object$args$transform.rho)))){
                vec.transform <- na.omit(c(sigma=object$args$transform.sigma,k=object$args$transform.k,rho=object$args$transform.rho))
                if(NROW(object.print)==1){
                    short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
                    txt <- unique(short2text[intersect(names(short2text),names(object.print))])
                }else{
                    short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
                    txt <- unique(short2text[intersect(names(short2text),names(object.print))])
                }
                substr(txt[1], 1, 1) <- toupper(substr(txt[1], 1, 1))
                cat("  ",paste(txt,collapse = ", ")," have been transformed",sep="")
                if(print>0.5){
                    cat(" (",paste0(paste(names(vec.transform),collapse = "/")," parameters with ",paste(vec.transform,collapse="/")),"). \n", sep ="")
                }
            }

            if(any(c("p.value", "lower", "upper") %in% columns.indiv)){
                if("p.value" %in% columns.indiv == FALSE){
                    txt.cip <- "P-values"
                }else if("lower" %in% columns.indiv == FALSE && "upper" %in% columns.indiv == FALSE){
                    txt.cip <- "CIs"
                }else{
                    txt.cip <- "CIs/p-values"
                }
                if(univariate.method[1] == "none"){ ## always only one hypothesis in each global test
                    ## cat("  (",txt.cip," not adjusted for multiple comparisons) \n", sep="")
                }else if(length(n.hypoPerTest)==1){ ## only one global test
                    if(univariate.method[1] == "bonferroni"){
                        cat("  (",txt.cip," adjusted for multiple comparisons -- Bonferroni)\n", sep="")
                    }else if(univariate.method[1] %in% c("single-step", "single-step2")){
                        cat("  (",txt.cip," adjusted for multiple comparisons -- max-test adjustment)\n", sep="")
                    }else{
                        cat(paste0("  (",txt.cip," adjusted for multiple comparisons -- ",univariate.method[1],")\n", sep=""),sep="")
                    }
                }else{
                    if(univariate.method[1] == "bonferroni"){
                        cat("  (",txt.cip," adjusted for multiple comparisons within each global test -- bonferroni) \n", sep="")
                    }else if(univariate.method[1] %in% c("single-step","single-step2")){
                        cat("  (",txt.cip," adjusted for multiple comparisons within each global test -- max-test adjustment) \n", sep="")
                    }else{
                        cat(paste0("(",txt.cip," adjusted for multiple comparisons within each global test -- ",univariate.method[1],") \n", sep=""),sep="")
                    }
                }

                if(univariate.method[1] == "single-step"){
                    if(first){cat("\n");first <- FALSE}
                    error <- max(c(0,abs(attr(table.univariate,"error")[table.multivariate$type.original==iType])), na.rm = TRUE)
                    if(error > 1e-12){
                        txt.error <- paste0("Error when computing the adjusted ",txt.cip," by numerical integration: ", signif(error, digits = 5))
                        if(!is.null(seed)){
                            txt.error <- paste0(txt.error," (seed ",seed,")")
                        }
                        cat(txt.error,".\n",sep="")
                    }
                }else if(univariate.method[1] == "single-step2"){
                    if(first){cat("\n");first <- FALSE}
                    txt.sample <- paste("Adjusted ",txt.cip," computed using ",attr(table.univariate,"n.sample")," samples", sep = "")
                    if(!is.null(seed)){
                        txt.sample <- paste0(txt.error," (seed ",seed,")")
                    }
                    cat(txt.sample,".\n")                            
                }
            }
            cat("\n")
        }
        
   
}


## * summary.LRT_lmm
summary.LRT_lmm <- function(object, digits = 3, columns = NULL, ...){

    ## ** normalize input
    if(length(digits)==1){
        digits <- rep(digits,2)
    }
    valid.columns <- c("null","logLikH0","logLikH1","statistic","df","p.value","")

    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(is.null(columns) || identical(columns,"all")){
        columns <- setdiff(valid.columns,"null")
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

    out <- as.data.frame(object)
    out[names(out)[names(out) %in% columns == FALSE]] <- NULL
    rownames(out) <- ""
    toPrint <- capture.output(stats::printCoefmat(out, digits = digits[1],
                                                  has.Pvalue = "p.value" %in% columns,
                                                  P.values = "p.value" %in% columns,
                                                  eps.Pvalue = 10^{-digits[2]},
                                                  signif.legend = FALSE))
    sapply(toPrint, function(iTxt){cat("  ",iTxt,"\n",sep="")})
    cat(" ",paste(rep("-", max(nchar(toPrint))),collapse=""),"\n")
    if("" %in% columns){
        cat("  Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n")
    }
    cat("\n")

    ## ** export
    return(invisible(out))
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
##' @param ... other arguments are passed to \code{\link{summary.anova_lmm}}.
##'

## * summary.mlmm (code)
##' @export
summary.mlmm <- function(object, digits = 3, method = NULL, print = NULL, hide.data = FALSE, hide.fit = FALSE, ...){

    if(is.null(method)){
        method <- "none"
    }
    if(is.null(print)){
        print <- c(FALSE,TRUE)
    }


    ## extract models
    ls.model <- object$model
    method.fit <- object$args$method.fit
    optimizer <- ls.model[[1]]$opt$name
    logLik <- sapply(ls.model, logLik)
    cv <- sapply(ls.model, function(iM){iM$opt$cv})
    n.iter <- sapply(ls.model, function(iM){iM$opt$n.iter})
    param.value <- coef(object, effects = "all")

    nparam.mu  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="mu")})
    nparam.sigma  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="sigma")})
    nparam.k  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="k")})
    nparam.rho  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="rho")})

    M.nobs <- do.call(rbind,lapply(ls.model, stats::nobs))
    vec.nobs <- apply(M.nobs,2,paste,collapse = ", ")
    nobsByCluster <- lapply(ls.model,function(iM){iM$design$cluster$nobs})
    n.cluster.original <- sapply(ls.model,function(iM){iM$design$cluster$n})
    n.cluster.design <- sapply(ls.model,function(iM){iM$design$cluster$n})

    call <- attr(object,"call")
    
    ## ** welcome message
    if(any(print)){
        cat("	Linear Mixed Models stratified according to \"",call$by,"\" \n\n",sep="")
    }

    ## ** data message    
    if(!hide.data){
        cat("Dataset:", deparse(call$data), "\n")
        cat("Strata : \"",paste(names(ls.model),collapse = "\", \""),"\"\n\n",sep="")
        if(any(M.nobs[,"missing"]>0)){
            if(any(n.cluster.original-n.cluster.design>0)){
                cat("  - ", vec.nobs["cluster"], " clusters were analyzed, ",paste(n.cluster.original-n.cluster.design,collapse=", ")," were excluded because of missing values \n" , sep = "")
            }else{
                cat("  - ", vec.nobs["cluster"], " clusters \n" , sep = "")
            }
            cat("  - ", paste(sapply(nobsByCluster,sum), collapse = ", "), " observations were analyzed, ",vec.nobs["missing"]," were excluded because of missing values \n",  sep = "")
        }else{
            cat("  - ", vec.nobs["cluster"], " clusters \n" , sep = "")
            cat("  - ", paste(sapply(nobsByCluster,sum), collapse = ", "), " observations \n",  sep = "")
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
        if(optimizer!="gls"){
            cat("  - convergence: ",paste(cv>0,collapse = ", ")," (",paste(n.iter,collapse = ", ")," iterations) \n", sep = "")
        }
        cat(" \n")
    }


    ## ** extract test
    if(any(print)){
        cat("Statistical inference \n")
        if("columns" %in% names(list(...))){
            out <- summary.anova_lmm(object, method = method, print = print/2, ...)
        }else{
            options <- LMMstar.options()
            out <- summary.anova_lmm(object, method = method, columns = options$columns.summary, print = print/2, ...)
        }
    }

    ## ** export
    return(invisible(out))
}

######################################################################
### summary.R ends here
 
