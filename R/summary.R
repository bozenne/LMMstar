 ### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: aug  2 2024 (10:42) 
##           By: Brice Ozenne
##     Update #: 1729
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.effect_lmm (code)
##' @export
summary.effect_lmm <- function(object, columns = NULL, print = TRUE, ...){

    ## ** normalize user input
    if("print" %in% names(match.call())==FALSE && all(is.na(object$multivariate$df.num))){
        print <- c(0,1.5)        
    }
    
    if("columns" %in% names(match.call())==FALSE){
        if(object$args$type=="outcome"){
            object$univariate$p.value <- NULL
            columns <- c("n","estimate","se","df","lower","upper")
        }else{
            columns <- c("n","estimate","se","df","lower","upper","p.value")
        }
    }

    ## ** prepare
    outcome.txt <- switch(object$args$type,
                          "outcome" = "outcome",
                          "change" = "change in outcome",
                          "auc" = "area under the outcome curve",
                          "auc-b" = "area under the outcome curve above baseline")
    contrast.txt <- switch(object$args$effect,
                           "identity" = "Average",
                           "difference" = "Difference in average")

    if(object$args$type=="change" & length(object$args$ref.repetition[[1]]==1) & max(print) > 0.7){
        if(length(attr(object$args$time,"original"))>1){
            vec.ref <- unlist(attr(object$args$ref.repetition,"original"))
            reference.txt <- paste0("  Reference repetition: ",paste(paste0(names(vec.ref),"=",vec.ref), collapse=", "),"\n",sep="")
        }else{
            reference.txt <- paste0("  Reference repetition: ",object$args$time,"=",object$args$ref.repetition[[1]],"\n",sep="")
        }            
    }else{
        reference.txt <- NULL
    }

    ## ** display
    if(max(print)>0){
        ## Note: print calls summary with argument print 0.5
        if(is.null(object$args$variable) || max(print) < 0.7){
            cat("\t\t",contrast.txt," counterfactual ",outcome.txt,"\n\n", sep = "")
        }else{
            cat("\t\t",contrast.txt," counterfactual ",outcome.txt,"\n\t\t w.r.t \'",object$args$variable,"\' values \n\n", sep = "")
        }
    }

    ## ** contrast
    contrast <- coef(object, type = "contrast")
    if(max(print)>=1){
        cat("\tPlanned contrast: \n")
        contrast.print <- contrast
        rownames(contrast.print) <- paste0("   ",rownames(contrast.print))
        print(contrast.print[,colSums(abs(contrast))>0,drop=FALSE])
        cat("\n")
    }



    if(print[1]==0){
        if(NROW(object$univariate)==1){
            cat("\tUnivariate Wald test: \n")
        }else if(NROW(object$univariate)>1){
            cat("\tUnivariate Wald tests: \n")
        }
        out <- summary.Wald_lmm(object, print = c(0,0.4), columns = columns, ...)
        if(!is.null(reference.txt)){
            cat(reference.txt)
        }
        cat("\n")
    }else{
        if(!is.null(reference.txt)){
            attr(print,"message") <- reference
        }
        out <- summary.Wald_lmm(object, print = print, columns = columns, ...)
    }

    ## ** export
    return(invisible(contrast))
}

## * summary.lmm (documentation)
##' @title Summary Output for a Linear Mixed Model
##' @description Summary output for a linear mixed model fitted with \code{lmm}.
##'
##' @param object [lmm] output of the \code{lmm} function.
##' @param level [numeric,0-1] confidence level for the confidence intervals.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' Can also be \code{2} compute the degrees of freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
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
summary.lmm <- function(object, level = 0.95, robust = FALSE, df = NULL,
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
    options <- LMMstar.options()

    n.strata <- object$strata$n
    U.strata <- object$strata$levels
    
    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** df
    if(is.null(df)){
        df <- !is.null(object$df)
    }
    
    ## *** columns
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
    
    ## *** type.cor
    if(!is.null(type.cor)){
        type.cor <- match.arg(type.cor, c("matrix","param"))
    }

    ## *** hide
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
                              df = df,
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
        cat("  ","Null hypothesis: ",table$null,"\n\n",sep="")
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
        name.param <- unique(unlist(object$univariate$parameter))
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



## * summary.resample
##' @export
summary.resample <- function(object, digits = 3, ...){

    object.type <- object$args$type
    n.sample <- object$args$n.sample
    
    ## ** check cv
    n.cv <- sum(object$cv)

    ## ** compute CI/pvalue
    table.object <- model.tables(object, ...)

    txt.type <- switch(object.type,
                       "perm-var" = "permutation",
                       "perm-res" = "permutation",
                       "boot" = "bootstrap")
    txt.method <- switch(object.type,
                         "perm-var" = paste0("p-value computed using ", attr(table.object,"method")," method"),
                         "perm-res" = paste0("p-value computed using ", attr(table.object,"method")," method"),
                         "boot" = paste0("lower,upper,p-value computed using ", attr(table.object,"method")," method"))

    ## ** display
    if(object.type == "perm-var"){
        cat("\tPermutation test (",n.sample," samples) \n\n", sep = "")
    }else if(object.type == "perm-res"){
        cat("\tResidual Permutation test (",n.sample," samples) \n\n", sep = "")
    }else if(object.type == "boot"){
        cat("\tNon-parametric bootstrap (",n.sample," samples)\n\n", sep = "")
    }

    rownames(table.object) <- paste0("   ",rownames(table.object))
    base::print.data.frame(table.object, digits = digits)

    cat("   ",rep("-",ncharTable(table.object, digits = digits)),"\n",sep="")
    cat("  ",txt.method,"\n")
    if(n.sample!=n.cv){
        cat(paste0("   Based on ",n.cv," ",txt.type," samples - ",round((1-n.cv/n.sample)*100, digits = digits),"% failed\n"))
    }else{
        cat(paste0("   All ",txt.type," samples were successful\n"))
    }


    cat("\n")
    return(invisible(NULL))
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

    ## ** extract from object
    options <- LMMstar.options()
    pool.method <- options$pool.method
    type.information <- object$object$type.information
    df <- object$args$df
    robust <- object$args$robust

    transform.sigma <- object$args$transform.sigma
    transform.k <- object$args$transform.k
    transform.rho <- object$args$transform.rho

    ## ** normalize input
    ## *** print
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

    ## *** columns
    valid.columns <- c("null","type","estimate","se","statistic","df","quantile","lower","upper","p.value","")
    if(identical(columns,"all")){
        columns.multivariate <- valid.columns
        columns.univariate <- valid.columns
    }else{
        if(is.null(columns) || !is.null(names(columns))){
            columns.univariate <- options$columns.anova
            columns.multivariate <- union("statistic", options$columns.anova)
        }else{
            columns <- tolower(columns)
            if(any(columns %in% valid.columns == FALSE) && any(columns %in% names(object$univariate) == FALSE)){
                stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse ="\" \""),"\" for argument \'columns\'. \n",
                     "Valid values: \"",paste(setdiff(valid.columns, columns), collapse ="\" \""),"\".\n")
            }
            if(!is.null(columns) && any(names(columns) %in% c("add","remove") == FALSE)){
                stop("Incorrect names for argument \'columns\': should be \"add\" or \"remove\". \n")
            }
        }

        if(!is.null(columns)){
            if(is.null(names(columns))){
                columns.univariate <- columns
                columns.multivariate <- columns
            }else{
                columns.univariate <- setdiff(union(columns.univariate, unname(columns[names(columns)=="add"])),
                                              unname(columns[names(columns)=="remove"]))
                columns.multivariate <- setdiff(union(columns.multivariate, unname(columns[names(columns)=="add"])),
                                                unname(columns[names(columns)=="remove"]))
            }
        }
    }
    columns.multivariate <- setdiff(columns.multivariate, c("estimate", "se", "quantile", "lower", "upper"))

    if(length(columns.univariate)==0 || !object$args$univariate){
        print.univariate <- FALSE
    }
    if(length(columns.multivariate)==0 || !object$args$multivariate){
        print.multivariate <- FALSE
    }
    if(print.multivariate == FALSE && print.univariate == FALSE){
        stop("Nothing to print. \n")
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

    ## ** ensure reproducibility
    if(!is.null(seed)){
        old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit( assign(".Random.seed", old.seed, envir = .GlobalEnv, inherits = FALSE) )
        set.seed(seed)
    }

    ## ** extract information

    ## *** multivariate tests
    if(print.multivariate>0){

        table.multivariate <- object$multivariate[,setdiff(columns.multivariate,c("type","")),drop=FALSE]
        
        ## add type
        if(object$args$type=="auto"){            
            if("type" %in% columns.multivariate == FALSE){
                nchar.type <- nchar(object$multivariate$term)
                maxchar.type <- max(nchar.type)
                vec.white <- sapply(maxchar.type-nchar.type, function(iN){paste(rep(" ", iN), collapse = "")})
                name.type <- sapply(object$multivariate$type, switch,
                                    "mu" = "mean",
                                    "k" = "variance",
                                    "rho" = "correlation")
                rowname.multivariate <- paste0(ifelse(duplicated(name.type),"",paste0(name.type,": ")),paste(object$multivariate$term,vec.white,sep=""))
                
                nchar2.type <- nchar(rowname.multivariate)
                maxchar2.type <- max(nchar2.type)
                vec2.white <- sapply(maxchar2.type-nchar2.type, function(iN){paste(rep(" ", iN), collapse = "")})
                rownames(table.multivariate) <- paste(vec2.white,rowname.multivariate,sep="")                                
            }else{
                rownames(table.multivariate) <- object$multivariate$term
            }
        }else{            
            rownames(table.multivariate) <- "all"
        }

        ## restaure attributes (i.e. message)
        attr(table.multivariate,"se") <- attr(object$multivariate,"se")
        attr(table.multivariate,"df") <- attr(object$multivariate,"df")
        
        if(print.multivariate>0.5){
            if(NROW(table.multivariate)==1){
                cat("\t\tMultivariate Wald test \n\n")
            }else{
                cat("\t\tMultivariate Wald tests \n\n")
            }
        }

        .printStatTable(table = table.multivariate, robust = robust, df = df, level = NULL, type.information = type.information,
                        method.p.adjust = NULL,
                        backtransform = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                        columns = setdiff(columns.multivariate,"type"), col.df = c("df.num","df.denom"), name.statistic = c("Chi2-statistic","F-statistic"),
                        digits = digits, digits.df = c(df.num = 0, df.denom = digits.df), digits.p.value = digits.p.value,
                        decoration = legend, legend = legend)
        if(any(!is.na(object$multivariate$message))){
            cat(paste(unique(stats::na.omit(object$multivariate$message)), collapse = "\n"))
        }
        cat("\n")
    }

    ## *** univariate tests
    if(print.univariate>0){
        
        table.univariate <- confint(object, columns = union(setdiff(columns.univariate,""),c("type","term","name")), ...)
        if(is.null(columns) && all(is.na(table.univariate$lower)) && all(is.na(table.univariate$upper))){
            columns.univariate <- setdiff(columns.univariate, c("lower","upper"))
        }
        error <- attr(table.univariate,"error")
        n.sample <- attr(table.univariate,"n.sample")
        method.p.adjust <- attr(table.univariate,"method")
        level <- attr(table.univariate,"level")
        backtransform <- attr(table.univariate,"backtransform")

        ## group of test relative to which multiple comparison adjustment is performed
        if(object$args$type=="auto" && length(unique(table.univariate$name))>1){
            if(length(unique(table.univariate$type))==1){
                factor.p.adjust <- "covariate"
            }else if(length(unique(table.univariate$term))==1){
                factor.p.adjust <- "parameter type"
            }else{
                factor.p.adjust <- "covariate and parameter type"
            }
        }else{
            factor.p.adjust <- NULL
        }

        ## type of parameter: mean/variance/correlation
        if(object$args$type=="auto" && "type" %in% columns.univariate == FALSE){
            nchar.type <- nchar(rownames(table.univariate))
            maxchar.type <- max(nchar.type)
            vec.white <- sapply(maxchar.type-nchar.type, function(iN){paste(rep(" ", iN), collapse = "")})
            name.type <- sapply(table.univariate$type, switch,
                                "mu" = "mean",
                                "k" = "variance",
                                "rho" = "correlation")
            rownames.univariate <- paste0(ifelse(duplicated(name.type),"",paste0(name.type,": ")),paste(rownames(table.univariate),vec.white,sep=""))

            nchar2.type <- nchar(rownames.univariate)
            maxchar2.type <- max(nchar2.type)
            vec2.white <- sapply(maxchar2.type-nchar2.type, function(iN){paste(rep(" ", iN), collapse = "")})
            rownames(table.univariate) <- paste(vec2.white,rownames.univariate,sep="")                                
        }

        if(print.univariate>0.5){
            if(NROW(table.univariate)==1){
                cat("\t\tUnivariate Wald test \n\n")
            }else{
                cat("\t\tUnivariate Wald tests \n\n")
            }
        }
        if(any(c("pool.gls","pool.gls1") %in% method.p.adjust)){
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
                        method.p.adjust = method.p.adjust, factor.p.adjust = factor.p.adjust, error.p.adjust = error, pool.method = pool.method, seed = seed, n.sample = n.sample,
                        backtransform = backtransform, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                        columns = setdiff(columns.univariate,"type"), col.df = c("df"), name.statistic = c("t-statistic","z-statistic"),
                        digits = digits, digits.df = digits.df, digits.p.value = digits.p.value2,
                        decoration = legend, legend = legend)

        if(!is.null(attr(print,"message"))){
            cat(attr(print,"message"))
        }
        if(print.univariate>0.5){
            cat("\n")
        }


    }

    ## ** export
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
##' @param pool.method [character vector] adjustment method for multiple comparisons consisting in pooling.
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
                            method.p.adjust = NULL, factor.p.adjust, error.p.adjust, pool.method = NULL, seed, n.sample,
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
    columns.num <- setdiff(names(table)[sapply(table,is.numeric)],c("p.value",col.df))
    for(iCol in columns.num){ 
        table.print[[iCol]] <- as.character(round(table[[iCol]], digits = digits))
        iIndex.NNA <- which(!is.na(table.print[[iCol]]))
        if(any(table[iIndex.NNA,iCol]!=0 & table.print[iIndex.NNA,iCol] == "0")){
            iIndex.0 <- iIndex.NNA[which(table[iIndex.NNA,iCol]!=0 & table.print[iIndex.NNA,iCol] == "0")]
            table.print[iIndex.0,iCol] <- paste0("<0.",paste0(rep(0,digits-1),collapse=""),"1")
        }
    }
    if("p.value" %in% names(table.print)){
        table.print$p.value <- as.character(format.pval(table.print$p.value, digits = digits.p.value[1], eps = digits.p.value[2]))
    }
    if(!is.null(col.df) && all(col.df %in% columns)){
        if(identical(col.df,"df")){
            table.print$df <- as.character(round(table.print$df, digits = digits.df))
        }else{
            table.print[[col.df[1]]] <- paste0("(",apply(do.call(cbind,lapply(col.df, function(iCol){ ## iCol <- col.df[2]
                if(all(is.na(table.print[[iCol]]))){
                    iDF <- as.character(table.print[[iCol]])
                }else{
                    iDF <- formatC(table.print[[iCol]], digits = digits.df[[iCol]], format = "f")
                }
                nchar.iDF <- nchar(iDF, allowNA = TRUE, keepNA = FALSE)
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
            cat(space,"  :  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1.\n",sep="")
        }

        ## *** type of degree of freeedom 
        if(identical(df,TRUE) && "df" %in% columns){
            if(!is.null(attr(table,"df"))){
                cat(space,"df: ",attr(table,"df"),". \n", sep = "")
            }else if(robust %in% 0:1){
                cat(space,"df: Satterthwaite approximation w.r.t. model-based se. \n", sep = "")
            }else if(robust == 2){
                cat(space,"df: Satterthwaite approximation w.r.t. robust se. \n", sep = "")
            }
        }

        ## *** type of standard error 
        if("se" %in% columns){
            if(robust>0){
                txt.se <- paste0("Robust based on the ",type.information," information")
            }else if(robust==0){
                txt.se <- paste0("Model-based based on the ",type.information," information")
            }
            if(!is.null(attr(table,"se"))){
                txt.se <- paste0(txt.se,attr(table,"se"))
            }
            cat(space,"se: ",txt.se,". \n", sep = "")
        }

        ## *** estimate
        if(length(method.p.adjust)==1 && method.p.adjust %in% pool.method && "estimate" %in% columns){

            if(NROW(table)>1){
                txt.strata <- paste("within ",factor.p.adjust," ",sep = "")
            }else{
                txt.strata <- ""
            }
            txt.estimate <- switch(method.p.adjust,
                                   "average" = paste0("averaged estimates",txt.strata),
                                   "pool.fixse" = paste0("pooled estimates",txt.strata," using inverse variance weights"),
                                   "pool.se" = paste0("pooled estimates",txt.strata," using inverse variance weights"),
                                   "pool.gls" = paste0("pooled estimates",txt.strata," using GLS weights"),
                                   "pool.gls1"= paste0("pooled estimates",txt.strata," using constrained GLS weights"),
                                   "pool.rubin" = paste0("pooled estimates",txt.strata," using Rubin's rule")
                                   )
            cat(space,"estimate: ",txt.estimate)

        }

        ## *** Adjustment for multiple testing
        method.p.adjust2 <- setdiff(method.p.adjust, pool.method)

        if(!is.null(level) && !is.na(level)){
            display.cip <- intersect(c("lower","upper","p.value"),columns)
        }else{
            display.cip <- intersect(c("p.value"),columns)

        }


        if((df && name.statistic[2] == "F-statistic") || (!df && name.statistic[1] == "Chi2-statistic")){

            if("df.num" %in% names(table) == FALSE){
                ## nothing: cannot say what is going without the df.num column
            }else if(NROW(table)==1 && table$df.num==1){
                ## nothing: no need for adjustment
            }else if(NROW(table)==1 && table$df.num>1){
                cat(space,"Multiple testing adjustment: joint test.\n", sep = "")
            }else if(NROW(table)>1 && all(table$df.num==1)){
                cat(space,"No adjustment for multiple testing.\n", sep = "")
            }else if(NROW(table)>1 && any(table$df.num>1)){
                test.type <- grepl("mean:",rownames(table),fixed=TRUE) + grepl("variance:",rownames(table),fixed=TRUE) + grepl("correlation:",rownames(table),fixed=TRUE)
                if(sum(test.type)>1 && sum(test.type)==NROW(table)){
                    cat(space,"Multiple testing adjustment: within parameter type using joint test.\n", sep = "")
                }else if(sum(test.type)>1){
                    cat(space,"Multiple testing adjustment: within covariate and parameter type using joint test.\n", sep = "")
                }else{
                    cat(space,"Multiple testing adjustment: within covariate using joint test.\n", sep = "")
                }
            }
        }else if((is.null(method.p.adjust) || (length(method.p.adjust2)==1 && method.p.adjust2 == "none")) && length(display.cip)>0 && NROW(table)>1){
            ## no adjustment
            cat(space,"No adjustment for multiple testing.\n", sep = "")
        }else if(length(method.p.adjust2)>=1 && length(display.cip)>0 && NROW(table)>1){

            ## adjustment
            if(any(c("single-step", "single-step2") %in% method.p.adjust)){

                name.adjmethod <- "max test"

                txt.adjustment <- NULL ## paste(display.cip, collapse = "/")
                if(!is.null(seed)){
                    txt.adjustment <- paste(c(txt.adjustment,paste0("RNG seed ",seed)), collapse=", ")
                }
                if("single-step" %in% method.p.adjust && !is.null(error.p.adjust) && any(!is.na(error.p.adjust)) && any(abs(stats::na.omit(error.p.adjust))>1e-12)){
                    txt.adjustment <- paste(c(txt.adjustment, paste0("numerical intergration error ", signif(max(error.p.adjust, na.rm=TRUE), digits = digits.p.value[1]))), collapse=", ")
                }
                if("single-step2" %in% method.p.adjust){
                    txt.adjustment <- paste(c(txt.adjustment, paste0(n.sample," samples")), collapse=", ")
                }
                if(!is.null(txt.adjustment) && nchar(txt.adjustment)>0){
                    txt.adjustment <- paste0(" (",txt.adjustment,")")
                }
                
            }else{
                name.adjmethod <- method.p.adjust
                txt.adjustment <- ""
            }

            if(!is.null(factor.p.adjust) && nchar(factor.p.adjust)>0){
                txt.adjustment <- paste0("within ", factor.p.adjust," using ",name.adjmethod,txt.adjustment)
            }else{
                txt.adjustment <- paste0("using ",name.adjmethod,txt.adjustment)
            }

            if(!is.null(txt.adjustment)){
                cat(space,"Multiple testing adjustment: ",txt.adjustment,".\n",sep="")
            }

        }

        ## *** backtransformation
        if(!is.null(backtransform) && any(!is.na(backtransform$FUN))){

            vec.backtransform <- backtransform[!is.na(backtransform$FUN),]
         
            cat(space,"Back-transformation: ",paste0(paste(rownames(vec.backtransform),collapse = "/")," parameters with ",paste(vec.backtransform$FUN,collapse="/")),".\n",
                ## " (",paste(intersect(c("estimate","se","lower","upper"),columns),collapse = "/"),"). \n",
                sep ="")

        }

    }
    
    ## ** export
    return(invisible(NULL))

}


######################################################################
### summary.R ends here
 
   
