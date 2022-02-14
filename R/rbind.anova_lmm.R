### rbind.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: feb 14 2022 (12:20) 
##           By: Brice Ozenne
##     Update #: 37
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * rbind.anova_lmm
##' @title Linear Hypothesis Testing Across Linear Mixed Models
##' @description Linear hypothesis testing accross linear mixed model.
##'
##' @param model a \code{anova_lmm} object (output of \code{anova} applied to a \code{lmm} object)
##' @param ...  possibly other \code{anova_lmm} objects
##' @param sep [character] character used to separate the outcome and the covariate when naming the tests.
##'
##' @examples
##' ## simulate data
##' set.seed(10)
##' dL <- sampleRem(1e2, n.times = 3, format = "long")
##'
##' ## estimate mixed models
##' e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, data = dL)
##' e.lmm2 <- lmm(Y ~ X1+X8+X9, repetition = ~visit|id, data = dL)
##'
##' ## select null hypotheses
##' AAA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"))
##' BBB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"))
##'
##' ## combine
##' ZZZ <- rbind(AAA,BBB)
##' summary(ZZZ)
##' @export
rbind.anova_lmm <- function(model, ..., sep = ": "){
    default <- LMMstar.options()

    ## ** check user input
    dots <- list(...)
    if(any(sapply(dots,inherits,"anova_lmm")==FALSE)){
        stop("Extra arguments should inherit from anova_lmm. \n")
    }
    ls.object <- c(list(model),dots)
    if(any(sapply(ls.object,function(iO){"all" %in% names(iO)})==FALSE)){
        stop("All argument should correspond to user specified hypothesis, i.e. call anova with argument linfct. \n")
    }
    if(any(sapply(ls.object,function(iO){!is.null(attr(iO$all,"glht"))})==FALSE)){
        stop("All argument should contain a \"glht\" object, i.e. call anova with argument ci=TRUE. \n")
    }
    ls.glht <- lapply(ls.object, function(iO){attr(iO$all,"glht")[[1]]})
    if(any(sapply(ls.glht, is.null))){
        stop("Could not extract glht object. \n Make sure that argument \'ci\' is TRUE when calling anova. \n")
    }

    ## ** Extract elements from anova object
    ls.C <- lapply(ls.glht,"[[","linfct")
    ls.rhs <- lapply(ls.glht,"[[","rhs")
    ls.coef <- lapply(ls.glht,"[[","coef")
    ls.lmm <- lapply(ls.glht,"[[","model")
    ls.df <- lapply(ls.object,function(iO){attr(iO$all,"CI")[[1]]$df})
    ls.alternative <- lapply(ls.glht,"[[","alternative")
    if(any(unlist(ls.alternative) != ls.alternative[[1]])){
        stop("Element \'alternative\' should take the same value for all glht objects. \n")
    }
    ls.robust <- lapply(ls.glht,function(iO){attr(iO$vcov,"robust")})
    if(any(unlist(ls.robust) != ls.robust[[1]])){
        stop("Element \'robust\' should take the same value for all glht objects. \n")
    }
    robust <- unique(unlist(ls.robust))
    ls.transform.sigma <- lapply(ls.object[[1]],function(iO){if(is.null(iO$call$transform.sigma)){default$transform.sigma}else{iO$call$transform.sigma}})
    if(any(unlist(ls.transform.sigma) != ls.transform.sigma[[1]])){
        stop("Element \'transform.sigma\' should take the same value for all glht objects. \n")
    }
    transform.sigma <- unique(unlist(ls.transform.sigma))
    ls.transform.k <- lapply(ls.object[[1]],function(iO){if(is.null(iO$call$transform.k)){default$transform.k}else{iO$call$transform.k}})
    if(any(unlist(ls.transform.k) != ls.transform.k[[1]])){
        stop("Element \'transform.k\' should take the same value for all glht objects. \n")
    }
    transform.k <- unique(unlist(ls.transform.k))
    ls.transform.rho <- lapply(ls.object[[1]],function(iO){if(is.null(iO$call$transform.rho)){default$transform.rho}else{iO$call$transform.rho}})
    if(any(unlist(ls.transform.rho) != ls.transform.rho[[1]])){
        stop("Element \'transform.rho\' should take the same value for all glht objects. \n")
    }
    transform.rho <- unique(unlist(ls.transform.rho))
    ls.method.fit <- lapply(ls.lmm,"[[","method.fit")
    if(any(unlist(ls.method.fit) != ls.method.fit[[1]])){
        stop("Element \'ls.method.fit\' should take the same value for all glht objects. \n")
    }
    method.fit <- unique(unlist(ls.method.fit))
    vec.outcome <- unname(unlist(sapply(ls.lmm,"[[","outcome")))

    names(ls.C) <- vec.outcome
    names(ls.coef) <- vec.outcome
    names(ls.lmm) <- vec.outcome

    ## ** extract iid
    ls.iid <- lapply(ls.lmm, function(iO){
        iid(iO, effects = if(method.fit=="REML"){"mean"}else{"all"}, robust = robust,
            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        
    })
    names(ls.iid) <- vec.outcome

    ## ** build glht object
    out <- list(model = ls.lmm, robust = robust)

    out$linfct <- as.matrix(do.call(Matrix::bdiag, ls.C))
    rownames(out$linfct) <- unlist(lapply(1:length(vec.outcome), function(iO){paste0(vec.outcome[iO],sep,rownames(ls.C[[iO]]))}))
    colnames(out$linfct) <- unlist(lapply(1:length(vec.outcome), function(iO){paste0(vec.outcome[iO],sep,colnames(ls.C[[iO]]))}))

    out$rhs <- unlist(ls.rhs)
    out$coef <- unlist(lapply(1:length(vec.outcome), function(iO){stats::setNames(ls.coef[[iO]],paste0(vec.outcome[iO],sep,names(ls.coef[[iO]])))}))

    iIID <- do.call(cbind,ls.iid)
    out$vcov <- crossprod(iIID[rowSums(is.na(iIID))==0,,drop=FALSE])
    rownames(out$vcov) <- unlist(lapply(1:length(vec.outcome), function(iO){paste0(vec.outcome[iO],sep,colnames(ls.iid[[iO]]))}))
    colnames(out$vcov) <- unlist(lapply(1:length(vec.outcome), function(iO){paste0(vec.outcome[iO],sep,colnames(ls.iid[[iO]]))}))

    if(method.fit=="REML"){
        if(any(abs(out$linfct[,setdiff(colnames(out$linfct),colnames(out$vcov))]>1e-10))){
            stop("Cannot test covariance structure across models when using REML. \n",
                 "Consider setting argument \'method.fit\' to \"ML\" when calling lmm. \n")
        }
        rm.name <- unique(setdiff(colnames(out$linfct),colnames(out$vcov)))
        keep.col <- setdiff(1:NCOL(out$linfct),which(colnames(out$linfct) %in% rm.name))
        out$linfct <- out$linfct[,keep.col,drop=FALSE]
        out$coef <- out$coef[keep.col]
    }
    
    out$df <- unlist(ls.df)
    out$alternative <- ls.alternative[[1]]
    class(out) <- "glht"

    ## ** global test
    outSimp <- .simplifyContrast(out$linfct, out$rhs) ## remove extra lines
    iC.vcov.C_M1 <- try(solve(outSimp$C %*% out$vcov %*% t(outSimp$C)), silent = TRUE)
    if(inherits(iC.vcov.C_M1,"try-error")){
        iStat <- NA
        iDf <- c(outSimp$dim,Inf)
        attr(iStat,"error") <- "\n  Could not invert the covariance matrix for the proposed contrast."
    }else{
        iStat <- as.double(t(outSimp$C %*% out$coef - outSimp$rhs) %*% iC.vcov.C_M1 %*% (outSimp$C %*% out$coef - outSimp$rhs))/outSimp$dim 
        iDf <- c(outSimp$dim,Inf)
    }

    
    ## ** export
    out2 <- list(all = data.frame("null" = paste(rownames(out$linfct),collapse=", "),
                                  "statistic" = iStat,
                                  "df.num" = iDf[1],
                                  "df.denom" = iDf[2],
                                  "p.value" = 1 - stats::pf(iStat, df1 = iDf[1], df2 = iDf[2]),
                                  stringsAsFactors = FALSE),
                 call = lapply(ls.object,"[[","call"))

    attr(out2$all,"CI") <- list(data.frame(estimate = as.double(out$linfct %*% out$coef),
                                           se = sqrt(diag(out$linfct %*% out$vcov %*% t(out$linfct))),
                                           df = out$df,
                                           statistic = NA,
                                           lower = NA,
                                           upper = NA,
                                           null = out$rhs,
                                           p.value = NA,
                                           stringsAsFactors = FALSE))
    attr(out2$all,"CI")[[1]]$statistic <- (attr(out2$all,"CI")[[1]]$estimate-out$rhs)/attr(out2$all,"CI")[[1]]$se

    out$df <- ceiling(stats::median(out$df))
    attr(out2$all,"glht") <- list(out)
    
    attr(out2, "test") <- "Wald"
    attr(out2, "robust") <- robust
    class(out2) <- append("anova_lmm",class(out))
    return(out2)
}


##----------------------------------------------------------------------
### rbind.anova_lmm.R ends here
