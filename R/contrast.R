### contrast.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2024 (09:37) 
## Version: 
## Last-Updated: Jul 23 2024 (10:00) 
##           By: Brice Ozenne
##     Update #: 256
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simplifyContrast
##' @description remove contrasts making the contrast matrix singular.
##' Used in anova.lmm when doing multivariate Wald tests
##' 
##' @noRd
simplifyContrast <- function(object, rhs, tol = 1e-10, trace = TRUE){
    object.eigen <- eigen(tcrossprod(object))
    n.zero <- sum(abs(object.eigen$values) < tol)
    if(length(rhs)==1 && NROW(object)>1){
        rhs <- rep(rhs, NROW(object))
    }

    if(n.zero==0){return(list(C = object, rhs = rhs, dim = NROW(object), rm = 0))}
    
    keep.lines <- 1:NROW(object)
    for(iLine in NROW(object):1){ ## iLine <- 3
        iN.zero <- sum(abs(eigen(tcrossprod(object[setdiff(keep.lines,iLine),,drop=FALSE]))$values) < tol)
        if(iN.zero<n.zero){
            keep.lines <- setdiff(keep.lines,iLine)
            n.zero <- iN.zero
        }
        if(n.zero==0){
            
            if(trace){
                name.rm <- rownames(object)[-keep.lines]
                if(length(name.rm)==1){
                    message("Singular contrast matrix: contrast \"",name.rm,"\" has been removed. \n")
                }else{
                    message("Singular contrast matrix: contrasts \"",paste(name.rm,collapse= "\" \""),"\" have been removed. \n")
                }
            }
            return(list(C = object[keep.lines,,drop=FALSE], rhs = rhs[keep.lines], dim = length(keep.lines), rm = NROW(object)-length(keep.lines)))
        }
    }

    ## n.zero>0 so failure
    return(NULL)
}

## * equation2contrast (documentation)
##' @description Find the contrast matrix corresponding to an equation
##' Used by anova.lmm.
##' @noRd
##'
##' @param object [character vector] variable or equation to be converted to a contrast matrix.
##' @param name.coef [character vector] name of all model coefficients
##' @param X [list of matrix] design matrices for the mean and variance coefficients (cannot handle correlation design matrix).
##' @param name.arg [character] name of the argument to be display in error messages.
##' @param error.message [character] additional text for the error message.
##'
##' @examples
##' 
##' data(gastricbypassL, package = "LMMstar")
##' e.lmmUN <- lmm(glucagonAUC ~ visit + weight, repetition = ~time|id, data = gastricbypassL)
##'
##' #### only variable name ###
##' X <- model.matrix(e.lmmUN, simplify = FALSE)
##' equation2contrast("visit",
##'                   name.coef = names(coef(e.lmmUN, effects = "all")),
##'                   X = X$mean,
##'                   name.arg = "effects")
##' 
##' equation2contrast(c("visit","time"),
##'                   name.coef = names(coef(e.lmmUN, effects = "all")),
##'                   X = list(mean = X$mean, var = X$vcov$var$X),
##'                   name.arg = "effects")
##' 
##' #### only equation ###
##' equation2contrast(c("visit3-visit2=0"),
##'                   name.coef = names(coef(e.lmmUN, effects = "all")),
##'                   X = list(mean = X$mean, var = X$vcov$var$X),
##'                   name.arg = "effects")
##' 
##' equation2contrast(c(" visit3=0","5*rho(-13,-1) + rho(-13,1)=0","5.1 * rho(-13,-1) - log(3)*rho(-1,1)=0"),
##'                   name.coef = names(coef(e.lmmUN, effects = "all")),
##'                   X = list(mean = X$mean, var = X$vcov$var$X),
##'                   name.arg = "effects")


## * equation2contrast (code)
equation2contrast <- function(object, name.coef, X, 
                              name.arg){

    n.coef <- length(name.coef)

    ## ** normalize user input
    ## remove leading and ending white space
    ## remove breakline in case of too long character string
    object <- gsub("\n","",trimws(object, which = "both"), fixed = TRUE)

    rescue <- FALSE
    if(!is.null(attr(name.coef,"rescue"))){
        if(all(attr(name.coef,"rescue") %in% name.coef)){
            attr(name.coef,"rescue") <- NULL
        }else{
            attr(name.coef,"rescueOnly") <- setdiff(attr(name.coef,"rescue"),name.coef)
            attr(name.coef,"originalOnly") <- setdiff(name.coef,attr(name.coef,"rescue"))
        }
    }
    
    ## ** distinguish effects (X1+visit2=0) from variable name (visit)
    if(!is.character(object)){
        stop("Incorrect argument \'object\': should be a character vector. \n")
    }    
    if(any(duplicated(name.coef))){
        stop("Cannot generate a unique contrast matrix due to duplicated names among model coefficients. \n")
    }    
    test.equal <- grepl("=",object)
    equation <- object[test.equal]
    n.equation <- length(equation)
    variable <- object[!test.equal]
    n.variable <- length(variable)

    object.contrast <- stats::setNames(vector(mode = "list", length = length(object)), object)
    object.rhs <- stats::setNames(vector(mode = "list", length = length(object)), object)
    
    ## ** convert variables into contrast matrix
    if(n.variable>0){
        
        if(!is.list(X)){
            X <- list(mean = X)
        }else if(is.null(names(X)) || any(duplicated(names(X)))){
            stop("Incorrect argument \'X\': must be a named list with distinct names. \n")
        }else if(any(names(X) %in% c("mean","var","cor") == FALSE)){
            stop("Incorrect names for argument \'X\': should be \"mean\", \"var\", or \"cor\". \n")
        }

        ## extract all variables and add the X name (e.g. mean, variance, correlation)
        ls.X.variable <- lapply(1:length(X), function(iI){
            stats::setNames(attr(X[[iI]],"variable"),rep(names(X)[iI],length(attr(X[[iI]],"variable"))))
        })
        X.variable <- unlist(ls.X.variable)
        ## when duplicated variables favor the first contrast matrix, e.g. mean
        UX.variable <- X.variable[!duplicated(X.variable)] ## do not use unique as it removes the names
        
        if(all(variable %in% UX.variable)){
            
            variable2 <- UX.variable[UX.variable %in% variable]

            ls.contrastNull <- stats::setNames(lapply(1:length(variable2), function(iVar){ ## iVar <- 2
                iName <- names(variable2)[iVar] ## type of coefficient
                iNull <- switch(iName, "mean" = 0, "var" = 1, "cor" = 0) ## select null
                iX <- X[[iName]] ## select X
                iAssign <- which(attr(iX,"variable") == variable2[iVar]) ## position of the variable among all variables
                iLabel <- colnames(iX)[attr(iX,"assign") %in% iAssign] ## term(s) of the variable
                iN.label <- length(iLabel)
                iC <- matrix(0, nrow = iN.label, ncol = n.coef, dimnames = list(iLabel,name.coef)) ## create contrast matrix
                iC[iLabel,iLabel] <- diag(1, iN.label)
                return(list(contrast = iC, null = rep(iNull, iN.label)))
            }), variable2)

            object.contrast[] <- lapply(ls.contrastNull, "[[","contrast")
            object.rhs[] <- lapply(ls.contrastNull, "[[","null")

        }else{

            if(any(sapply(variable, function(iVar){sapply(name.coef, grepl, x = iVar)}))){
                test.variable <- sapply(variable, function(iVar){any(sapply(name.coef, grepl, x = iVar))})
                stop("Incorrect argument \'",name.arg,"\': e.g. when using \"",variable[test.variable][1],"\" the right-hand side of the equation is missing. \n",
                     "Consider using: \"",variable[test.variable][1],"=0\" instead. \n")
            }else if(!is.null(attr(name.coef, "rescue")) && any(sapply(variable, function(iVar){sapply(attr(name.coef,"rescueOnly"), grepl, x = iVar)}))){
                test.variable <- sapply(variable, function(iVar){any(sapply(attr(name.coef,"rescueOnly"), grepl, x = iVar))})
                stop("Incorrect argument \'",name.arg,"\': e.g. when using \"",variable[test.variable][1],"\" the right-hand side of the equation is missing. \n",
                     "Consider using: \"",variable[test.variable][1],"=0\" instead. \n")
            }else{
                stop("Incorrect variable value for argument \'",name.arg,"\': \"",paste(setdiff(variable,UX.variable), collapse = "\", \""),"\" \n",
                     "Possible variables: \"",paste(UX.variable, collapse = "\", \""),"\". \n")
            }
        }
    }

    ## ** convert simple equation (beta=0) into contrast matrix
    ## inspired from multcomp:::chrlinfct2matrix
    ## multcomp:::chrlinfct2matrix(equation, name.coef)
    if(n.equation>0){

        ## *** initialize contrast matrix
        equation.contrast <- matrix(0, nrow = n.equation, ncol = n.coef, dimnames = list(equation, name.coef))
        object.contrast[equation] <- lapply(equation, function(iE){equation.contrast[iE,,drop=FALSE]})

        ## *** split based on the equal sign
        splitEquation.equal <- strsplit(equation, split = "=", fixed = TRUE)
        if(any(lengths(splitEquation.equal)>2)){
            stop("Incorrect argument \'",name.arg,"\': should not contain more than one equal sign per element. \n",
                 "Problematic example: \"",equation[lengths(splitEquation.equal)>2][1],"\". \n")
        }       

        ## *** get right hand side
        object.rhs[equation] <- stats::setNames(lapply(splitEquation.equal, function(iE){as.numeric(iE[2])}),equation)

        ## *** get left hand side
        equation.lhs <- trimws(sapply(splitEquation.equal, "[", 1), which = "both")

        ## *** take care of the easy case where there is 'just' the coefficient
        index.coef <- which(equation.lhs %in% name.coef)
        ## check whether rescue coefficient is more approrpiate
        if(!is.null(attr(name.coef,"rescue")) && any(equation.lhs %in% attr(name.coef,"rescueOnly"))){
            if(all(equation.lhs[index.coef] %in% attr(name.coef,"rescue"))){
                rescue <- TRUE
                name.coef <- attr(name.coef,"rescue")
                index.coef <- which(equation.lhs %in% name.coef)
                equation.contrast <- matrix(0, nrow = n.equation, ncol = n.coef, dimnames = list(equation, name.coef))
                object.contrast[equation] <- lapply(equation, function(iE){equation.contrast[iE,,drop=FALSE]})
            }else{
                stop("Incorrect argument \'",name.arg,"\': contains untransformed and transformed parameter names. \n",
                     "Problematic example: \"",intersect(equation.lhs %in% attr(name.coef,"originalOnly"))[1],"\" and \"",intersect(equation.lhs %in% attr(name.coef,"rescueOnly"))[1],"\". \n")
            }
        }
        if(length(index.coef)>0){
            for(iC in index.coef){ ## iC <- 2
                object.contrast[[equation[iC]]][1,equation.lhs[iC]] <- 1
            }
        }
    }
    
    ## ** handle complex equation (5*beta-3*alpha=0) into a contrast matrix:
    if(n.equation>0){
        index.keep <- setdiff(1:n.equation, index.coef)
        equationComplex <- equation[index.keep]
        equationComplex.lhs <- equation.lhs[index.keep]
        n.equationComplex <- length(equationComplex)
    }else{
        n.equationComplex <- 0
    }
    if(n.equationComplex>0){
        ## *** standardize coefficient into letters paramAB, paramAC, ...
        stdname.coef <- paste0("param",sapply(strsplit(addLeading0(1:n.coef),split = "", fixed = TRUE), function(iVec){
            paste(LETTERS[as.numeric(iVec)+1], collapse = "")
        }))
        ## make sure it does no correspond to existing names by adding as many param as necessary
        if(any(grepl("param",name.coef,fixed = TRUE))){
            stdname.coef <- paste0(rep("param",max(lengths(strsplit(name.coef, split = "param",fixed = TRUE)))-1), stdname.coef)
        }
        names(stdname.coef) <- name.coef
        
        ## find which parameter appear in which equation     
        Mcoef.test <- do.call(rbind,lapply(name.coef, function(iCoef){sapply(equationComplex.lhs, grepl, pattern = iCoef, fixed = TRUE)}))
        rownames(Mcoef.test) <- name.coef
        ## check whether the 'transformed' parameter names are more appropriate
        if(!is.null(attr(name.coef,"rescue"))){
            Moriginal.test <- Mcoef.test[attr(name.coef,"originalOnly"),,drop=FALSE]
            Mrescue.test <- do.call(cbind,lapply(equationComplex.lhs, FUN = function(iEq){sapply(attr(name.coef,"rescueOnly"), FUN = grepl, x = iEq)}))
            if(any(Mrescue.test) && length(index.keep)>0 && any(equation.lhs[index.coef] %in% attr(name.coef,"rescue") == FALSE)){
                stop("Incorrect argument \'",name.arg,"\': contains untransformed and transformed parameter names. \n",
                     "Problematic example: \"",equation[index.coef][1],"\" and \"",equation[colSums(Mrescue.test)>0][1],"\". \n")
            }else if(all(colSums(Mrescue.test) >= colSums(Moriginal.test))){
                rescue <- TRUE
                name.coef <- attr(name.coef,"rescue")
                equation.contrast <- matrix(0, nrow = n.equation, ncol = n.coef, dimnames = list(equation, name.coef))
                object.contrast[equation] <- lapply(equation, function(iE){equation.contrast[iE,,drop=FALSE]})       
                names(stdname.coef) <- name.coef
                Mcoef.test <- do.call(rbind,lapply(name.coef, function(iCoef){sapply(equationComplex.lhs, grepl, pattern = iCoef, fixed = TRUE)}))
                rownames(Mcoef.test) <- name.coef
            }else if(all(colSums(Mrescue.test) <= colSums(Moriginal.test))){
                attr(name.coef,"rescue") <- NULL
                attr(name.coef,"rescueOnly") <- NULL
                attr(name.coef,"originalOnly") <- NULL
            }else{
                stop("Incorrect argument \'",name.arg,"\': contains untransformed and transformed parameter names. \n",
                     "Problematic example: \"",equation[colSums(Mrescue.test)<colSums(Moriginal.test)][1],"\" and \"",equation[colSums(Mrescue.test)>colSums(Moriginal.test)][1],"\". \n")
            }
        }

        if(any(colSums(Mcoef.test)==0)){
            stop("Equation \"",paste(equationComplex[colSums(Mcoef.test)==0], collapse = "\", \""),"\" does not contain recognizable model parameters. \n")
        }

        ## standardize equation
        equationStd.lhs <- sapply(equationComplex.lhs, FUN = function(iEq){ ## iEq <- equation.lhs[2]
            iCoef.all <- names(which(Mcoef.test[,iEq]))
            iCoef.all.order <- iCoef.all[order(nchar(iCoef.all), decreasing = TRUE)] ## substitute longer string first to avoid confusion, e.g. between sigma.12 and sigma.1
            for(iCoef in iCoef.all.order){ ## iCoef <- iCoef.all[2]
                iEq <- gsub(pattern = iCoef, replacement = stdname.coef[iCoef], x = iEq, fixed = TRUE)
            }
            return(iEq)
        })

        ## **** option 1: via multcomp
        ## object.contrast[equationComplex] <- lapply(1:n.equationComplex, function(iE){ ## iE <- 1
        ##     iLS <- multcomp:::expression2coef(parse(text = paste0(equationStd.lhs[iE],"=",0)), stdname.coef)
        ##     iContrast <- object.contrast[[equationComplex[iE]]]
        ##     iContrast[1,match(iLS$names,stdname.coef)] <- iLS$coef
        ##     return(iContrast)
        ## })
        
        ## **** option 2: by hand
        equationStd.lhs <- gsub(pattern = " ",replacement = "", x = equationStd.lhs, fixed = TRUE)

        ## futher standardize equation ... + ... + ...
        equationStd2.lhs <- sapply(1:n.equationComplex, function(iE){ ## iE <- 1
            iEq <- equationStd.lhs[iE]
            iCoef.all <- names(which(Mcoef.test[,equationComplex.lhs[iE]]))
            ## add missing 1*
            for(iCoef in iCoef.all){ ## iCoef <- "paramBD"
                ## if leading without multiplier add *1, e.g. age+... -> 1*age+... or age-... -> 1*age+-...
                iEq <- gsub(pattern = paste0("^",stdname.coef[iCoef],"\\+"), replacement = paste0("1\\*",stdname.coef[iCoef],"\\+"), x = iEq, fixed = FALSE)
                iEq <- gsub(pattern = paste0("^",stdname.coef[iCoef],"\\-"), replacement = paste0("1\\*",stdname.coef[iCoef],"\\-"), x = iEq, fixed = FALSE)
                ## if in the middle without multiplier or division then add 1*, e.g. +age+... -> +1*age+... or +age-... -> +1*age+-...
                iEq <- gsub(pattern = paste0("+",stdname.coef[iCoef],"+"), replacement = paste0("+1*",stdname.coef[iCoef],"+"), x = iEq, fixed = TRUE)
                iEq <- gsub(pattern = paste0("+",stdname.coef[iCoef],"-"), replacement = paste0("+-1*",stdname.coef[iCoef],"-"), x = iEq, fixed = TRUE)
                iEq <- gsub(pattern = paste0("-",stdname.coef[iCoef],"+"), replacement = paste0("+1*",stdname.coef[iCoef],"+"), x = iEq, fixed = TRUE)
                iEq <- gsub(pattern = paste0("-",stdname.coef[iCoef],"-"), replacement = paste0("+-1*",stdname.coef[iCoef],"-"), x = iEq, fixed = TRUE)
                ## if in the middle with multiplier and - add +
                iEq <- gsub(pattern = paste0("*",stdname.coef[iCoef],"-"), replacement = paste0("*",stdname.coef[iCoef],"+-"), x = iEq, fixed = TRUE)
                iEq <- gsub(pattern = paste0("/",stdname.coef[iCoef],"-"), replacement = paste0("/",stdname.coef[iCoef],"+-"), x = iEq, fixed = TRUE)
                ## if ending without multiplier or division then add 1*
                iEq <- gsub(pattern = paste0("\\+",stdname.coef[iCoef],"$"), replacement = paste0("\\+1\\*",stdname.coef[iCoef]), x = iEq, fixed = FALSE)
                iEq <- gsub(pattern = paste0("\\-",stdname.coef[iCoef],"$"), replacement = paste0("\\-1\\*",stdname.coef[iCoef]), x = iEq, fixed = FALSE)
            }
            return(iEq)            
        })

        ## decode
        for(iE in 1:n.equationComplex){
            iEq.split <- strsplit(x = equationStd2.lhs[iE], split = "+", fixed = TRUE)[[1]]
            ## identify sign (TRUE -, FALSE +)
            iMinus.split <- sapply(iEq.split,grepl, pattern = "-")
            ## identify coefficient
            iCoef.all <- names(which(Mcoef.test[,equationComplex.lhs[iE]]))
            iLs.termcoef <- lapply(iEq.split, function(iEq){
                iCoef.all[sapply(stdname.coef[iCoef.all], grepl, x = iEq, fixed = TRUE)]
            })
            ## check single coefficient per term
            if(any(lengths(iLs.termcoef)==0)){
                stop("Unknown element \"",paste(iEq.split[which(lengths(iLs.termcoef)==0)], collapse = "\", \""),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                     "Constant terms are not recognized on the left hand side of the equation. \n")
            }else if(any(lengths(iLs.termcoef)>1)){
                pb <- iEq.split[which(lengths(iLs.termcoef)>1)]
                for(iCoef in iCoef.all){
                    pb <- gsub(pattern = stdname.coef[iCoef], replacement = iCoef, pb)
                }
                ## check whether the problem comes from the multiplier not being on the same side of the coefficient,
                ## e.g. "theta * 5 - 6 * rho" will not get standardized into "theta * 5 +- 6 * rho"
                iEq.split2 <- strsplit(x = iEq.split, split = "-")[[1]]
                iLs.termcoef2 <- lapply(iEq.split2, function(iEq){iCoef.all[sapply(stdname.coef[iCoef.all], grepl, x = iEq, fixed = TRUE)]})
                if(any(lengths(iLs.termcoef2)>1)){
                    stop("Unknown element \"",paste(pb, collapse = "\", \""),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                         "Terms containing multiple parameters, e.g. products, are not recognized on the left hand side of the equation. \n")
                }else if(any(!grepl(pattern = "\\*|\\/", setdiff(iEq.split2,stdname.coef)))){
                    stop("Unknown element \"",paste(pb, collapse = "\", \""),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                         "Missing operator (* or /) between the constant and the model parameter. \n")
                }else{
                    stop("Unknown element \"",paste(pb, collapse = "\", \""),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                         "Constant multiplying the model parameters should always be on the same side of the model parameter (always before or always after). \n")
                }
            }
            iCoef.split <- unlist(iLs.termcoef)
            if(any(duplicated(iCoef.split))){                
                stop("Multiple occurence of the model parameter \"",paste(unique(iCoef.split[duplicated(iCoef.split)]), collapse ="\", \""),"\" in the equation \"",equation[iE],"\". \n")
            }
            ## identify constant
            iFactor.split <- mapply(iCoef = iCoef.split, iSplit = iEq.split, iSign = iMinus.split, function(iCoef,iSplit,iSign){
                ## remove parameter name
                iSplitRed <- gsub(pattern = stdname.coef[iCoef], replacement = "", x = iSplit, fixed = TRUE)
                if(iSign){ ## remove sign
                    iSplitRed <- gsub(pattern = "-", replacement = "", x = iSplitRed, fixed = TRUE)
                }
                ## identify operator
                iOperator <- NULL
                if(grepl("^\\*",iSplitRed,fixed = FALSE)){
                    iOperator <- 1
                    iSplitRed <- gsub(pattern = "^\\*", replacement = "", iSplitRed)
                }else if(grepl("^\\/",iSplitRed,fixed = FALSE)){
                    iOperator <- 2
                    iSplitRed <- gsub(pattern = "^\\/", replacement = "", iSplitRed)
                }
                if(grepl("\\*$",iSplitRed,fixed = FALSE)){
                    iOperator <- c(iOperator,1)
                    iSplitRed <- gsub(pattern = "\\*$", replacement = "", iSplitRed)
                }else if(grepl("\\/$",iSplitRed,fixed = FALSE)){
                    iOperator <- c(iOperator,2)
                    iSplitRed <- gsub(pattern = "\\/$", replacement = "", iSplitRed)
                }
                if(length(iOperator)==0){
                    if(length(strsplit(x = iSplit, split = stdname.coef[iCoef], fixed = TRUE)[[1]])>1){
                        stop("Unknown element \"",gsub(pattern = stdname.coef[iCoef], replacement = iCoef, x = iSplit),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                             "Cannot handle constant on both sides of the model parameter. \n")
                    }else{
                        stop("Unknown element \"",gsub(pattern = stdname.coef[iCoef], replacement = iCoef, x = iSplit),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                             "Missing operator (* or /) between the constant and the model parameter. \n")
                    }                    
                }else if(length(iOperator)>1){
                    stop("Unknown element \"",gsub(pattern = stdname.coef[iCoef], replacement = iCoef, x = iSplit),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                         "Too many operator (* or /) next to the model parameter. \n")
                }
                iOut <- try(eval(parse(text = iSplitRed)), silent = TRUE)
                if(inherits(iOut, "try-error") | !is.numeric(iOut)){
                    stop("Unknown element \"",gsub(pattern = stdname.coef[iCoef], replacement = iCoef, x = iSplit),"\" in the left hand side of the equation in argument \'",name.arg,"\'. \n",
                         "Cannot convert constant to numeric \n")
                }
                if(iOperator==2){
                    iOut <- 1/iOut
                }
                if(iSign){
                    iOut <- -iOut
                }
                return(iOut)
            })
            ## update contrast matrix
            object.contrast[[equationComplex[[iE]]]][1,iCoef.split] <- as.double(iFactor.split)
        }
    }

    ## ** export
    out <- list(contrast = do.call(rbind,object.contrast),
                rhs = do.call(c,object.rhs),
                rescue = rescue)
    return(out)
}


## * contrast2name
##' @description Find a name corresponding to each contrast.
##' Used by rbind.Wald_lmm.
##'
##' @param contrast [matrix] constrast matrix.
##' @param MARGIN [1 or 2] typically 1 when contrasts are encoded per row over model parameters in columns.
##' @param simplify [logical] should the result be output as matrix instead of a list?
##' @param ignore.value [logical] should the value of the contrast be ignored?
##' @param collapse [logical] 
##' 
##' @noRd
##'
##' @examples
##' M <- matrix(0,3,3, dimnames = list(NULL, c("A","B","C")))
##' M[1,1] <- M[2,1] <- M[2,2] <- M[3,3] <- 1
##' 
##' contrast2name(M)
##' contrast2name(M, collapse = TRUE)
contrast2name <- function(contrast, MARGIN = 1, simplify = collapse, ignore.value = FALSE, collapse = FALSE, sep = ", "){

    ## ** combine
    ls.out <- apply(contrast, MARGIN = MARGIN, FUN = function(iVec){
        iTest.n0 <- which(iVec!=0)
        iOut <- data.frame(name = names(iTest.n0),
                           value = unname(iVec[iTest.n0]))
        if(ignore.value){
            iOut$value <- NULL
            if(collapse){
                iOut <- paste(iOut$name, collapse = sep)
            }
        }else if(collapse){
            iValue <- iOut$value
            iValue[iOut$value==1] <- ""
            iValue[iOut$value==-1] <- "-"
            iOut <- paste(paste0(iValue,iOut$name), collapse = sep)
        }
        return(iOut)
    }, simplify = FALSE)

    ## ** export
    if(simplify){
        if(collapse){
            out <- unlist(ls.out)
        }else{
            out <- do.call(rbind,lapply(1:length(ls.out), FUN = function(iOut){
                cbind(index = iOut, ls.out[[iOut]])
            }))
        }
    }else{
        out <- ls.out
    }
    return(out)
}

##----------------------------------------------------------------------
### contrast.R ends here
