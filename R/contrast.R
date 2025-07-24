### contrast.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2024 (09:37) 
## Version: 
## Last-Updated: jul 24 2025 (16:02) 
##           By: Brice Ozenne
##     Update #: 383
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
    object.eigen <- eigen(tcrossprod(object), symmetric = TRUE)
    n.zero <- sum(abs(object.eigen$values) < tol)
    if(length(rhs)==1 && NROW(object)>1){
        rhs <- rep(rhs, NROW(object))
    }

    if(n.zero==0){return(list(C = object, rhs = rhs, dim = NROW(object), rm = 0))}

    keep.lines <- 1:NROW(object)
    for(iLine in NROW(object):1){ ## iLine <- 1
        iN.zero <- sum(abs(eigen(tcrossprod(object[setdiff(keep.lines,iLine),,drop=FALSE]), symmetric = TRUE)$values) < tol)
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
##' @description Find the contrast matrix corresponding to an equation
##' Used by anova.lmm.
##' @noRd
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
        if(is.null(names(equation))){
            names(equation) <- equation
        }
        equation.contrast <- matrix(0, nrow = n.equation, ncol = n.coef, dimnames = list(names(equation), name.coef))
        object.contrast[equation] <- lapply(names(equation), function(iE){equation.contrast[iE,,drop=FALSE]})

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
                equation.contrast <- matrix(0, nrow = n.equation, ncol = n.coef, dimnames = list(names(equation), name.coef))
                for(iE in equation){
                    colnames(object.contrast[equation]) <- name.coef
                }
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
            }else if(any(Mrescue.test) && all(colSums(Mrescue.test) >= colSums(Moriginal.test))){
                rescue <- TRUE
                name.coef <- attr(name.coef,"rescue")
                for(iE in equation){
                    colnames(object.contrast[equation]) <- name.coef
                }
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
        equationStd.lhs <- sapply(names(equationComplex), FUN = function(iEq){ ## iEq <- equationComplex[1]
            iCoef.all <- names(which(Mcoef.test[,iEq]))
            iCoef.all.order <- iCoef.all[order(nchar(iCoef.all), decreasing = TRUE)] ## substitute longer string first to avoid confusion, e.g. between sigma.12 and sigma.1
            iEq.lhs <- equationComplex.lhs[iEq]
            for(iCoef in iCoef.all.order){ ## iCoef <- iCoef.all[2]
                iEq.lhs <- gsub(pattern = iCoef, replacement = stdname.coef[iCoef], x = iEq.lhs, fixed = TRUE)
            }
            return(iEq.lhs)
        }, USE.NAMES = FALSE)

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
            iCoef.all <- names(which(Mcoef.test[,names(equationComplex)[iE]]))
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
            iCoef.all <- names(which(Mcoef.test[,names(equationComplex)[iE]]))
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
            if(all(iFactor.split==0)){
                stop("Incorrect value for argument \'",name.arg,"\': cannot have an equation where all coefficients equal 0. \n")
            }
            object.contrast[[equationComplex[[iE]]]][1,iCoef.split] <- as.double(iFactor.split)
        }
    }

    ## ** export
    out <- list(contrast = do.call(rbind,object.contrast),
                rhs = do.call(c,object.rhs),
                rescue = rescue)
    return(out)
}

## * effects2contrast 
##' @description Generate a contrast matrix and a right-hand side based on the effects argument
##' Used by anova.lmm.
##' @noRd
effects2contrast <- function(object, effects, rhs, options){

    effects.ref <- c("all","mean","fixed","variance","correlation")
    transform.sigma <- object$reparametrize$transform.sigma
    transform.k <- object$reparametrize$transform.k
    transform.rho <- object$reparametrize$transform.rho

    ## ** extract from object
    object.coef <- stats::model.tables(object, effects = "param",
                                       transform.sigma = transform.sigma,
                                       transform.k = transform.k,
                                       transform.rho = transform.rho,
                                       transform.names = TRUE)
    name.coef <- object.coef$name
    n.coef <- length(name.coef)
    type.coef <- stats::setNames(object.coef$type, name.coef)

    name.coef.rescue <- name.coef
    attr(name.coef.rescue, "rescue") <- object.coef$trans.name

    if(is.character(effects)){
        object.X <- stats::model.matrix(object, effects = "all", simplify = 0.5, options = options)
        
        strata.var <- object$strata$var
        strata.levels <- object$strata$levels        
        name.coef.sigma <- object.coef$name[object.coef$type == "sigma"]
        if(all(is.na(attr(strata.var,"original")))){
            name.strata.sigma <- NULL
        }else{
            name.strata.sigma <- paste0(strata.var, strata.levels[object.coef[match(name.coef.sigma,object.coef$name),"index.strata"]])
        }
    }

    ## ** normalize input
    if(is.character(effects)){

        if(all(tolower(effects) %in% effects.ref)){
            effects <- tolower(effects)

            if("all" %in% effects){
                if(length(effects) == 1){
                    effects <- c("mean","variance","correlation")
                }else{
                    stop("When argument \'effect\' contains \"all\" it should be of length 1. \n")
                }
            }else if("fixed" %in% effects){
                effects[effects=="fixed"] <- "mean"
            }
            if(transform.k %in% c("sd","var","logsd","logvar")){
                stop("Cannot use \'transform.k\' equals \"sd\", \"var\", \"logsd\", or \"logvar\". \n",
                     "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
            }
            if(all(attr(object.X$mean,"assign")==0)){
                effects <- setdiff(effects,"mean")
            }
            if("k" %in% type.coef == FALSE){
                effects <- setdiff(effects,"variance")
            }
            if("rho" %in% type.coef == FALSE){
                effects <- setdiff(effects,"correlation")
            }
            if(length(effects)==0){
                message("No variable relative to which model parameters should be tested. \n")
                return(invisible(NULL))
            }
        }

    }else if(is.matrix(effects)){

        ## effects
        if(any(is.na(effects))){
            stop("When a matrix, argument \'effects\' should no contain NA values. \n")
        }
        if(any(rowSums(effects!=0)==0)){
            stop("When a matrix, argument \'effects\' should no contain rows with only 0. \n")
        }
        if(is.null(rownames(effects))){
            stop("When a matrix, argument \'effects\' should have row names. \n")
        }
        if(any(duplicated(rownames(effects)))){
            stop("When a matrix, argument \'effects\' should not have duplicated row names. \n")
        }
        if(is.null(colnames(effects))){
            if(NCOL(effects)!=n.coef){
                stop("When a matrix, argument \'effects\' should have column names or as many columns as model parameters (here ",n.coef,"). \n")
            }
            colnames(effects) <- name.coef
        }else{
            if(any(duplicated(colnames(effects)))){
                stop("When a matrix, argument \'effects\' should not have duplicated column names. \n")
            }            
            if(any(colnames(effects) %in% name.coef == FALSE)){
                if(any(colnames(effects) %in% attr(name.coef.rescue,"rescue") == FALSE)){
                    stop("When a matrix, argument \'effects\' should have column names matching the names of the model parameters. \n",
                         "Valid names: \"",paste(name.coef, collapse = "\", \""),"\". \n")
                }else{
                    colnames(effects) <- name.coef.rescue[match(colnames(effects), attr(name.coef.rescue,"rescue"))]
                }
            }
            if(NCOL(effects)!=n.coef){
                effectsSave <- effects
                effects <- matrix(0, nrow = NROW(effectsSave), ncol = n.coef,
                                  dimnames = list(rownames(effectsSave),name.coef))
                effects[,colnames(effectsSave)] <- effectsSave
            }
            
        }

        ## rhs
        if(is.null(rhs)){
            effects.type <- do.call(rbind,apply(effects, MARGIN = 1, function(iRow){tapply(iRow!=0,type.coef,sum)}, simplify=FALSE))
            if("k" %in% names(effects.type) == FALSE){effects.type <- cbind(effects.type, k = 0)} ## case of a homoschedastic model

            if(any(effects.type[,"sigma"]!=0) || (any(effects.type[,"k"]!=0) && transform.k %in% c("sd","var","logsd","logvar"))){
                stop("Unable to decide on a value for the right-hand side of the null hypothesis\n",
                     "due to the presence of parameters of type \"sigma\" in the same null hypothesis. \n",
                     "Consider specifying the argument \'rhs\'. \n")
            }
            if(any(rowSums(effects.type>0)>1)){
                stop("Unable to decide on a value for the right-hand side of the null hypothesis\n",
                     "due to multiple types of parameters in the same null hypothesis. \n",
                     "Consider specifying the argument \'rhs\'. \n")
            }
        }else if(!is.numeric(rhs) || !is.vector(rhs)){
            stop("Argument \'rhs\' should be a numeric vector. \n")
        }else if(!is.null(names(rhs)) && any(rownames(effects) %in% names(rhs) == FALSE)){
            stop("Argument \'rhs\' should have the same names as the row names of argument \'effects\'. \n")
        }else if(length(rhs)!=1 && length(rhs)!=NROW(effects)){
            stop("Incorrect length for argument \'rhs\': should have length the number of rows of the contrast matrix (here ",NROW(effects),"). \n")
        }

    }else if(!inherits(effects,"mcp")){
        stop("Incorrect argument 'effects': can be \"mean\", \"variance\", \"correlation\", \"all\", \n", 
             "or an equation such compatible with the argument 'linfct' of multcomp::glht \n ", 
             "or a contrast matrix. \n", "or covariate names \n ")
    }

    ## ** Generate contrast matrix
    if(all(effects %in% effects.ref)){
        ## *** Case 1: effects refer to a type of parameter

        ls.contrast <- list()
        ls.null <- list()

        if("mean" %in% effects){
            
            ## names of the terms in the design matrix
            terms.mean <- attr(object.X$mean,"term.labels")

            ## contrast matrix
            ls.contrast$mu <- stats::setNames(lapply(1:length(terms.mean), function(iT){ ## iT <- 1
                iIndex <- which(attr(object.X$mean,"assign")==iT)
                iCoef <- colnames(object.X$mean)[iIndex]
                iC <- matrix(0, nrow = length(iIndex), ncol = n.coef,
                             dimnames = list(paste0(iCoef,"=0"), name.coef))
                iC[,iCoef] <- diag(1, nrow = length(iIndex))
                return(iC)                
            }), terms.mean)
            ls.null$mu <- lapply(ls.contrast$mu, function(iC){stats::setNames(rep(0, NROW(iC)),rownames(iC))})
        }

        if("variance" %in% effects){
            
            ## names of the k coefficients
            name.coef.k <- name.coef[type.coef == "k"]
            n.coef.k <- length(name.coef.k)

            ## terms
            if(is.null(name.strata.sigma)){
                ## no strata: test each term
                terms.var <- unique(attr(object.X$var,"term.labels") [match(name.coef.k, colnames(object.X$var))])
            }else if(all(attr(object.X$var,"variable") %in% strata.var | attr(object.X$var,"variable") %in% attr(strata.var,"original"))){
                ## no term: test strata
                terms.var <- rep(paste(attr(strata.var,"original"),collapse=", "), length(name.strata.sigma))
            }else{
                ## term and strata: test term(s) within strata
                terms.var <- sapply(name.strata.sigma, function(iStrata){paste(iStrata,setdiff(attr(object.X$var,"variable"), strata.var), sep = ":")})
            }
            Uterms.var <- unique(terms.var)

            ## contrast matrix 
            ls.contrast$k <- stats::setNames(lapply(Uterms.var, function(iTerm){ ## iTerm <- Uterms.var[1]
                iCoef <- name.coef.k[terms.var==iTerm]
                iC <- matrix(0, nrow = length(iCoef), ncol = n.coef,
                             dimnames = list(paste0(iCoef,"=1"), name.coef))
                iC[,iCoef] <- diag(1, nrow = length(iCoef))
                return(iC)                
            }),terms.var)
            ls.null$k <- lapply(ls.contrast$k, function(iC){stats::setNames(rep(1, NROW(iC)),rownames(iC))})
        }

        if("correlation" %in% effects){

            ## names of the rho coefficients
            name.coef.rho <- name.coef[type.coef == "rho"]
            n.coef.rho <- length(name.coef.rho)
            if(is.null(name.strata.sigma)){
                ## no strata: test each term
                terms.cor <- unique(attr(object.X$cor,"term.labels"))
            }else if(all(attr(object.X$cor,"variable") %in% strata.var | attr(object.X$cor,"variable") %in% attr(strata.var,"original"))){
                ## no term: test strata
                terms.cor <- rep(paste(attr(strata.var,"original"),collapse=", "), length(name.strata.sigma))
            }else{
                ## term and strata: test term(s) within strata
                terms.cor <- sapply(name.strata.sigma, function(iStrata){paste(iStrata,setdiff(attr(object.X$cor,"variable"), strata.var), sep = ":")})
            }
            Uterms.cor <- unique(terms.cor)

            ## contrast matrix 
            ls.contrast$rho <- stats::setNames(lapply(Uterms.cor, function(iTerm){ ## iTerm <- Uterms.cor[1]
                iCoef <- name.coef.rho[terms.cor==iTerm]
                iC <- matrix(0, nrow = length(iCoef), ncol = n.coef,
                             dimnames = list(paste0(iCoef,"=0"), name.coef))
                iC[,iCoef] <- diag(1, nrow = length(iCoef))
                return(iC)                
            }),Uterms.cor)
            ls.null$rho <- lapply(ls.contrast$rho, function(iC){stats::setNames(rep(0, NROW(iC)),rownames(iC))})
        }
        
        backtransform <- TRUE

    }else if(inherits(effects,"mcp")){
        ## *** Case 2: effects refer to a mcp object (multcomp package)

        out.m2c <- try(multcomp::glht(object, linfct = effects), silent = TRUE)
        
        if(inherits(out.m2c,"try-error")){
            valid.type <- eval(formals(multcomp::contrMat)$type)
            if(length(effects)==1 && length(effects[[1]])==1 && is.character(effects[[1]]) && effects[[1]] %in% valid.type == FALSE){
                stop("Incorrect argument \'effects\': running the mcp2matrix function from the mulcomp package lead to the following error: \n",
                     out.m2c,
                     "Valid contrast: \"",paste(valid.type, collapse ="\", \""),"\". \n")
            }else{
                stop("Incorrect argument \'effects\': running the mcp2matrix function from the mulcomp package lead to the following error: \n",
                     out.m2c)
            }
        }
        
        ls.contrast <- list(user = list(user = matrix(0, nrow = NROW(out.m2c$linfct), ncol = n.coef,
                                                      dimnames = list(rownames(out.m2c$linfct), name.coef))))
        ls.contrast$user$user[,colnames(out.m2c$linfct)] <- out.m2c$linfct
        ls.null  <- list(user = list(user = stats::setNames(out.m2c$rhs, rownames(out.m2c$linfct))))
        backtransform <- TRUE

    }else if(is.matrix(effects)){
        ## *** Case 3: effects refer to contrast matrix
        
        if(is.null(colnames(effects))){
            colnames(effects) <- name.coef
            rescue <- FALSE
        }else{
            ## it has been checked before that when some names dot not belong to the untransformed names, all names belong to the transformed names
            rescue <- any(colnames(effects) %in% name.coef == FALSE)

            effects.save <- effects
            effects <- matrix(0, nrow = NROW(effects.save), ncol = length(name.coef), dimnames = list(rownames(effects.save),name.coef))
            effects[,colnames(effects.save)] <- effects.save
        }

        ## normalize rhs
        if(is.null(rhs)){
            effects.type2 <- colnames(effects.type)[apply(effects.type!=0,1,which)]
            default.null <- c("mu" = 0,
                              "k" = ifelse(transform.k %in% c("log","logsquare"), 0, 1),
                              "rho" = 0)
            rhs <- stats::setNames(default.null[effects.type2], rownames(effects))
        }else if(length(rhs)==1){
            rhs <- stats::setNames(rep(rhs, NROW(effects)), rownames(effects))
        }else if(is.null(names(rhs))){
            names(rhs) <- rownames(effects)
        }else{
            rhs <- rhs[rownames(effects)]
        }

        ## collect
        if(rescue==TRUE && (!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho))){ ## the user specifically requests the transformed scale
            backtransform <- FALSE
        }else{
            backtransform <- TRUE
        }
        ls.contrast <- list(user = list(user = effects))
        ls.null  <- list(user = list(user = rhs))
        
    }else if(is.character(effects)){
        ## *** Case 4: effects refer to variable names or equations
        out.eq2c <- equation2contrast(effects,
                                      name.coef = name.coef.rescue,
                                      X = list(mean = object.X$mean, var = object.X$var),
                                      name.arg = "effects")

        if(out.eq2c$rescue==TRUE){
            colnames(out.eq2c$contrast) <- name.coef
            if(!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho)){ ## the user specifically requests the transformed scale
                backtransform <- FALSE
            }else{
                backtransform <- TRUE
            }
        }else{
            backtransform <- TRUE
        }
        ls.contrast <- list(user = list(user = out.eq2c$contrast))
        ls.null  <- list(user = list(user = out.eq2c$rhs))

    }

    ## ** export
    return(list(contrast = ls.contrast, null = ls.null, backtransform = backtransform))
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
