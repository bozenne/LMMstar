### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: okt  7 2020 (11:58) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmm (documentation)
##' @title Linear mixed model
##' @description Fit a linear mixed model using either a compound symmetry structure or an unstructured covariance matrix.
##' This is essentially an interface to the \code{nlme::gls} function.
##'
##' @param formula [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param covariance [formula] Specify the model for the covariance.
##' No left hand side. On the right hand side,
##' either only the grouping variable (when specifying a compound symmetry structure), e.g. ~1|id,
##' or the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param ... passed to \code{nlme::gls}.

## * lmm (examples)
##' @examples
##' ## simulate data in the wide format
##' library(lava)
##' m <- lvm(c(Y1,Y2,Y3,Y4) ~ age + gender)
##' categorical(m, labels = c("male","female")) <- ~gender
##' transform(m, id~gender) <- function(x){1:NROW(x)}
##' distribution(m, ~age) <- gaussian.lvm(mean = 50, sd = 10)
##'
##' set.seed(10)
##' dW <- lava::sim(m, 1e2)
##'
##' ## move to the long format
##' name.varying <- paste0("Y",1:4)
##' dL <- reshape(dW, direction  = "long",
##'               idvar = c("id","age","gender"),
##'               varying = name.varying,
##'               v.names = "Y",
##'               timevar = "visit")
##' rownames(dL) <- NULL
##' dL$visit <- factor(dL$visit,
##'                    levels = 1:length(name.varying),
##'                    labels = name.varying)
##' head(dL)
##' 
##' ## fit mixed model
##' e0.lmm <- lmm(Y ~ visit + age + gender, covariance = ~1|id, data = dL)
##' summary(e0.lmm)
##' nlme:::summary.gls(e0.lmm)
##'
##' e.lmm <- lmm(Y ~ visit + age + gender, covariance = ~visit|id, data = dL)
##' cat(attr(e.lmm,"code")) ## code used to fit the model
##' head(attr(e.lmm,"data")) ## data used to fit the model
##' summary(e.lmm)

## * lmm (code)
##' @export
lmm <- function(formula, covariance, data, df = FALSE, ...){

    ## ** check and normalize user imput
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' must be of class formula \n",
             "Something like: outcome ~ fixedEffect1 + fixedEffect2 \n")
    }
    name.mean <- all.vars(formula)
    if(any(name.mean %in% names(data) == FALSE)){
        invalid <- name.mean[name.mean %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    if(!inherits(covariance,"formula")){
        stop("Argument \'covariance\' must be of class formula. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }
    if(length(lhs.vars(covariance))!=0){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "Should not have any variable on the left hand side. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }

    if(!grepl("|",deparse(covariance),fixed = TRUE)){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "No | symbol found so no grouping variable could be defined. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }

    if(length(grepl("|",deparse(covariance),fixed = TRUE))>1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "The symbol | should only appear once. \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }
    res.split <- strsplit(deparse(covariance),"|", fixed = TRUE)[[1]]
    name.cluster <- trimws(res.split[2], which = "both")
    name.repetition <- all.vars(stats::as.formula(res.split[1]))
    if(length(name.repetition)>1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "Should have at most one variable before the grouping symbol (|). \n",
             "Shoud be something like: ~ 1|id (compound symmetry) or ~ time|id (unstructured). \n")
    }

    ## ** fit mixed model
    if(length(name.repetition)==0){
        form.cor <- stats::as.formula(paste0("~1|",name.cluster))
        
        e.lmm <- eval(parse(text = paste0("nlme::gls(",deparse(formula),",
                     correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                     data = data,
                     ...)")))
    }else{
        name.repetition.index <- paste0(name.repetition,".index")
        if(name.repetition.index %in% names(data)){
            stop("Incorrect specification of argument \'data\'. \n",
                 "The variable ",name.repetition.index," is used internally but already exists in \'data\' \n")
        }
        data[[name.repetition.index]] <- as.numeric(as.factor(data[[name.repetition]]))

        form.cor <- stats::as.formula(paste0("~",name.repetition.index,"|",name.cluster))
        form.var <- stats::as.formula(paste0("~1|",name.repetition))
        e.lmm <- eval(parse(text = paste0("nlme::gls(",deparse(formula),",
                     correlation = nlme::corSymm(form = ",deparse(form.cor),"),
                     weights = nlme::varIdent(form = ",deparse(form.var),"),
                     data = data,
                     ...)")))
    }

    ## ** small sample correction
    if(df){
        stop("not implemented yet!")
        ## sCorrect(e.lmm) <- TRUE
    }

    ## ** export
    attr(e.lmm, "code") <- paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(e.lmm$call), collapse = ""))),"\n")
    attr(e.lmm, "data") <- data
    class(e.lmm) <- append("lmm",class(e.lmm))
    return(e.lmm)
}


######################################################################
### lmm.R ends here
