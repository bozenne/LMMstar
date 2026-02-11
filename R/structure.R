### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: feb 11 2026 (13:16) 
##           By: Brice Ozenne
##     Update #: 1486
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ID (identity, documentation)
##' @title Identity Structure
##' @description Variance-covariance structure where the residuals are independent and identically distributed.
##' Can be stratified on a categorical variable.
##' 
##' @param formula [formula] variables to stratify the residual variance (left hand side, like sex~1).
##' @param var.time not used.
##'
##' @details A typical formula would be \code{~1}. It will model: \itemize{
##' \item \eqn{\sigma^2_{j}=\sigma^2}: a variance constant over repetitions.
##' \item \eqn{\rho_{j,j'}=0}: no correlation.
##' }
##' With 4 repetitions (\eqn{(j,j') \in \{1,2,3,4\}^2} and \eqn{j\neq j'}) this leads to the following residuals variance-covariance matrix:
##' \deqn{ \begin{bmatrix} \sigma^2 & 0 & 0 & 0 \\ 0 & \sigma^2 & 0 & 0 \\ 0 & 0 & \sigma^2 & 0 \\ 0 & 0 & 0 & \sigma^2 \end{bmatrix}}
##' It is possibly to stratify covariance pattern on a categorical variable (e.g. sex) using a formula such as \code{~sex}. This will lead to \eqn{\sigma^2_{j,sex}=\sigma_{sex}^2} and \eqn{\rho_{j,j',sex}=0}, i.e.:
##' \deqn{ \begin{bmatrix} \sigma_{\text{female}}^2 & 0 & 0 & 0 \\ 0 & \sigma_{\text{female}}^2 & 0 & 0 \\ 0 & 0 & \sigma_{\text{female}}^2 & 0 \\ 0 & 0 & 0 & \sigma_{\text{female}}^2 \end{bmatrix}
##' \text{ and }
##'       \begin{bmatrix} \sigma_{\text{male}}^2 & 0 & 0 & 0 \\ 0 & \sigma_{\text{male}}^2 & 0 & 0 \\ 0 & 0 & \sigma_{\text{male}}^2 & 0 \\ 0 & 0 & 0 & \sigma_{\text{male}}^2 \end{bmatrix}
##' }
##'
##' @return An object of class \code{ID} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @keywords multivariate
##' 
##' @examples
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$sex <- factor(as.numeric(gastricbypassL$id) %% 2,
##'                              labels = c("female", "male"))
##'
##' ## default: single variance parameter
##' eID.default <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = ID, data = gastricbypassL)
##' sigma(eID.default)
##' model.tables(eID.default, effects = "variance")
##' 
##' ## stratified: one variance parameter per strata
##' eID.strata <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = ID(sex~1), data = gastricbypassL)
##' sigma(eID.strata)
##' model.tables(eID.strata, effects = "variance")
##' 
##' ## same as
##' eID.reg <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = ID(~sex), data = gastricbypassL)
##' sigma(eID.reg)
##' model.tables(eID.reg, effects = "variance")
##' }

## * ID (identity, code)
##' @export
ID <- function(formula, var.time){

    ## ** check input
    ## *** formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' should be a formula. \n")
    }

    ## *** var.time
    if(!missing(var.time)){
        if(any(inherits(var.time,"formula"))){ ## possible user mistake ID(~1,~1) instead of ID(list(~1,~1))
            stop("Argument \'var.time\' should not be a formula. \n",
                 "Consider using a list to collect the formula for the variance and correlation structure. \n")
        }else{
            warning("Argument \'var.time\' ignored with ID structure. \n")
        }
    }

    

    ## ** identify strata variable and covariate variables
    outCov <- .formulaStructure(formula, add.X = NULL, strata.X = FALSE, correlation = FALSE)
    if(length(outCov$X.var)>0 || length(outCov$X.cor)>0){ 
        stop("Argument \'formula\' should not contain variable(s) on the right-hand side of the formula. \n",
             "Structure ID can only be stratified (left hand side of the formula). \n")
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = NA,
                                  variance = NA,
                                  correlation = NA,
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula.var,
                               correlation = NULL),
                class = c(variance = "ID", correlation = NA))

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("ID",class(out))
    return(out)
}

## * IND (independence, documentation)
##' @title Independence Structure
##' @description Variance-covariance structure where the residuals are independent but may have different variance.
##' Can be stratified on a categorical variable.
##'
##' @param formula [formula] variables influencing the residual variance,
##' using either as a multiplicative factor (right hand side) or stratification (left hand side) to model their effect.
##' @param var.time Missing or \code{NULL} to enable or disable repetition-specific variance parameters
##'
##' @details A typical formula would be \code{~time}. It will model: \itemize{
##' \item \eqn{\sigma^2_{j}=\sigma^2 k^2_j}: a variance specific to each repetition. This is achieved using a baseline standard deviation parameter (\eqn{\sigma}) and a multiplier parameter (\eqn{k_j}) for each repetition \eqn{j}, with \eqn{k_1=1}.
##' \item \eqn{\rho_{j,j'}=0}: no correlation
##' }
##' With 4 repetitions (\eqn{(j,j') \in \{1,2,3,4\}^2} and \eqn{j\neq j'}) this leads to the following residuals variance-covariance matrix:
##' \deqn{ \begin{bmatrix} \sigma^2 & 0 & 0 & 0 \\ 0 & k_2^2\sigma^2 & 0 & 0 \\ 0 & 0 & k_3^2\sigma^2 & 0 \\ 0 & 0 & 0 & k^2_4\sigma^2 \end{bmatrix}}
##' It is possibly to stratify the covariance pattern on a categorical variable (e.g. sex) using a formula such as \code{sex~1}. This will lead to:
##' \deqn{ \begin{bmatrix} \sigma_{\text{female}}^2 & 0 & 0 & 0 \\ 0 & k_{2:\text{female}}^2\sigma_{\text{female}}^2 & 0 & 0 \\ 0 & 0 & k_{3:\text{female}}^2\sigma_{\text{female}}^2 & 0 \\ 0 & 0 & 0 & k^2_{4:\text{female}}\sigma_{\text{female}}^2 \end{bmatrix}
##' \text{and}
##' \begin{bmatrix} \sigma_{\text{male}}^2 & 0 & 0 & 0 \\ 0 & k_{2:\text{male}}^2\sigma_{\text{male}}^2 & 0 & 0 \\ 0 & 0 & k_{3:\text{male}}^2\sigma_{\text{male}}^2 & 0 \\ 0 & 0 & 0 & k^2_{4:\text{male}}\sigma_{\text{male}}^2 \end{bmatrix}
##' }
##' Instead of stratifying one can also make the multiplier parameter specific to arbitrary (categorical) variables, e.g.  \code{~sex}:
##' \deqn{ \begin{bmatrix} \sigma^2 & 0 & 0 & 0 \\ 0 & k_{\text{female}:2}^2\sigma^2 & 0 & 0 \\ 0 & 0 & k_{\text{female}:3}^2\sigma^2 & 0 \\ 0 & 0 & 0 & k^2_{\text{female}:4}\sigma^2 \end{bmatrix}
##' \text{and}
##' \begin{bmatrix} \sigma^2 & 0 & 0 & 0 \\ 0 & k_{\text{male}:2}^2\sigma^2 & 0 & 0 \\ 0 & 0 & k_{\text{male}:3}^2\sigma^2 & 0 \\ 0 & 0 & 0 & k^2_{\text{male}:4}\sigma^2 \end{bmatrix}
##' }
##' which is just a re-parametrisation the stratified structure since, by default, the variance is taken dependent of a time effect.
##' Removing the time effect by setting \code{var.time} to \code{NULL}) leads to:
##' \deqn{ \begin{bmatrix} \sigma^2 & 0 & 0 & 0 \\ 0 & \sigma^2 & 0 & 0 \\ 0 & 0 & \sigma^2 & 0 \\ 0 & 0 & 0 & \sigma^2 \end{bmatrix}
##' \text{and}
##' \begin{bmatrix} \sigma^2 k_{\text{male}}^2 & 0 & 0 & 0 \\ 0 & \sigma^2 k_{\text{male}}^2 & 0 & 0 \\ 0 & 0 & \sigma^2 k_{\text{male}}^2 & 0 \\ 0 & 0 & 0 & \sigma^2 k_{\text{male}}^2 \end{bmatrix}
##' }
##' 
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @keywords multivariate
##' 
##' @examples
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$sex <- factor(as.numeric(gastricbypassL$id) %% 2,
##'                              labels = c("female", "male"))
##'
##' ## default: repetition specific variance
##' eIND.default <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = IND, data = gastricbypassL)
##' sigma(eIND.default)
##' paramIND.default <- model.tables(eIND.default, effects = "variance")
##' paramIND.default
##' paramIND.default[1,"estimate"]^2 * c(1,paramIND.default[-1,"estimate"]^2)
##' 
##' ## stratified: repetition and strata specific variance parameter
##' eIND.strata <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                    structure = IND(sex~1), data = gastricbypassL)
##' sigma(eIND.strata)
##' model.tables(eIND.strata, effects = "variance")
##' 
##' ## same as
##' eIND.reg <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = IND(~sex), data = gastricbypassL)
##' sigma(eIND.reg) 
##' model.tables(eIND.reg, effects = "variance") 
##'
##' ## to not only use sex dependent variance
##' eIND.reg0 <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = IND(~sex, var.time = NULL), data = gastricbypassL)
##' sigma(eIND.reg0)
##' model.tables(eIND.reg0, effects = "variance") 
##' }
##' 

## * IND (independence, code)
##' @export
IND <- function(formula, var.time){

    ## ** check input
    ## *** formula
    if(inherits(formula,"formula")){
        stop("Argument \'formula\' should be a formula. \n")
    }

    ## *** var.time
    if(!missing(var.time)){
       if(any(inherits(var.time,"formula"))){ ## possible user mistake IND(~1,~1) instead of InD(list(~1,~1))
           stop("Argument \'var.time\' should not be a formula. \n",
                "Consider using a list to collect the formula for the variance and correlation structure. \n")
       }else if(!is.null(var.time)){
           stop("Argument \'var.time\' should either be missing or NULL to enable or disable automatically repetition-specific variance parameters. \n")
       }
    }

    ## ** identify strata variable and covariate variables
    outCov <- .formulaStructure(formula, add.X = NULL, strata.X = FALSE, correlation = FALSE)

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = NA,
                                  variance = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  correlation = NA,
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula.var,
                               correlation = NULL),
                class = c(var = "IND", cor = NA))

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("IND",class(out))
    return(out)
}


## * CS (compound symmetry, documentation)
##' @title Compound Symmetry Structure
##' @description Variance-covariance structure where the variance and correlation of the residuals is constant within covariate levels.
##' In presence of a categorical covariate varying within cluster, the correlation structure is structured in blocks. For instance with 2 categories: \itemize{
##' \item \eqn{A}: pairs of observations both at level 0. 
##' \item \eqn{B}: pairs of observations both at level 1. 
##' \item \eqn{C}: pairs of observations, one with level 0 and the other with level 1.
##' }
##'
##' @param formula [formula or list of 2 formula] variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side).
##' @param var.time [character] time variable relative to which the correlation structure is constructed when \code{cross} is \code{"TOEPLITZ"}, \code{"DUN"}, or \code{"UN"}.
##' @param twin [logical] in presence of a covariate varying within-cluster, should block \eqn{A} and \eqn{B} be identical?
##' \code{TRUE} is analogous to nested random effects.
##' @param cross [character] in presence of a covariate varying within-cluster, structure for block \eqn{C}: \code{"ID"}, \code{"CS"}, \code{"TOEPLITZ"}, or \code{"UN"}.
##' \code{"ID"} is analogous to cross random effects.
##'
##' @details A typical formula would be \code{~1}. It will model: \itemize{
##' \item \eqn{\sigma^2_{j}=\sigma^2}: a variance constant over repetitions.
##' \item \eqn{\rho_{j,j'}=\rho}: a correlation constant over repetitions.
##' }
##' With 4 repetitions (\eqn{(j,j') \in \{1,2,3,4\}^2} and \eqn{j\neq j'}) this leads to the following residuals variance-covariance matrix:
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho\sigma^2 & \rho\sigma^2 & \rho\sigma^2 \\ \rho\sigma^2 & \sigma^2 & \rho\sigma^2 & \rho\sigma^2 \\ \rho\sigma^2 & \rho\sigma^2 & \sigma^2 & \rho\sigma^2 \\ \rho\sigma^2 & \rho\sigma^2 & \rho\sigma^2 & \sigma^2 \end{bmatrix}}
##' With positive correlation, this is a reparametrisation of a random intercept model. \cr
##' 
##' \bold{Covariate constant within-cluster}: it is possibly to stratify the covariance pattern on a categorical variable (e.g. sex) using a formula such as \code{sex~1}. This will lead to:
##' \deqn{ \begin{bmatrix} \sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 \\ \rho_{\text{female}}\sigma_{\text{female}}^2 & \sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 \\ \rho_{\text{female}}\sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 & \sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 \\ \rho_{\text{female}}\sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 & \rho_{\text{female}}\sigma_{\text{female}}^2 & \sigma_{\text{female}}^2 \end{bmatrix}
##' \text{and}
##' \begin{bmatrix} \sigma_{\text{male}}^2 & \rho_{\text{male}} \sigma^2_{\text{male}} & \rho_{\text{male}} \sigma^2_{\text{male}} & \rho_{\text{male}} \sigma^2_{\text{male}} \\ \rho_{\text{male}} \sigma^2_{\text{male}} & \sigma_{\text{male}}^2 & \rho_{\text{male}} \sigma^2_{\text{male}} & \rho_{\text{male}} \sigma^2_{\text{male}} \\ \rho_{\text{male}} \sigma^2_{\text{male}} & \rho_{\text{male}} \sigma^2_{\text{male}} & \sigma_{\text{male}}^2 & \rho_{\text{male}} \sigma^2_{\text{male}} \\ \rho_{\text{male}} \sigma^2_{\text{male}} & \rho_{\text{male}} \sigma^2_{\text{male}} & \rho_{\text{male}} \sigma^2_{\text{male}} & \sigma_{\text{male}}^2 \end{bmatrix}
##' } 
##' By default, the argument \code{twin=TRUE} meaning that \code{~sex} will only make the correlation sex dependent:
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 \\ \rho_{\text{female}}\sigma^2 & \sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 \\ \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \sigma^2 & \rho_{\text{female}}\sigma^2 \\ \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \sigma^2 \end{bmatrix}
##' \text{and}
##' \begin{bmatrix} \sigma^2  & \rho_{\text{male}} \sigma^2  & \rho_{\text{male}} \sigma^2  & \rho_{\text{male}} \sigma^2  \\ \rho_{\text{male}} \sigma^2  & \sigma^2  & \rho_{\text{male}} \sigma^2  & \rho_{\text{male}} \sigma^2  \\ \rho_{\text{male}} \sigma^2  & \rho_{\text{male}} \sigma^2  & \sigma^2  & \rho_{\text{male}} \sigma^2  \\ \rho_{\text{male}} \sigma^2  & \rho_{\text{male}} \sigma^2  & \rho_{\text{male}} \sigma^2  & \sigma^2  \end{bmatrix}
##' }
##' To obtain a sex-dependent variance one can either set the argument \code{twin=FALSE} or specify a separate formula for the variance and correlation using a list: \code{list(~sex,~sex)}
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 \\ \rho_{\text{female}}\sigma^2 & \sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 \\ \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \sigma^2 & \rho_{\text{female}}\sigma^2 \\ \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \rho_{\text{female}}\sigma^2 & \sigma^2 \end{bmatrix}
##' \text{and}
##' \begin{bmatrix} \sigma^2 k_{\text{male}}^2 & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} \\ \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \sigma^2 k_{\text{male}}^2 & \rho_{\text{male}} \sigma^2 k_{\text{male}}^2 & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} \\ \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \sigma^2 k^2_{\text{male}} & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} \\ \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \rho_{\text{male}} \sigma^2 k^2_{\text{male}} & \sigma^2 k^2_{\text{male}} \end{bmatrix}
##' }
##' This is a reparametrization of the stratified CS structure. \cr \cr
##'
##' \bold{Covariate varying within-cluster}: it is not possible to stratify on covariates whose value is not constant across repetition within clusters. But one can make the variance and correlation dependent on such variable. By default (\code{twin=TRUE}) two correlation parameters are used: one when the variable values are same between the two repetitions (\eqn{\rho_W}) and another when the variable values differ between the two repetitions (\eqn{\rho_B}), e.g. if baseline refers to the firt two repetitions then using the formula \code{~baseline} will model:
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho_W\sigma^2 & \rho_B\sigma^2 & \rho_B\sigma^2 \\ \rho_W\sigma^2 & \sigma^2 & \rho_B\sigma^2 & \rho_B\sigma^2 \\ \rho_B\sigma^2 & \rho_B\sigma^2 & \sigma^2 & \rho_W\sigma^2 \\ \rho_B\sigma^2 & \rho_B\sigma^2 & \rho_W\sigma^2 & \sigma^2 \end{bmatrix}}
##' With positive correlation, this is reparametrisation of a nested random effect model. \cr
##'
##' Setting \code{twin=FALSE} enables to obtain a correlation parameter specific to each variable value and each difference in variable value:
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho_{\text{baseline}}\sigma^2 & \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} & \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} \\ \rho_{\text{baseline}}\sigma^2 & \sigma^2 & \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} & \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} \\ \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} & \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} & \sigma^2 k^2_{\text{follow-up}} & \rho_{\text{follow-up}}\sigma^2 k^2_{\text{follow-up}} \\ \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} & \rho_{\text{baseline,follow-up}}\sigma^2 k_{\text{follow-up}} & \rho_{\text{follow-up}}\sigma^2 k^2_{\text{follow-up}} & \sigma^2 k^2_{\text{follow-up}} \end{bmatrix}} 
##'
##' Setting \code{cross="ID"} only considers correlation between observations with identical covariate values. Consider clusters of two subjects and modeling: within-subject correlation and within-time correlation for subjects from the same cluster but assuming independence between observations measured at different timepoints on different subjects:
##' \deqn{ \begin{bmatrix}
##' \sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{time}}\sigma^2  & 0 & 0 & 0 \\
##' \rho_{\text{subject}}\sigma^2 & \sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & 0 & \rho_{\text{time}}\sigma^2 & 0 & 0 \\
##' \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \sigma^2 & \rho_{\text{subject}}\sigma^2 & 0 & 0  & \rho_{\text{time}}\sigma^2 & 0 \\
##' \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \sigma^2 & 0 & 0 & 0 & \rho_{\text{time}}\sigma^2 \\
##' \rho_{\text{time}}\sigma^2 & 0 & 0 & 0 & \sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2\\
##' 0 & \rho_{\text{time}}\sigma^2 & 0 & 0 & \rho_{\text{subject}}\sigma^2 & \sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 \\
##' 0 & 0  & \rho_{\text{time}}\sigma^2 & 0 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \sigma^2 & \rho_{\text{subject}}\sigma^2\\
##' 0 & 0 & 0 & \rho_{\text{time}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \rho_{\text{subject}}\sigma^2 & \sigma^2 \\
##' \end{bmatrix}} 
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @keywords multivariate
##' 
##' @examples
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$sex <- factor(as.numeric(gastricbypassL$id) %% 2,
##'                              labels = c("female", "male"))
##' gastricbypassL$baseline <- factor(gastricbypassL$time < 0,
##'                              labels = c("baseline", "follow-up"))
##' gastricbypassL$family <- paste0("F",(as.numeric(gastricbypassL$id)-1) %/% 2)
##'
##' #### default: constant variance and correlation ####
##' eCS.default <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = CS, data = gastricbypassL)
##' sigma(eCS.default)
##' paramCS.default <- model.tables(eCS.default, effects = "all")
##' paramCS.default
##' paramCS.default["sigma","estimate"]^2 * paramCS.default["rho(id)","estimate"]
##'
##' #### Covariate constant within-cluster #### 
##' ## stratified: strata specific variance and correlation
##' eCS.strata <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = CS(sex~1), data = gastricbypassL)
##' sigma(eCS.strata)
##' model.tables(eCS.strata, effects = "all")
##' 
##' ## covariate: constant variance and sex-dependent correlation
##' eCS.sexCor <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                   structure = CS(~sex), data = gastricbypassL)
##' sigma(eCS.sexCor)
##' 
##' ## covariate: sex-dependent variance and correlation
##' eCS.sex <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                structure = CS(~sex, twin = TRUE), data = gastricbypassL)
##' sigma(eCS.sex)
##' 
##' eCS.sex2 <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                structure = CS(list(~sex, ~sex)), data = gastricbypassL)
##' sigma(eCS.sex2)
##'
##' #### Covariate varying within-cluster ####
##'
##' ## twin: within (rho id/basline) and between (rho id) correlation
##' eBCS.default <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = CS(~baseline), data = gastricbypassL)
##' sigma(eBCS.default)
##' model.tables(eBCS.default, effects = "all")
##' 
##' ## heterogeneous: within (rho id/basline) and between (rho id) correlation
##' eBCS.hetero <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                     structure = CS(~baseline, twin = FALSE), data = gastricbypassL)
##' sigma(eBCS.hetero)
##' model.tables(eBCS.hetero, effects = "all")
##' 
##' #### Cross compound symmetry ####
##' eCS.cross <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                  structure = CS(~baseline, cross = "ID"),
##'                  data = gastricbypassL)
##' sigma(eCS.cross)
##' model.tables(eCS.cross, effects = "all")
##' 
##' ## artificially create a two level structure family-id
##' ## index family member: first, second
##' gastricbypassL <- gastricbypassL[order(gastricbypassL$id),]
##' gastricbypassL$member <- repetition(~id|family, data = gastricbypassL,
##'                                     type = "consecutive", label.rep = c("I","II"))
##' gastricbypassL$visit.member <- interaction(gastricbypassL[,c("visit","member")])
##'
##' ## constant correlation within subject, constant correlation within time,
##' ## no correlation within family for different subjects at different times
##' eCS.family <- lmm(glucagonAUC ~ visit, repetition = ~visit.member|family,
##'                   structure = CS(~id+visit, cross = "ID"),
##'                   data = gastricbypassL)
##' sigma(eCS.family)
##' model.tables(eCS.family, effects = "all")
##' 
##' eHCS.family <- lmm(glucagonAUC ~ visit, repetition = ~visit.member|family,
##'                   structure = CS(~member+visit, twin = TRUE, cross = "ID"),
##'                   data = gastricbypassL) ## lack of convergence
##' sigma(eHCS.family)
##' model.tables(eHCS.family, effects = "all")
##' }
##' 

## * CS (compound symmetry, code)
##' @export
CS <- function(formula, var.time, twin = TRUE, cross = "CS"){
browser()
    ## ** normalize input
    ## *** formula
    if(inherits(formula,"formula") && !(is.list(formula) && length(formula)= 2 && all(sapply(formula,inherits,"formula")))){
        stop("Argument \'formula\' should be a formula or a list of 2 formula. \n")
    }
    if(is.list(formula) && length(formula)==2){
        if(!is.null(names(formula))){
            names(formula) <- c("variance","correlation")
        }else if(any(duplicated(names(formula)))){
            stop("Elements in the list specified by argument \'formula\' should not have duplicated names. \n",
                 "Current names: \"",paste(names(formula),collapse ="\", \""),"\".\n")
        }else if(any(names(formula) %in% c("variance","correlation") == FALSE)){
            stop("Elements in the list specified by argument \'formula\' should be named \"variance\" or \"correlation\". \n",
                 "Current names: \"",paste(names(formula),collapse ="\", \""),"\".\n")
        }
    }

    ## *** var.time
    if(!missing(var.time) && inherits(var.time,"formula")){ ## possible user mistake CS(~1,~1) instead of CS(list(~1,~1))
        stop("Argument \'var.time\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## *** twin
    if(is.null(twin)){ ## default value for twin depending on how much information is provided by the user (one or two formula)
        twin <- inherits(formula,"formula") ## only automatically TRUE when variance model not specified 
    }

    ## *** cross
    cross <- match.arg(cross, choices = c("ID","CS","TOEPLTIZ","DUN","UN"))
    browser()
    if(all.vars(formula))
    if(!missing(var.time)){
        if(cross %in% c("ID","CS")){
            warning("Argument \'var.time\' ignored with CS structure unless it has argument \'cross\' being \"TOEPLTIZ\", \"DUN\", \"UN\". \n")
        }else if(cross %in% c("TOEPLTIZ","DUN","UN") && is.null(var.time)){
            stop("Argument \'var.time\' cannot be NULL with argument \'cross\' being \"TOEPLTIZ\", \"DUN\", \"UN\". \n",
                 "Should be missing or a character indicating the time variable. \n")
        }
    }
    
    ## ** normalize input
    ## *** enforce twin
    if(inherits(formula,"formula")){
        if(twin){
            formula <- list(variance = stats::update(formula,".~1"), ## remove covariate in variance model
                            correlation = formula)
        }else{
            formula <- list(variance = formula,
                            correlation = formula)
        }
    }

    ## *** add time
    if(!missing(var.time) & cross %in% c("TOEPLTIZ","DUN","UN")){
        if(length(var.time)>1)

            if(var.time %in% all.vars(formula$correlation)){
                stop("The formula relative to the correlation model in argument \'formula\' already contains the time variable \"",var.time,"\" specified in argument \'var.time\'. \n")
            }else if(attr(var.time,"original") %in% all.vars(formula[[2]])){
                stop("The formula relative to the correlation model in argument \'formula\' already contains the time variable \"",attr(var.time,"original"),"\" specified in argument \'var.time\'. \n")
            }

        if(!is.null(attr(var.time,"original"))){
            formula$correlation <- stats::update(formula$correlation, paste0("~.+",attr(var.time,"original")))
        }else{
            formula$correlation <- stats::update(formula$correlation, paste0("~.+",var.time))
        }
    }
    
    ## ** from formula to covariate names
    outCov <- .formulaStructure(formula, add.X = NULL, strata.X = FALSE, correlation = TRUE)
    if(length(outCov$X.cor)==0){
        twin <- TRUE
    }

    if(!missing(var.cluster) && inherits(var.cluster,"formula")){ ## possible user mistake CS(~1,~1) instead of CS(list(~1,~1))
        stop("Argument \'var.cluster\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }
    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = NA,
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                twin = twin,
                cross = cross,
                class = c(var = "IND", cor = "CS"))

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("CS",class(out))
    return(out)
}

## * RE (random effect, documentation)
##' @title Random Effect Structure
##' @description Variance-covariance structure parametrized via random effects.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side).##' 
##' @param var.time [character] time variable.
##' @param ranef [list] characteristics of the random effects
##'
##' @details A typical formula would be \code{~1}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @keywords multivariate
##' 
##' @examples
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$sex <- factor(as.numeric(gastricbypassL$id) %% 2,
##'                              labels = c("female", "male"))
##' gastricbypassL$baseline <- factor(gastricbypassL$time < 0,
##'                              labels = c("baseline", "follow-up"))
##' gastricbypassL$family <- paste0("F",(as.numeric(gastricbypassL$id)-1) %/% 2)
##'
##' #### Random intercept ####
##' eRI <- lmm(glucagonAUC ~ visit + (1|id), data = gastricbypassL)
##' sigma(eRI)
##' nlme::ranef(eRI, effects = "variance")
##' model.tables(eRI, effects = "all")
##'
##' #### Nested random effects ####
##' eNRI <- lmm(glucagonAUC ~ visit + (1|id/baseline), data = gastricbypassL)
##' sigma(eNRI)
##' ## nlme::ranef(eNRI, effects = "variance")
##' model.tables(eNRI, effects = "all")
##'
##' ## more than two blocks
##' data("sleepL", package = "LMMstar")
##' eNRI3 <- lmm(signal.34 ~ day + (1|id/day), data = sleepL)
##' sigma(eNRI3)
##' 
##' #### Cross random effects ####
##' gastricbypassL <- gastricbypassL[order(gastricbypassL$id),]
##' gastricbypassL$member <- repetition(~id|family, data = gastricbypassL,
##'                                     type = "consecutive", label.rep = c("I","II"))
##' gastricbypassL$visit.member <- interaction(gastricbypassL[,c("visit","member")])
##' 
##' eCRI <- lmm(glucagonAUC ~ visit + (1|id) + (1|visit),
##'             repetition = ~visit.member|family, data = gastricbypassL)
##' ## nlme::ranef(eCRI, effects = "variance")
##' sigma(eCRI)
##' model.tables(eCRI, effects = "all")
##'
##' data(Penicillin, package = "lme4")
##' eCRI.2 <- lmm(diameter ~ 1 + (1|plate) + (1|sample), data = Penicillin)
##' nlme::ranef(eCRI.2, effects = "variance")
##' model.tables(eCRI.2, effects = "all")
##' }

## * RE (random effect, code)
##' @export
RE <- function(formula, var.time, ranef = NULL){

    ## ** normalize input
    ## ranef
    test.ranef <- is.null(ranef)

    ## var.cluster
    if(!missing(var.cluster) && inherits(var.cluster,"formula")){ ## possible user mistake RE(~1,~1) instead of RE(list(~1,~1))
        stop("Argument \'var.cluster\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' must inherits from formula. \n")
    }
    detail.formula <- formula2var(formula)

    ## convert ~1|id to ~(1|id)
    if(detail.formula$special=="repetition"){
        if(length(detail.formula$vars$repetition)>2 || length(detail.formula$terms$all)!=1){
            stop("Something when wrong when identifying the random effects. \n",
                 "Consider wrapping the random effect term in the formula by parentheses, e.g. ~ (1|id). \n")
        }
        newformula <- stats::as.formula(paste0(paste(detail.formula$vars$response, collapse = "+"),"~(",detail.formula$terms$repetition,")"))
        detail.formula <- formula2var(newformula)
    }

    if(detail.formula$special=="none"){
        if(length(detail.formula$vars$regressor)>0){
            if(length(detail.formula$vars$response)>0){
                stop("The strata variable in argument \'formula\' should be specified on the left or right hand side (not both). \n")
            }
            var.strata <- detail.formula$vars$regressor
        }else{
            var.strata <- detail.formula$vars$response
        }
    }else if(detail.formula$special=="ranef"){
        var.strata <- detail.formula$vars$response        
    }else{
        stop("Incorrect argument \'formula\' for structure RE. \n",
             "Should be something like ~strata or strata ~ (1|id).")
    }

    ## ** prepare variance structure
    if(length(var.strata)>0){
        ff.var <- stats::as.formula(paste0(var.strata,"~1"))
    }else{
        ff.var <- ~1
    }

    ## ** prepare correlation / random effect structure
    if(test.ranef && detail.formula$special=="ranef"){
        ranef <- detail.formula$ranef
        if(missing(var.cluster) && ranef$cross == FALSE){
            var.cluster <- ranef$cluster
            attr(var.cluster,"original") <- ranef$cluster
        }
    }
    
    if(is.null(ranef) || (ranef$cross == FALSE && ranef$nested == FALSE)){ ## no random effect or random intercept
        ff.cor <- ff.var   
        group <- NULL
    }else if(ranef$cross == FALSE && ranef$nested == TRUE){
        Ovar.cluster <- attr(var.cluster,"original")
        if(!is.null(Ovar.cluster) && Ovar.cluster %in% ranef$vars){
            ## remove cluster level from formula 
            ff.cor <- stats::as.formula(paste0(var.strata,"~",paste0(setdiff(ranef$vars,Ovar.cluster), collapse = "+")))
        }else{
            stop("Something went wrong when identifying the nested random effects. \n",
                 "Cluster variable ",attr(var.cluster,"original")," not found among the variables defining the random effects \"",paste(ranef$vars, collapse = "\" \""),"\".\n")
        }
    }else if(ranef$cross == TRUE && ranef$nested == FALSE){
        ## crossed random effect: typically each random effect has been specified (e.g. (1|id) and (1|scan))
        ## but not the overall level (e.g. hospital where specific id and scans are performed)
        ff.cor <- stats::as.formula(paste0(var.strata,"~",paste0(ranef$vars, collapse = "+")))
    }else if(ranef$cross == TRUE && ranef$nested == TRUE){
        stop("Random effect structre with both nested and cross random effect not (yet!) possible. \n")
    }

    ## ** create structure
    out <- CS(list(variance = ff.var, correlation = ff.cor),
              var.cluster = var.cluster, var.time = var.time,
              twin = TRUE,
              cross = ranef$cross)
    out$ranef <- ranef

    ## ** export
    out$call <- match.call()
    out$class <- "RE"
    class(out) <- append("RE",class(out))
    return(out)
}


## * TOEPLITZ (Toeplitz, documentation)
##' @title Toeplitz Structure
##' @description Variance-covariance structure where the correlation depends on time elapsed between two repetitions.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side). 
##' @param var.cluster [character] cluster variable.
##' @param type [character] degree of flexibility of the covariance structure within covariate.
##' For a binary covariate the correlation structure can be decomposed into 4 blocks: within level 0 of the covariate (\eqn{A}),
##' C within level 1 (\eqn{A}), and B between level 0 and 1 (\eqn{A}). The following parametrisations are available:
##' \itemize{ 
##' \item \code{"UN"}: unstructured matrix except for the diagonal elements of C which are constrained to be equal.
##' \item \code{"LAG"}: Toeplitz structure within A, B, and C, i.e. correlation specific to each time lag and covariate level.
##' \item \code{"CS"}: block-specific value except for C which has a different value for its diagonal elements.
##' }
##' @param var.time [character] time variable.
##' @param add.time Should the default formula (i.e. when \code{NULL}) contain a time effect.
##'
##' @details \bold{formula}: there can only be at most one covariate for the correlation structure.
##' A typical formula would be \code{~1}, indicating a variance constant over time and a correlation specific to each gap time.
##'
##' 
##' @return An object of class \code{TOEPLITZ} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @keywords multivariate
##' 
##' @details A typical formula would be \code{~1}. It will model: \itemize{
##' \item \eqn{\sigma^2_{j}=\sigma^2 k^2_j}: a variance specific to each repetition. This is achieved using a baseline standard deviation parameter (\eqn{\sigma}) and a multiplier parameter (\eqn{k_j}) for each repetition \eqn{j}, with \eqn{k_1=1}.
##' \item \eqn{\rho_{j,j'}=\rho_{|j-j'|}}: a correlation specific to the difference in index between repetitions.
##' }
##' With 4 repetitions (\eqn{(j,j') \in \{1,2,3,4\}^2} and \eqn{j\neq j'}) this leads to the following residuals variance-covariance matrix:
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho_1\sigma^2 k_2 & \rho_2\sigma^2 k_3 & \rho_3\sigma^2 k_4 \\
##'                        \rho_1\sigma^2 k_2 & \sigma^2 k_2^2 & \rho_1\sigma^2 k_2 k_3 & \rho_2\sigma^2 k_2 k_4 \\
##'                        \rho_2\sigma^2 k_3 & \rho_1\sigma^2 k_2 k_3 & \sigma^2 k_3^2 & \rho_1\sigma^2 k_3 k_4 \\
##'                        \rho_3\sigma^2 k_4 & \rho_2\sigma^2 k_2 k_4 & \rho_1 \sigma^2 k_3 k_4 & \sigma^2 k_4^2
##' \end{bmatrix}}
##' WARNING: the variance-covariance pattern is not invariant under reparametrisation of the repetition variable. For instance with \eqn{(j,j') \in \{0,1,2,5\}^2}:
##' \deqn{ \begin{bmatrix} \sigma^2 & \rho_1\sigma^2 k_1 & \rho_2\sigma^2 k_2 & \rho_5\sigma^2 k_5 \\
##'                        \rho_1\sigma^2 k_1 & \sigma^2 k_1^2 & \rho_1\sigma^2 k_1 k_2 & \rho_4\sigma^2 k_1 k_5 \\
##'                        \rho_2\sigma^2 k_2 & \rho_1\sigma^2 k_1 k_2 & \sigma^2 k_2^2 & \rho_3\sigma^2 k_2 k_5 \\
##'                        \rho_5\sigma^2 k_5 & \rho_4\sigma^2 k_1 k_5 & \rho_3 \sigma^2 k_2 k_5 & \sigma^2 k_5^2
##' \end{bmatrix}}
##' 
##' \bold{Covariate constant within-cluster}: it is possibly to stratify the covariance pattern on a categorical variable (e.g. sex) specifying it in the left hand side of the formula, e.g. \code{sex~1}.
##' \deqn{ \begin{bmatrix} \sigma_{\text{female}}^2 & \rho_{1\text{female}}\sigma_{\text{female}}^2 k_{2,\text{female}} & \rho_{2,\text{female}}\sigma_{\text{female}}^2 k_{3,\text{female}} & \rho_{3,\text{female}}\sigma_{\text{female}}^2 k_{4,\text{female}} \\
##'                        \rho_{1,\text{female}}\sigma_{\text{female}}^2 k_{2,\text{female}} & \sigma_{\text{female}}^2 k_{2,\text{female}}^2 & \rho_{1,\text{female}}\sigma_{\text{female}}^2 k_{2,\text{female}} k_{3,\text{female}} & \rho_{2,\text{female}}\sigma_{\text{female}}^2 k_{2,\text{female}} k_{4,\text{female}} \\
##'                        \rho_{2,\text{female}}\sigma_{\text{female}}^2 k_{3,\text{female}} & \rho_{1,\text{female}}\sigma_{\text{female}}^2 k_{2,\text{female}} k_{3,\text{female}} & \sigma_{\text{female}}^2 k_{3,\text{female}}^2 & \rho_{1,\text{female}}\sigma_{\text{female}}^2 k_{3,\text{female}} k_{4,\text{female}} \\
##'                        \rho_{3,\text{female}}\sigma_{\text{female}}^2 k_{4,\text{female}} & \rho_{2,\text{female}}\sigma_{\text{female}}^2 k_{2,\text{female}} k_{4,\text{female}} & \rho_{1,\text{female}} \sigma_{\text{female}}^2 k_{3,\text{female}} k_{4,\text{female}} & \sigma_{\text{female}}^2 k_{4,\text{female}}^2
##' \end{bmatrix}
##' \text{ and }
##' \begin{bmatrix} \sigma_{\text{male}}^2 & \rho_{1\text{male}}\sigma_{\text{male}}^2 k_{2,\text{male}} & \rho_{2,\text{male}}\sigma_{\text{male}}^2 k_{3,\text{male}} & \rho_{3,\text{male}}\sigma_{\text{male}}^2 k_{4,\text{male}} \\
##'                        \rho_{1,\text{male}}\sigma_{\text{male}}^2 k_{2,\text{male}} & \sigma_{\text{male}}^2 k_{2,\text{male}}^2 & \rho_{1,\text{male}}\sigma_{\text{male}}^2 k_{2,\text{male}} k_{3,\text{male}} & \rho_{2,\text{male}}\sigma_{\text{male}}^2 k_{2,\text{male}} k_{4,\text{male}} \\
##'                        \rho_{2,\text{male}}\sigma_{\text{male}}^2 k_{3,\text{male}} & \rho_{1,\text{male}}\sigma_{\text{male}}^2 k_{2,\text{male}} k_{3,\text{male}} & \sigma_{\text{male}}^2 k_{3,\text{male}}^2 & \rho_{1,\text{male}}\sigma_{\text{male}}^2 k_{3,\text{male}} k_{4,\text{male}} \\
##'                        \rho_{3,\text{male}}\sigma_{\text{male}}^2 k_{4,\text{male}} & \rho_{2,\text{male}}\sigma_{\text{male}}^2 k_{2,\text{male}} k_{4,\text{male}} & \rho_{1,\text{male}} \sigma_{\text{male}}^2 k_{3,\text{male}} k_{4,\text{male}} & \sigma_{\text{male}}^2 k_{4,\text{male}}^2
##' \end{bmatrix}
##' }
##' 
##' \bold{Covariate varying within-cluster}: it is not possible to stratify on covariates whose value is not constant across repetition within clusters.
##' Such variables should be cateorical covariate and specified on the right-hand side of the formula, e.g. \code{~baseline} to divide the variance-covariance matrix by block.
##' Each block can be a compound symmetry matrix (\code{type=="LAG"}), e.g. consider 3 baseline and 3 follow-up observations per subject:
##' \deqn{ \begin{bmatrix} \sigma^2           & \rho_X \sigma^2   & \rho_X \sigma^2       & \rho \sigma^2 k_Y \\
##'                        \rho_X\sigma^2     & \sigma^2          & \rho \sigma^2 k_Y     & \rho \sigma^2 k_Y \\
##'                        \rho \sigma^2 k_Y  & \rho \sigma^2 k_Y & \sigma^2 k_Y^2        & \rho_Y \sigma^2 k_Y^2 \\
##'                        \rho \sigma^2 k_Y  & \rho \sigma^2 k_Y & \rho_Y \sigma^2 k_Y^2 & \sigma^2 k_Y^2
##' \end{bmatrix}}
##' or compound symmetry matrix (\code{type=="CS"}), e.g.:
##' 
##' or unstructure (\code{type=="UN"}), e.g.:
##' 
##' @examples
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$sex <- factor(as.numeric(gastricbypassL$id) %% 2,
##'                              labels = c("female", "male"))
##' gastricbypassL$baseline <- factor(gastricbypassL$time < 0,
##'                              labels = c("baseline", "follow-up"))
##' 
##' data("sleepL", package = "LMMstar")
##' sleepL$sex <- factor(as.numeric(sleepL$id) %% 2,
##'                      labels = c("female", "male"))
##' sleepL$baseline <- factor(sleepL$day==1, c(TRUE,FALSE))
##' sleepL$visit <- repetition(~1|id, data = sleepL, format = "numeric")
##' sleepL <- sleepL[sleepL$visit<=6,]
##' 
##' #### default ####
##' ## time as integer: 1, 2, 3 (traditional TOEPLITZ structure)
##' eTOEint.default <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                        structure = "TOEPLITZ", data = gastricbypassL)
##' sigma(eTOEint.default)
##' cov2cor(sigma(eTOEint.default))
##' model.tables(eTOEint.default, effects = "all")
##' 
##' ## time as numeric: here -13 vs -1 and -1 vs 1 are one index appart
##' ##                  but lead to different correlation coefficients
##' ##                  rho(12) and rho(2)
##' eTOEnum.default <- lmm(glucagonAUC ~ visit, repetition = ~time|id,
##'                     structure = "TOEPLITZ", data = gastricbypassL)
##' sigma(eTOEnum.default)
##' cov2cor(sigma(eTOEnum.default))
##' model.tables(eTOEnum.default, effects = "all")
##' 
##' ## another example with more timepoints
##' eTOEint6.default <- lmm(signal.98 ~ baseline, repetition = ~visit|id,
##'                         structure = "TOEPLITZ", data = sleepL)
##' sigma(eTOEint6.default)
##' cov2cor(sigma(eTOEint6.default))
##' model.tables(eTOEint6.default, effects = "all")
##' 
##' #### Covariate constant within-cluster ####
##' eTOE.strata <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                     structure = TOEPLITZ(sex~1), data = gastricbypassL)
##' sigma(eTOE.strata)
##' model.tables(eTOE.strata, effects = "all")
##'
##' #### Covariate varying within-cluster ####
##' 
##' ## block compound symmetry structure
##' eTOE.CS <- lmm(signal.98 ~ baseline, repetition = ~visit|id,
##'                structure = TOEPLITZ(~baseline, type = "CS"), data = sleepL)
##' sigma(eTOE.CS)
##' cov2cor(sigma(eTOE.CS))
##' paramTOE.CS <- model.tables(eTOE.CS, effects = "all")
##' paramTOE.CS
##' paramTOE.CS["rho1","estimate"]*paramTOE.CS["sigma","estimate"]^2
##' paramTOE.CS["rho2","estimate"]*prod(paramTOE.CS[c("sigma","k.FALSE"),"estimate"]^2)
##' prod(paramTOE.CS[c("rho(1,2,dt=1)","k.FALSE"),"estimate"])*paramTOE.CS[c("sigma"),"estimate"]^2
##'
##' ## block toeplitz structure
##' eTOE.LAG <- lmm(signal.98 ~ baseline, repetition = ~visit|id,
##'                structure = TOEPLITZ(~baseline, type = "LAG"), data = sleepL)
##' sigma(eTOE.LAG, cluster = 1)
##' cov2cor(sigma(eTOE.LAG, cluster = 1))
##' model.tables(eTOE.LAG, effects = "all")
##'
##' ## block unstructured structure
##' eTOE.UN <- lmm(glucagonAUC ~ visit, repetition = ~visit|id,
##'                structure = TOEPLITZ(~baseline, type = "UN"), data = gastricbypassL)
##' sigma(eTOE.UN)
##' cov2cor(sigma(eTOE.UN))
##' model.tables(eTOE.UN, effects = "all")
##'
##' 
##' }

## * TOEPLITZ (Toeplitz, code)
##' @export
TOEPLITZ <- function(formula, var.cluster, type = "LAG", var.time, add.time){

    ## ** normalize input
    if(is.null(type)){
        type <- "LAG"
        toeplitz.block <- FALSE
    }else{
        type <- match.arg(type, c("UN","LAG","CS"))
        toeplitz.block <- NULL
    }
    if(!missing(add.time)){

        if(is.character(add.time)){            
            if(type == "UN" || type == "LAG"){
                add.X <- list(variance = add.time,
                              correlation = add.time)
            }else if(type == "CS"){
                if(length(add.time)>1){
                    add.X <- list(variance = utils::head(add.time,1),
                                  correlation = add.time)
                }else{
                    add.X <- list(variance = NULL,
                                  correlation = add.time)
                }
            }
        }else if(add.time){
            if(type == "UN" || type == "LAG"){
                add.X <- list(variance = var.time,
                              correlation = var.time)
            }else if(type == "CS"){
                if(length(var.time)>1){
                    add.X <- list(variance = utils::head(var.time,1),
                                  correlation = var.time)
                }else{
                    add.X <- list(variance = NULL,
                                  correlation = var.time)
                }
            }
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }

    if(!missing(formula) && inherits(formula,"formula")){
        ## if(!is.null(add.X$variance)){
        ##     formula <- list(variance = ~1,
        ##                     correlation = formula)
        ## }else{        
        formula <- list(variance = formula,
                        correlation = formula)
        ## }
    }
    
    if(!missing(formula)){
        outCov <- .formulaStructure(formula, add.X = add.X, strata.X = FALSE, correlation = TRUE)

        ## NOTE: if length(outCov$X.cor)==0 i.e. not time has (yet) be defined, no error is output since time will be automatically added later on
        if(length(outCov$X.cor)>2){
            stop("TOEPLITZ covariance structure does not support more than 2 covariates for the correlation structure. \n")
        }else if(is.null(toeplitz.block)){
            toeplitz.block <- length(outCov$X.cor) > 1
        }
    }else{
        outCov <- NULL
    }

    if(!missing(var.cluster) && inherits(var.cluster,"formula")){ ## possible user mistake TOEPLITZ(~1,~1) instead of TOEPLITZ(list(~1,~1))
        stop("Argument \'var.cluster\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                type = type,
                block = toeplitz.block,
                class = c(var = "IND", cor = "TOEPLITZ"))

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("TOEPLITZ",class(out))
    return(out)
}

## * UN (unstructured, documentation)
##' @title Unstructured Structure 
##' @description Variance-covariance structure where the residuals have time-specific variance and correlation.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the covariance structure.
##' @param var.time [character] time variable.
##'
##' @details A typical formula would be \code{~1}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##'
##' @return An object of class \code{UN} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @keywords multivariate
##' 
##' @examples
##' UN(NULL, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(~sex, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(sex ~ 1, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(list(~sex,~1), var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(list(sex~age,sex~1), var.cluster = "id", var.time = "time", add.time = TRUE)
##' 

## * UN (unstructured, code)
##' @export
UN <- function(formula, var.time){

    ## ** normalize input
    if(!missing(add.time)){
        if(is.character(add.time)){
            add.X <- list(variance = add.time,
                          correlation = add.time)
        }else if(add.time){
            add.X <- list(variance = var.time,
                          correlation = var.time)
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }

    outCov <- .formulaStructure(formula, add.X = add.X, strata.X = inherits(formula,"formula"), correlation = 2) ## 2 for fully stratified structure

    if(!missing(var.cluster) && inherits(var.cluster,"formula")){ ## possible user mistake UN(~1,~1) instead of UN(list(~1,~1))
        stop("Argument \'var.cluster\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = I(list(outCov$X.var)),
                                  cor = I(list(outCov$X.cor)),
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                c(var = "IND", cor = "UN"))

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * LV (latent variable, documentation)

## ## * EXP (exponential)
## ##' @title Exponential Structure
## ##' @description Variance-covariance structure where the residuals have a correlation decreasing exponentially,
## ##' Can be stratified on a categorical variable.
## ##'
## ##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
## ##' and variables influencing the residual variance and correlation (right hand side).
## ##' @param var.cluster [character] cluster variable.
## ##' @param var.time [character] time variable.
## ##' @param nugget [logical] whether a nugget effect is present.
## ##' @param add.time not used.
## ##'
## ##' @details A typical formula would be \code{~1}, indicating a variance constant over time and correlation with exponential decrease over time.
## ##'
## ##' Inspired from \code{nlme::corExp} where if \eqn{K} denotes the nugget effect and \eqn{\rho} the time effect,
## ##' the correlation between two observations with a time gap \eqn{dt} is \eqn{exp(-\rho dt)} when no nugget effect is present and \eqn{(1-K) exp(-\rho dt)} when a nugget effect is assumed. 
## ##'
## ##' @return An object of class \code{EXP} that can be passed to the argument \code{structure} of the \code{lmm} function.
## ##' 
## ##' @keywords multivariate
## ##' 
## ##' @examples
## ##' EXP(var.cluster = "id", var.time = "time", add.time = TRUE)
## ##' EXP(~space, var.cluster = "id", var.time = "time", add.time = TRUE)
## ##' EXP(list(~space,~space), var.cluster = "id", var.time = "time", add.time = TRUE)
## ##' 
## ##' @export
## EXP <- function(formula, var.cluster, var.time, nugget = FALSE, add.time){

##     if(missing(formula) || is.null(formula)){
##         outCov <- .formulaStructure(list(~1,stats::as.formula(paste0("~",var.time))), heterogeneous = nugget)
##     }else if(is.list(formula)){
##         outCov <- .formulaStructure(formula, heterogeneous = nugget)
##     }else if(!missing(add.time) && (is.character(add.time) || identical(add.time,TRUE)) && length(all.vars(stats::update(formula,0~.)))==0){
##         if(is.character(add.time)){
##             var.time <- add.time
##         }
##         if(attr(stats::terms(formula),"response")==1){ # with strata
##             ff <- stats::as.formula(paste0(all.vars(formula),"~",var.time))
##         }else{
##             ff <- stats::as.formula(paste0("~",var.time))
##         }
##         outCov <- .formulaStructure(list(formula,ff), heterogeneous = nugget)
##     }else{
##         if(attr(stats::terms(formula),"response")==1){ # with strata
##             outCov <- .formulaStructure(list(stats::as.formula(paste0(all.vars(formula),"~1")),formula), heterogeneous = nugget)
##         }else{
##             outCov <- .formulaStructure(list(~1,formula), heterogeneous = nugget)
##         }
##     }

##     out <- list(call = match.call(),
##                 name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
##                                   strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
##                                   time = if(!missing(var.time)){var.time}else{NA},
##                                   var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
##                                   cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
##                                   stringsAsFactors = FALSE),
##                 formula = list(var = outCov$formula.var,
##                                cor = outCov$formula.cor),
##                 heterogeneous = nugget,
##                 class = "EXP")

##     ## export
##     class(out) <- append("structure",class(out))
##     class(out) <- append("EXP",class(out))
##     return(out)
## }


## * CUSTOM (user-specified, documentation)
##' @title Custom Structure
##' @description Variance-covariance structure specified by the user.
##' 
##' @param formula formula indicating variables influencing the residual variance and correlation (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param FCT.sigma [function] take as argument the model parameters, time, and design matrix.
##' Output the vector of residuals standard deviations.
##' @param dFCT.sigma [list of vectors] list whose elements are the first derivative of argument \code{FCT.sigma}. 
##' @param d2FCT.sigma [list of vectors] list whose elements are the second derivative of argument \code{FCT.sigma} (no cross-terms).
##' @param init.sigma [numeric vector] initial value for the variance parameters.
##' @param FCT.rho [function] take as argument the model parameters, time, and design matrix.
##' Output the matrix of residuals correlation.
##' @param dFCT.rho [list of matrices] list whose elements are the first derivative of argument \code{FCT.rho}. 
##' @param d2FCT.rho [list of matrices] list whose elements are the second derivative of argument \code{FCT.rho} (no cross-terms).
##' @param init.rho [numeric vector] initial value for the correlation parameters.
##' @param add.time not used.
##'
##' @return An object of class \code{CUSTOM} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @keywords multivariate
##' 
##' @examples
##' 
##' ## Compound symmetry structure
##' CUSTOM(~1,
##'        FCT.sigma = function(p,n.time,X){rep(p,n.time)},
##'        init.sigma = c("sigma"=1),
##'        dFCT.sigma = function(p,n.time,X){list(sigma = rep(1,n.time))},  
##'        d2FCT.sigma = function(p,n.time,X){list(sigma = rep(0,n.time))},  
##'        FCT.rho = function(p,n.time,X){
##'            matrix(p,n.time,n.time)+diag(1-p,n.time,n.time)
##'        },
##'        init.rho = c("rho"=0.5),
##'        dFCT.rho = function(p,n.time,X){
##'             list(rho = matrix(1,n.time,n.time)-diag(1,n.time,n.time))
##'        },
##'        d2FCT.rho = function(p,n.time,X){list(rho = matrix(0,n.time,n.time))}
##' )
##' 
##' ## 2 block structure
##' rho.2block <- function(p,n.time,X){
##'    rho <- matrix(0, nrow = n.time, ncol = n.time)
##'    rho[1,2] <- rho[2,1] <- rho[4,5] <- rho[5,4] <- p["rho1"]
##'    rho[1,3] <- rho[3,1] <- rho[4,6] <- rho[6,4] <- p["rho2"]
##'    rho[2,3] <- rho[3,2] <- rho[5,6] <- rho[6,5] <- p["rho3"]
##'    rho[4:6,1:3] <- rho[1:3,4:6] <- p["rho4"]
##'    return(rho)
##' }
##' drho.2block <- function(p,n.time,X){
##'    drho <- list(rho1 = matrix(0, nrow = n.time, ncol = n.time),
##'                 rho2 = matrix(0, nrow = n.time, ncol = n.time),
##'                 rho3 = matrix(0, nrow = n.time, ncol = n.time),
##'                 rho4 = matrix(0, nrow = n.time, ncol = n.time))   
##'    drho$rho1[1,2] <- drho$rho1[2,1] <- drho$rho1[4,5] <- drho$rho1[5,4] <- 1
##'    drho$rho2[1,3] <- drho$rho2[3,1] <- drho$rho2[4,6] <- drho$rho2[6,4] <- 1
##'    drho$rho3[2,3] <- drho$rho3[3,2] <- drho$rho3[5,6] <- drho$rho3[6,5] <- 1
##'    drho$rho4[4:6,1:3] <- drho$rho4[1:3,4:6] <- 1
##'    return(drho)
##' }
##' d2rho.2block <- function(p,n.time,X){
##'    d2rho <- list(rho1 = matrix(0, nrow = n.time, ncol = n.time),
##'                  rho2 = matrix(0, nrow = n.time, ncol = n.time),
##'                  rho3 = matrix(0, nrow = n.time, ncol = n.time),
##'                  rho4 = matrix(0, nrow = n.time, ncol = n.time))   
##'    return(d2rho)
##' }
##'
##' CUSTOM(~variable,
##'        FCT.sigma = function(p,n.time,X){rep(p,n.time)},
##'        dFCT.sigma = function(p,n.time,X){list(sigma=rep(1,n.time))},
##'        d2FCT.sigma = function(p,n.time,X){list(sigma=rep(0,n.time))},
##'        init.sigma = c("sigma"=1),
##'        FCT.rho = rho.2block,
##'        dFCT.rho = drho.2block,
##'        d2FCT.rho = d2rho.2block,
##'        init.rho = c("rho1"=0.25,"rho2"=0.25,"rho3"=0.25,"rho4"=0.25))
##' 

## * CUSTOM (user-specified, code)
##' @export
CUSTOM <- function(formula, var.time,
                   FCT.sigma, dFCT.sigma = NULL, d2FCT.sigma = NULL, init.sigma,
                   FCT.rho, dFCT.rho = NULL, d2FCT.rho = NULL, init.rho, add.time){

    if(is.null(formula)){
        outCov <- .formulaStructure(~1, add.X = NULL, strata.X = FALSE, correlation = TRUE)
    }else{
        outCov <- .formulaStructure(formula, add.X = NULL, strata.X = FALSE, correlation = TRUE)
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = NA,
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                format.cortime = NA,
                class = c(var = "CUSTOM", cor = "CUSTOM"))

    ## param
    ## if(!is.na(out$name$strata)){
    ##     stop("CUSTOM structure cannot (yet?) handle strata. \n")
    ## }
    if(!missing(FCT.sigma)){
        if(any(names(formals(FCT.sigma)) %in% c("p","n.time","X","...") == FALSE)){
            stop("Incorrect argument in \'FCT.sigma\': can be \"p\", \"n.time\", or \"X\". \n")
        }
        if(!is.null(dFCT.sigma) && any(names(formals(dFCT.sigma)) %in% c("p","n.time","X","...") == FALSE)){
            stop("Incorrect argument in \'dFCT.sigma\': can be \"p\", \"n.time\", or \"X\". \n")
        }
        if(!is.null(d2FCT.sigma) && any(names(formals(d2FCT.sigma)) %in% c("p","n.time","X","...") == FALSE)){
            stop("Incorrect argument in \'d2FCT.sigma\': can be \"p\", \"n.time\", or \"X\". \n")
        }
        out$FCT.sigma <- FCT.sigma
        out$dFCT.sigma <- dFCT.sigma
        out$d2FCT.sigma <- d2FCT.sigma
        out$init.sigma <- init.sigma
        if(length(out$init.sigma)>0 && is.null(names(out$init.sigma))){
            names(out$init.sigma) <- paste0("sigma",1:length(out$init.sigma))
        }
    }else{
        out$FCT.sigma <- function(p,time,X){rep(p,length(time))}
        out$dFCT.sigma <- function(p,time,X){list(rep(1,length(time)))}
        out$d2FCT.sigma <- function(p,time,X){list(rep(0,length(time)))}
        out$init.sigma <- c(sigma = NA)        
    }
    if(!missing(FCT.rho)){
        if(any(names(formals(FCT.rho)) %in% c("p","n.time","X","...") == FALSE)){
            stop("Incorrect argument in \'FCT.rho\': can be \"p\", \"n.time\", or \"X\". \n")
        }
        if(!is.null(dFCT.rho) && any(names(formals(dFCT.rho)) %in% c("p","n.time","X","...") == FALSE)){
            stop("Incorrect argument in \'dFCT.rho\': can be \"p\", \"n.time\", or \"X\". \n")
        }
        if(!is.null(d2FCT.rho) && any(names(formals(d2FCT.rho)) %in% c("p","n.time","X","...") == FALSE)){
            stop("Incorrect argument in \'d2FCT.rho\': can be \"p\", \"n.time\", or \"X\". \n")
        }
        
        out$FCT.rho <- FCT.rho
        out$dFCT.rho <- dFCT.rho
        out$d2FCT.rho <- d2FCT.rho
        out$init.rho <- init.rho
        if(length(init.rho)>0 && is.null(names(out$init.rho))){
            names(out$init.rho) <- paste0("rho",1:length(out$init.rho))
        }
    }else{
        out$FCT.rho <- NULL
        out$dFCT.rho <- NULL
        out$d2FCT.rho <- NULL
        out$init.rho <- NULL
        out$formula$cor <- NULL
    }

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("CUSTOM",class(out))
    return(out)
}

## * helper
## ** .formulaStructure
##' @title Extract Variable From Formula For VCOV Structure
##' @description Extract the variables from the variance and correlation formula to be used to initialize the variance-covariance structure.
##' @noRd
##'
##' @param formula A formula or a list of two formulas.
##' @param add.X [character vector] additional covariates to be added to the variance and correlation structure.
##' @param strata.X [logical] all covariates should be used to stratify the variance and/or correlation structure.
##' @param correlation [logical] should a correlation structure be output? Otherwise no correlation parameter is considered.
##' If greater than 1 then all interaction terms are considered when having multiple variables, e.g. time:sex instead of time + sex.
##' 
##' @keywords internal
.formulaStructure <- function(formula, add.X, strata.X, correlation){

    ## ** normalize to formula format
    if(is.null(formula)){
        formula <- ~1
    }

    class.formula <- inherits(formula,"formula")
    if(class.formula){
        detail.formula0 <- formula2var(formula)

        formula <- list(variance = formula,
                        correlation = formula)
        detail.formula <- list(variance = detail.formula0,
                               correlation = detail.formula0)
    }else if(is.list(formula) && length(formula)==2 && all(sapply(formula,inherits,"formula"))){

        if(is.null(names(formula))){
            names(formula) <- c("variance","correlation")
        }else if(all(names(formula) %in% c("var","cor"))){
            names(formula)[names(formula)=="var"] <- "variance"
            names(formula)[names(formula)=="cor"] <- "correlation"
        }else if(all(names(formula) %in% c("variance","correlation"))){
            formula <- formula[c("variance","correlation")]
        }else{
            stop("Incorrect names associated to the argument \'formula\' for the residual variance-covariance structure. \n",
                 "Should be \"variance\" and \"correlation\" (or \"var\" and \"cor\"). \n")
        }
        detail.formula <- lapply(formula, formula2var)
    }else{
        stop("Incorrect argument \'formula\': should be a formula or a list of 2 formula (var, cor).\n")
    }

    if(is.character(add.X)){
        add.X <- list(variance = add.X,
                      correlation = add.X)
    }else if(is.list(add.X) && length(add.X)==2){        
        if(is.null(names(add.X))){
            names(add.X) <- c("variance","correlation")
        }else if(all(names(add.X) %in% c("var","cor"))){
            names(add.X)[names(add.X)=="var"] <- "variance"
            names(add.X)[names(add.X)=="cor"] <- "correlation"
        }else if(all(names(add.X) %in% c("variance","correlation"))){
            add.X <- add.X[c("variance","correlation")]
        }else{
            stop("Incorrect names associated to the argument \'add.X\' for the residual variance-covariance structure. \n",
                 "Should be \"variance\" and \"correlation\" (or \"var\" and \"cor\"). \n")
        }
    }else if(!is.null(add.X)){
        stop("Incorrect argument \'add.X\': should be a character or a list with 2 elements (var, cor).\n")
    }

    ## ** type
    if(any(sapply(detail.formula, function(iDetail){iDetail$special})!="none")){
        stop("Incorrect argument \'formula\': there should be no special character, e.g. no |. \n")
    }

    ## ** left hand side
    ls.var.strata <- lapply(detail.formula, function(iDetail){iDetail$vars$response})
    if(!identical(ls.var.strata[[1]],ls.var.strata[[2]])){
        stop("Incorrect argument \'formula\': strata variable differ between the correlation and variance structure. \n")
    }
    var.strata <- ls.var.strata[[1]]
    if(length(var.strata)==0){
        var.strata <- NULL
    }else if(length(var.strata)>1){
        stop("There should be at most one strata variable. \n")
    }

    ## ** right hand side    
    ls.var.X <- list(variance = unique(c(detail.formula$variance$vars$regressor, add.X$variance)),
                     correlation = unique(c(detail.formula$correlation$vars$regressor, add.X$correlation)))
    test.interaction <- sapply(formula, function(iF){
        any(attr(stats::delete.response(stats::terms(iF)),"order")>1)
    })
    if(any(test.interaction)){
        stop("Does not handle interactions in the formula. \n")
    }

    ## right hand side variables define strata
    if(length(var.strata)==0 && strata.X && length(unlist(ls.var.X))>0){
        if(length(ls.var.X$variance)>0 && length(ls.var.X$correlation)>0){
            if(!identical(detail.formula$variance$vars$regressor,detail.formula$correlation$vars$regressor)){
                stop("The strata variable should not differ between the variance and the correlation structure. \n")
            }else{
                var.strata <- detail.formula$variance$vars$regressor
                ls.var.X$variance <- setdiff(ls.var.X$variance,detail.formula$variance$vars$regressor)
                ls.var.X$correlation <- setdiff(ls.var.X$correlation,detail.formula$variance$vars$regressor)
            }
        }else if(length(ls.var.X$variance)>0){
            var.strata <- detail.formula$variance$vars$regressor
            ls.var.X$variance <- setdiff(ls.var.X$variance,detail.formula$variance$vars$regressor)                
        }else if(length(ls.var.X$correlation)>0){
            var.strata <- detail.formula$correlation$vars$regressor
            ls.var.X$correlation <- setdiff(ls.var.X$correlation,detail.formula$correlation$vars$regressor)                
        }

    }

    ## ** combine left and right hand side
    ## *** variance
    if(length(ls.var.X$variance)==0){
        X.var <- NULL
        if(length(var.strata)==0){            
            if(attr(stats::terms(formula$variance),"intercept")){
                formula.var <- ~1
            }else{
                formula.var <- ~0
            }
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
        X.var <- unname(ls.var.X$variance)
        if(length(var.strata)==0){
            formula.var <- stats::as.formula(paste("~",paste(ls.var.X$variance,collapse=":")))
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata,"+ ",paste(c(ls.var.X$variance, var.strata),collapse=":")))            
        }
    }

    ## *** correlation
    if(length(ls.var.X$correlation)==0){
        X.cor <- NULL
        if(correlation==FALSE){
            formula.cor <- NULL
        }else if(length(var.strata)==0){
            formula.cor <- ~1
        }else{
            formula.cor <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
        X.cor <- unname(ls.var.X$correlation)
        if(correlation==FALSE){
            formula.cor <- NULL
            if(!class.formula){
                warning("Variable(s) \"",paste(ls.var.X$correlation, collapse = "\", \""),"\" in the correlation structure are ignored. \n")
            }
        }else if(length(var.strata)==0){
            if(correlation>1){
                formula.cor <- stats::as.formula(paste("~0+",paste(ls.var.X$correlation,collapse=":")))
            }else{
                formula.cor <- stats::as.formula(paste("~0+",paste(ls.var.X$correlation,collapse="+")))
            }
        }else{
            if(correlation>1){
                formula.cor <- stats::as.formula(paste("~0 + ",paste(c(ls.var.X$correlation, var.strata), collapse=":")))
            }else{
                formula.cor <- stats::as.formula(paste("~0 + (",paste(ls.var.X$correlation, collapse="+"),"):",var.strata))
                ## NOTE: wrapping the non-strata variable in () ensures proper ordering of the column names (i.e. strata variable at the end)
                ##' df <- data.frame(day = c(1, 1, 2, 2, 1, 1, 2, 2),
                ##'                  sex = c("1", "1", "1", "1", "0", "0", "0", "0"),
                ##'                  session = c(1, 2, 3, 4, 1, 2, 3, 4))
                ##' X <- stats::model.matrix(~0 + session:sex + day:sex, df)
                ##' colnames(X) ## "session:sex0" "session:sex1" "sex0:day" "sex1:day"
            }
        }
    }

    ## ** export
    out <- list(strata = var.strata,
                X.var = X.var,
                X.cor = X.cor,
                formula.var = formula.var,
                formula.cor = formula.cor
                )
    return(out)
}

##----------------------------------------------------------------------
### structure.R ends here
