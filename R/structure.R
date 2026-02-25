### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: feb 18 2026 (17:23) 
##           By: Brice Ozenne
##     Update #: 1897
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
ID <- function(formula = ~1){

    ## ** check input
    ## *** formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' should be a formula. \n")
    }
    
    ## ** identify strata variable and covariate variables
    outCov <- .formulaStructure(list(variance = formula, correlation = NULL, correlation.cross = NULL), correlation = FALSE)
    if(length(outCov$name$variance)>0 || length(outCov$name$correlation)>0){ 
        stop("Argument \'formula\' should not contain variable(s) on the right-hand side of the formula. \n",
             "Structure ID can only be stratified (left hand side of the formula). \n")
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  ordering = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = NA,
                                  variance = if(length(outCov$name$variance)>0){I(list(outCov$name$variance))}else{NA},
                                  correlation = NA,
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula$variance,
                               correlation = NULL,
                               correlation.cross = NULL),
                class = c(variance = "ID", correlation = NA, correlation.cross = NA))

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
##' @param heterogeneous [logical] should a repetition-specific variance model be considered?
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
IND <- function(formula = ~1, heterogeneous =  TRUE){

    ## ** check input
    ## *** formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' should be a formula. \n")
    }

    ## *** var.time
    if(any(inherits(heterogeneous,"formula"))){ ## possible user mistake IND(~1,~1) instead of InD(list(~1,~1))
        stop("Argument \'heterogeneous\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }else if(length(heterogeneous) == 1 && heterogeneous %in% 0:1){
        update.time <- c(variance = heterogeneous, correlation = FALSE, correlation.cross = FALSE) 
    }else{
        stop("Argument \'heterogeneous\' should have length 1 and be TRUE or FALSE. \n")
    }

    ## ** identify strata variable and covariate variables
    outCov <- .formulaStructure(list(variance = formula, correlation = NULL, correlation.cross = NULL), correlation = FALSE)
    if(update.time["variance"]){outCov$formula$variance <- formula$variance}

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  ordering = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = NA,
                                  variance = if(length(outCov$name$variance)>0){I(list(outCov$name$variance))}else{NA},
                                  correlation = NA,                                  
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula$variance,
                               correlation = NULL,
                               correlation.cross = NULL),
                class = c(variance = "IND", correlation = NA, correlation.cross = NA))
    attr(out$formula,"update.time") <- update.time

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
##' @param formula [formula or list of 2 formula] left hand side: strata variable for the variance and correlation structure. 
##' right hand side: covariate for the variance structure (multiplicative effect) and correlation structure (blocks).
##' @param heterogeneous [logical] should a repetition-specific variance model be considered?
##' @param cross [character] in presence of a covariate varying within-cluster, structure for block \eqn{C}: \code{"ID"}, \code{"CS"}, \code{"TOEPLITZ"}, or \code{"UN"}.
##' @param twin [logical] in presence of a covariate varying within-cluster, should block \eqn{A} and \eqn{B} be identical?
##' @param time [character] time variable relative to which the correlation structure is constructed when \code{cross} is \code{"TOEPLITZ"}, \code{"DUN"}, or \code{"UN"}.
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
CS <- function(formula = ~1, heterogeneous = FALSE, cross = "CS", twin = TRUE, time){

    ## ** normalize input    
    ## *** formula
    if(inherits(formula,"formula")){
        formula <- list(variance = update(formula,~1), ## by default CS structure ignores covariate for the variance model
                        correlation = formula,
                        correlation.cross = formula)
    }else{
        formula <- .initFormulaStructure(formula)    
    }

    ## *** heterogeneous
    if(!missing(heterogeneous) && inherits(heterogeneous,"formula")){ ## possible user mistake CS(~1,~1) instead of CS(list(~1,~1))
        stop("Argument \'heterogeneous\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## *** cross
    cross <- match.arg(cross, choices = c("ID","CS","TOEPLITZ","DUN","UN"))
    

    ## *** twin
    if(length(twin)!=1){
        stop("Argument \'twin\' should have length 1. \n")
    }
    if(!is.logical(twin) && twin!=0 && twin!=1){
        stop("Argument \'twin\' be binary: TRUE or FALSE. \n")
    }

    ## *** time
    if(!missing(time)){
        if(!is.character(time) || length(time)!=1){
            stop("Argument \'time\' should be a character of length 1 when not missing. \n")
        }
    }

    ## ** normalize input
    
    ## *** add time
    if(!missing(time)){
        update.time <- c(variance = FALSE, correlation = FALSE, correlation.cross = FALSE)
        ls.time <- list(variance = NULL, correlation = NULL, correlation.cross = NULL)
        if(heterogeneous){ls.time$variance <- time}
        if(cross %in% c("TOEPLITZ","DUN","UN")){ls.time$correlation.cross <- time}    
    }else{
        time <- NA
        update.time <- c(variance = heterogeneous,
                         correlation = FALSE,
                         correlation.cross = cross %in% c("TOEPLITZ", "DUN", "UN")) 
        ls.time <- list(variance = NULL, correlation = NULL, correlation.cross = NULL)        
    }

    ## ** from formula to covariate names
    outCov <- .formulaStructure(formula, ls.time = ls.time, correlation = TRUE)
    if(length(setdiff(outCov$name$correlation, c(outCov$strata,time)))==0){ ## no block if no covariates
        cross <- NA
        outCov$name$correlation <- all.vars(outCov$formula$correlation)
        outCov$formula$correlation.cross <- NULL
    }else if(length(setdiff(outCov$name$correlation, time))>1){
        warning("CS covariance structure has not been developped for more 1 covariate defining blocks in the correlation structure. \n")
    }

    if(update.time["variance"]){outCov$formula$variance <- formula$variance}
    if(update.time["correlation"]){outCov$formula$correlation <- formula$correlation}
    if(update.time["correlation.cross"]){outCov$formula$correlation.cross <- formula$correlation.cross}

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  ordering = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = time,
                                  variance = if(length(outCov$name$variance)>0){I(list(outCov$name$variance))}else{NA},
                                  correlation = if(length(outCov$name$correlation)>0){I(list(outCov$name$correlation))}else{NA},                                  
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula$variance,
                               correlation = outCov$formula$correlation,
                               correlation.cross = outCov$formula$correlation.cross),
                twin = twin,
                class = c(variance = "IND", correlation = "CS", correlation.cross = cross))
    attr(out$formula,"update.time") <- update.time

    ## ** export
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
##' @param time [character] time variable.
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
RE <- function(formula, time, ranef = NULL){

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
              var.cluster = var.cluster, time = time,
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
##' In presence of a categorical covariate varying within cluster, the correlation structure is structured in blocks. For instance with 2 categories: \itemize{
##' \item \eqn{A}: pairs of observations both at level 0. 
##' \item \eqn{B}: pairs of observations both at level 1. 
##' \item \eqn{C}: pairs of observations, one with level 0 and the other with level 1.
##' }
##'
##' @param formula [formula or list of 2 formula] left hand side: strata variable for the variance and correlation structure. 
##' right hand side: covariate for the variance structure (multiplicative effect) and correlation structure (blocks).
##' @param heterogeneous [logical] should a repetition-specific variance model be considered?
##' @param cross [character] in presence of a covariate varying within-cluster, structure for block \eqn{C}: \code{"ID"}, \code{"CS"}, \code{"TOEPLITZ"}, or \code{"UN"}.
##' @param twin [logical] in presence of a covariate varying within-cluster, should block \eqn{A} and \eqn{B} be identical?
##' @param time [character] time variable relative to which the correlation structure is constructed when \code{cross} is \code{"TOEPLITZ"}, \code{"DUN"}, or \code{"UN"}.
##' 
##' @details \bold{formula}: a typical formula would be \code{~1}, indicating a variance constant over time and a correlation specific to each gap time.
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
TOEPLITZ <- function(formula = ~1, heterogeneous = TRUE, cross = "TOEPLITZ", twin = TRUE, time){

    ## ** normalize input    
    ## *** formula
    formula <- .initFormulaStructure(formula)    

    ## *** heterogeneous
    if(!missing(heterogeneous) && inherits(heterogeneous,"formula")){ ## possible user mistake TOEPLITZ(~1,~1) instead of TOEPLITZ(list(~1,~1))
        stop("Argument \'heterogeneous\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## *** cross
    cross <- match.arg(cross, choices = c("ID","CS","TOEPLITZ","DUN","UN"))
    
    ## *** twin
    if(length(twin)!=1){
        stop("Argument \'twin\' should have length 1. \n")
    }
    if(!is.logical(twin) && twin!=0 && twin!=1){
        stop("Argument \'twin\' be binary: TRUE or FALSE. \n")
    }

    ## *** time
    if(!missing(time)){
        if(!is.character(time) || length(time)!=1){
            stop("Argument \'time\' should be a character of length 1 when not missing. \n")
        }
    }

    ## ** normalize input

    ## *** add time
    if(!missing(time)){
        update.time <- c(variance = FALSE, correlation = FALSE, correlation.cross = FALSE)
        ls.time <- list(variance = NULL, correlation = time, correlation.cross = NULL)
        if(heterogeneous){ls.time$variance <- time}
        if(cross %in% c("TOEPLITZ","DUN","UN")){ls.time$correlation.cross <- time}    
    }else{
        time <- NA
        update.time <- c(variance = heterogeneous,
                         correlation = TRUE,
                         correlation.cross = cross %in% c("TOEPLITZ", "DUN", "UN")) 
        ls.time <- list(variance = NULL, correlation = NULL, correlation.cross = NULL)        
    }
    
    ## ** from formula to covariate names
    outCov <- .formulaStructure(formula, ls.time = ls.time, correlation = TRUE)
    if(length(setdiff(outCov$name$correlation, time))==0){ ## no block if no covariates
        cross <- NA
        outCov$name$correlation <- all.vars(outCov$formula$correlation)
        outCov$formula$correlation.cross <- NULL
    }else if(length(setdiff(outCov$name$correlation, time))>1){
        warning("TOEPLITZ covariance structure has not been developped for more 1 covariate defining blocks in the correlation structure. \n")
    }

    if(update.time["variance"]){outCov$formula$variance <- formula$variance}
    if(update.time["correlation"]){outCov$formula$correlation <- formula$correlation}
    if(update.time["correlation.cross"]){outCov$formula$correlation.cross <- formula$correlation.cross}

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  ordering = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = time,
                                  variance = if(length(outCov$name$variance)>0){I(list(outCov$name$variance))}else{NA},
                                  correlation = if(length(outCov$name$correlation)>0){I(list(outCov$name$correlation))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula$variance,
                               correlation = outCov$formula$correlation,
                               correlation.cross = outCov$formula$correlation.cross),
                twin = twin,
                class = c(variance = "IND", correlation = "TOEPLITZ", correlation.cross = cross))
    attr(out$formula,"update.time") <- update.time
    
    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("TOEPLITZ",class(out))
    return(out)
}

## * UN (unstructured, documentation)
##' @title Unstructured Structure 
##' @description Variance-covariance structure where the residuals have time-specific variance and correlation.
##' In presence of a categorical covariate varying within cluster, the correlation structure is structured in blocks. For instance with 2 categories: \itemize{
##' \item \eqn{A}: pairs of observations both at level 0. 
##' \item \eqn{B}: pairs of observations both at level 1. 
##' \item \eqn{C}: pairs of observations, one with level 0 and the other with level 1.
##' }
##'
##' @param formula [formula or list of 2 formula] left hand side: strata variable for the variance and correlation structure. 
##' right hand side: covariate for the variance structure (multiplicative effect) and correlation structure (blocks).
##' @param heterogeneous [logical] should a repetition-specific variance model be considered?
##' @param cross [character] in presence of a covariate varying within-cluster, structure for block \eqn{C}: \code{"ID"}, \code{"CS"}, \code{"TOEPLITZ"}.
##' @param twin [logical] in presence of a covariate varying within-cluster, should block \eqn{A} and \eqn{B} be identical?
##' @param time [character] time variable relative to which the correlation structure is constructed when \code{cross} is \code{"TOEPLITZ"}, \code{"DUN"}.
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
UN <- function(formula = ~1, heterogeneous = TRUE, cross = "DUN", twin = TRUE, time){

    ## ** normalize input    
    ## *** formula
    formula <- .initFormulaStructure(formula)    

    ## *** heterogeneous
    if(!missing(heterogeneous) && inherits(heterogeneous,"formula")){ ## possible user mistake TOEPLITZ(~1,~1) instead of TOEPLITZ(list(~1,~1))
        stop("Argument \'heterogeneous\' should not be a formula. \n",
             "Consider using a list to collect the formula for the variance and correlation structure. \n")
    }

    ## *** cross
    cross <- match.arg(cross, choices = c("ID","CS","TOEPLITZ","DUN"))

    ## *** twin
    if(length(twin)!=1){
        stop("Argument \'twin\' should have length 1. \n")
    }
    if(!is.logical(twin) && twin!=0 && twin!=1){
        stop("Argument \'twin\' be binary: TRUE or FALSE. \n")
    }

    
    ## ** normalize input

    ## *** add time
    if(!missing(time)){
        update.time <- c(variance = FALSE, correlation = FALSE, correlation.cross = FALSE)
        ls.time <- list(variance = NULL, correlation = time, correlation.cross = NULL)
        if(heterogeneous){ls.time$variance <- time}
        if(cross %in% c("TOEPLITZ","DUN","UN")){ls.time$correlation.cross <- time}    
    }else{
        time <- NA
        update.time <- c(variance = heterogeneous,
                         correlation = TRUE,
                         correlation.cross = cross %in% c("TOEPLITZ", "DUN", "UN")) 
        ls.time <- list(variance = NULL, correlation = NULL, correlation.cross = NULL)        
    }
    
    ## ** from formula to covariate names
    outCov <- .formulaStructure(formula, ls.time = ls.time, correlation = TRUE)
    if(length(setdiff(outCov$name$correlation, time))==0){ ## no block if no covariates
        cross <- NA
        outCov$name$correlation <- all.vars(outCov$formula$correlation)
        outCov$formula$correlation.cross <- NULL
    }else if(length(setdiff(outCov$name$correlation, time))>1){
        warning("UN covariance structure has not been developped for more 1 covariate defining blocks in the correlation structure. \n")
    }

    if(update.time["variance"]){outCov$formula$variance <- formula$variance}
    if(update.time["correlation"]){outCov$formula$correlation <- formula$correlation}
    if(update.time["correlation.cross"]){outCov$formula$correlation.cross <- formula$correlation.cross}

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  ordering = NA,
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = time,
                                  variance = if(length(outCov$name$variance)>0){I(list(outCov$name$variance))}else{NA},
                                  correlation = if(length(outCov$name$correlation)>0){I(list(outCov$name$correlation))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula$variance,
                               correlation = outCov$formula$correlation,
                               correlation.cross = outCov$formula$correlation.cross),
                twin = twin,
                class = c(variance = "IND", correlation = "UN", correlation.cross = cross))
    attr(out$formula,"update.time") <- update.time
        
    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * EXP (exponential)
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

## * AR1 (documentation)
## special but famous case of exponential?

## * LV (latent variable, documentation)
## * ANTE (ante-dependence, documentation)


## * CUSTOM (user-specified, documentation)
##' @title Custom Structure
##' @description Variance-covariance structure specified by the user.
##' 
##' @param formula [formula] variables to stratify the residual variance (left hand side, like sex~1).
##' @param var.time not used.
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
CUSTOM <- function(formula = ~1, var.time,
                   FCT.sigma, dFCT.sigma = NULL, d2FCT.sigma = NULL, init.sigma,
                   FCT.rho, dFCT.rho = NULL, d2FCT.rho = NULL, init.rho, add.time){


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
            warning("Argument \'var.time\' ignored with CUSTOM structure. \n")
        }
    }

    ## ** identify strata variable and covariate variables
    if(is.null(formula)){
        outCov <- .formulaStructure(list(variance = ~1, correlation = ~1), correlation = TRUE)
    }else{
        outCov <- .formulaStructure(list(variance = formula, correlation = formula), correlation = TRUE)
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = NA,
                                  ordering = NA,
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = NA,
                                  variance = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  correlation = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(variance = outCov$formula.var,
                               correlation = outCov$formula.cor),
                class = c(var = "CUSTOM", cor = "CUSTOM"))
    attr(out$formula,"update.time") <- FALSE

    ## ** add user-specified structure
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
## ** .initFormulaStructure
.initFormulaStructure <- function(formula){

    if(inherits(formula,"formula")){

        formula <- list(variance = formula,
                        correlation = formula,
                        correlation.cross = formula)
        
    }else if(is.list(formula)){

        if(any(sapply(formula,inherits,"formula")==FALSE)){
            stop("When a list, each element of argument \'formula\' should be a formula. \n")
        }
        
        if(any(duplicated(names(formula)))){
            stop("Elements in the list specified by argument \'formula\' should not have duplicated names. \n",
                 "Current names: \"",paste(names(formula),collapse ="\", \""),"\".\n")
        }
        
        if(length(formula)==2){
            
            if(is.null(names(formula))){
                names(formula) <- c("variance","correlation")
            }else if(any(names(formula) %in% c("variance","correlation") == FALSE)){
                stop("When a list of length 2, elements in the argument \'formula\' should be named \"variance\" or \"correlation\". \n",
                     "Current names: \"",paste(names(formula),collapse ="\", \""),"\".\n")
            }
            formula$correlation.cross <- formula$correlation
            
        }else if(length(formula==3)){
            
            if(is.null(names(formula))){
                names(formula) <- c("variance","correlation","correlation.cross")
            }else if(any(names(formula) %in% c("variance","correlation","correlation.cross") == FALSE)){
                stop("When a list of length 3, elements in the argument \'formula\' should be named \"variance\", \"correlation\", or \"correlation.cross\". \n",
                     "Current names: \"",paste(names(formula),collapse ="\", \""),"\".\n")
            }
            
        }else{
            
            stop("When a list, argument \'formula\' should have length 2 (variance, correlation) or 3 (variance, correlation, correlation.cross). \n")
            
        }
        
    }else{
        
        stop("Argument \'formula\' should be a formula or a list of formula for the variance and correlation structure. \n")
        
    }

    return(formula)
}

## ** .formulaStructure
##' @title Extract Variable From Formula For VCOV Structure
##' @description Extract the variables from the variance and correlation formula to be used to initialize the variance-covariance structure.
##' @noRd
##'
##' @param formula [list] A list with two formulas: one element called variance and a second one called correlation.
##' @param ls.time [list of length 2] additional covariate to be added to the variance and correlation structure.
##' @param correlation [logical] should a correlation structure be output? Otherwise no correlation parameter is considered.
##' If greater than 1 then all interaction terms are considered when having multiple variables, e.g. time:sex instead of time + sex.
##' 
##' @keywords internal
.formulaStructure <- function(formula, ls.time, correlation){

    ## *** add time
    if(!missing(ls.time) && !is.null(ls.time) && any(lengths(ls.time)>0)){

        if(any(lengths(ls.time)>1)){
            stop("Structure cannot handle multiple time variables. \n")
        }

        ## add time to the each structure
        for(iMoment in names(formula)){

            ## skip when no model for this moment
            if(is.null(formula[[iMoment]]) || length(ls.time[[iMoment]])==0){ next }
            ## use original time variable when possible
            if(!is.null(attr(ls.time[[iMoment]],"original"))){
                ls.time[[iMoment]] <- attr(ls.time[[iMoment]],"original")
            }
            ## check that the time variable is not already in the formula
            if(ls.time[[iMoment]] %in% all.vars(formula[[iMoment]])){
                stop("The formula relative to the ",iMoment," model in argument \'formula\' already contains the time variable \"",ls.time,"\" specified in argument \'ls.time\'. \n",
                     ifelse(iMoment=="variance","Consider setting the argument \'heterogeneous\' to FALSE when creating the covariance structure. \n",""))
            }
            formula[[iMoment]] <- stats::update(formula[[iMoment]], paste0("~.+",ls.time[[iMoment]]))
        }
    }

    ## *** move from formula format to characters
    detail.formula <- lapply(formula, formula2var)
    
    ## ** type
    if(any(sapply(detail.formula, function(iDetail){!is.null(iDetail$special) && iDetail$special != "none"}))){
        stop("Incorrect argument \'formula\': there should be no special character, e.g. no |. \n")
    }

    ## ** check same strata across all formulas (left hand side)
    test.strata <- lapply(detail.formula, function(iDetail){sort(iDetail$vars$response)})
    if(sum((!duplicated(test.strata))>1)){
        stop("Incorrect argument \'formula\': strata variable differ between the ",paste(names(detail.formula)[which(!duplicated(test.strata))], collapse=", ")," structures. \n")
    }
    var.strata <- detail.formula$variance$vars$response
    if(length(var.strata)==0){
        var.strata <- NULL
    }else if(length(var.strata)>1){
        stop("There should be at most one strata variable. \n")
    }

    ## ** right hand side
    test.interaction <- sapply(formula, function(iMoment){
        if(!is.null(formula$correlation)){
            return(any(attr(stats::terms(formula$variance),"order")>1))
        }else{
            return(FALSE)
        }
    })
    if(any(test.interaction)){
        stop("Cannot handle interaction(s) in the ",paste(names(detail.formula)[test.interaction>0], collapse=", ")," structure. \n")
    }

    ls.regressor <- lapply(detail.formula, function(iO){iO$vars$regressor})

    ## ** combine strata variable (left hand side) with covariates (right hand side)
    out <- list(strata = var.strata,
                name = list(variance = NULL, correlation = NULL),
                formula = list(variance = NULL, correlation = NULL, correlation.cross = NULL)
                )
    
    ## *** variance
    if(length(ls.regressor$variance)==0){
        if(length(var.strata)==0){            
            if(attr(stats::terms(formula$variance),"intercept")){
                out$formula$variance <- ~1
            }else{
                out$formula$variance <- ~0
            }
        }else{
            out$name$variance <- var.strata
            out$formula$variance <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
        if(length(var.strata)==0){
            out$name$variance <- unname(ls.regressor$variance)
            out$formula$variance <- stats::as.formula(paste("~",paste(ls.regressor$variance,collapse=":")))
        }else{
            out$name$variance <- c(unname(ls.regressor$variance),var.strata)
            out$formula$variance <- stats::as.formula(paste("~0+",var.strata,"+ ",paste(c(ls.regressor$variance, var.strata),collapse=":")))            
        }
    }

    ## *** correlation
    for(iMoment in c("correlation","correlation.cross")){ ## iMoment <- "correlation.cross"
            
        if(correlation && !is.null(formula[[iMoment]]) && length(ls.regressor[[iMoment]])==0){
            if(length(var.strata)==0){
                out$formula[[iMoment]] <- ~1
            }else{
                out$name$correlation <- union(out$name$correlation,var.strata)
                out$formula[[iMoment]] <- stats::as.formula(paste("~0+",var.strata))
            }
        }else if(correlation && !is.null(formula[[iMoment]])){
            if(length(var.strata)==0){
                out$name$correlation <- union(out$name$correlation,unname(ls.regressor[[iMoment]]))
                out$formula[[iMoment]] <- stats::as.formula(paste("~0+",paste(ls.regressor[[iMoment]],collapse="+")))
            }else{
                out$name$correlation <- union(out$name$correlation,c(unname(ls.regressor[[iMoment]]),var.strata))
                out$formula[[iMoment]] <- stats::as.formula(paste("~0 + (",paste(ls.regressor[[iMoment]], collapse="+"),"):",var.strata))
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
    return(out)
}

##----------------------------------------------------------------------
### structure.R ends here
