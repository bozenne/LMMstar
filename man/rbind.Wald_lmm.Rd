% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbind.R
\name{rbind.Wald_lmm}
\alias{rbind.Wald_lmm}
\title{Linear Hypothesis Testing Across Linear Mixed Models}
\usage{
\method{rbind}{Wald_lmm}(model, ..., effects = NULL, rhs = NULL, name = NULL, sep = ": ")
}
\arguments{
\item{model}{a \code{Wald_lmm} object (output of \code{anova} applied to a \code{lmm} object)}

\item{...}{possibly other \code{Wald_lmm} objects}

\item{effects}{[character or numeric matrix] how to combine the left-hand side of the hypotheses.
By default identity matrix but can also be \code{"Dunnett"},  \code{"Tukey"}, or  \code{"Sequen"} (see function \code{multcomp::contrMat} from the multcomp package).}

\item{rhs}{[numeric vector] the right hand side of the hypothesis. Should have the same length as the number of row of argument \code{effects}.}

\item{name}{[character vector or NULL] character used to identify each model in the output.
By default, use the name of the outcome of the model.}

\item{sep}{[character] character used to separate the outcome and the covariate when naming the tests.}
}
\description{
Linear hypothesis testing accross linear mixed model.
}
\details{
WARNING: in presence of measurements from the same cluster across several models,
the influence function is used to estimate the covariance between the model parameters.
This is why the (robust) standard errors may not match the (model-based) standard error from the linear mixed
Setting the argument \code{robust} to \code{FALSE} when calling \code{anova.lmm} will rescale the (robust) standard errors to mimic the original model-based standard errors.
}
\examples{
## simulate data
set.seed(10)
dL <- sampleRem(1e2, n.times = 3, format = "long")

## estimate mixed models
e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, data = dL,
              structure = "CS", df = FALSE)
e.lmm2 <- lmm(Y ~ X1+X8+X9, repetition = ~visit|id, data = dL,
              structure = "CS", df = FALSE)

model.tables(e.lmm1) ## model-based standard errors
model.tables(e.lmm1, robust = TRUE) ## robust standard errors

## select null hypotheses & combine (robust standard errors)
AAA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"))
BBB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"))

ZZZ <- rbind(AAA,BBB)

## select null hypotheses & combine (model-based like standard errors)
AA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"),
             robust = FALSE)
BB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"),
             robust = FALSE)
ZZ <- rbind(AA,BB)
}
\keyword{methods}
