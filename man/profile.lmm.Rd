% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile.R
\name{profile.lmm}
\alias{profile.lmm}
\title{Log-Likelihood Contour For a Linear Mixed Model}
\usage{
\method{profile}{lmm}(
  fitted,
  effects = NULL,
  profile.likelihood = FALSE,
  maxpts = NULL,
  level = 0.95,
  df = NULL,
  trace = FALSE,
  transform.sigma = NULL,
  transform.k = NULL,
  transform.rho = NULL,
  transform.names = TRUE,
  ...
)
}
\arguments{
\item{fitted}{a \code{lmm} object.}

\item{effects}{[character vector] name of the parameters to be constrained.
Alternatively can be the type of parameters, e.g. \code{"mean"}, \code{"variance"}, \code{"correlation"}, or \code{"all"}.}

\item{profile.likelihood}{[FALSE,TRUE,"ci"] Should the unconstrained parameter(s) be kept at their (MLE) value (\code{FALSE}),
or set to the value maximizing the constrained likelihood (\code{TRUE}),
or should confidence intervals be computed for the parameters defined in argument \code{effects}.}

\item{maxpts}{[integer, >0] number of points use to discretize the likelihood, \code{maxpts} points smaller than the MLE and \code{maxpts} points higher than the MLE.}

\item{level}{[numeric, 0-1] the confidence level of the confidence intervals.
Used to decide about the range of values for each parameter when argument \code{profile.likelihood} is \code{TRUE} or \code{FALSE}.}

\item{df}{[logical] Should a Student's t-distribution be used to model the distribution of the coefficients when evaluating the confidence intervals. Otherwise a normal distribution is used.
Ignored when \code{profile.likelihood} is \code{"ci"}.}

\item{trace}{[logical] Show the progress of the execution of the function.}

\item{transform.sigma}{[character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.}

\item{transform.k}{[character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.}

\item{transform.rho}{[character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.}

\item{transform.names}{[logical] Should the name of the coefficients be updated to reflect the transformation that has been used?}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
[profile.likelihood = TRUE/FALSE] A data.frame object containing the log-likelihood for various parameter values. \cr
[profile.likelihood = "ci"] A data.frame object containing the REML or ML estimated parameter (\code{"estimate"}),
lower and upper bound of the profile-likelihood confidence interval (\code{"lower"}, \code{"upper"}),
and the discrepancy in log-likelihood between the bounds found by the root finding algorithm and the requested difference in log-likelihood reflecting the confidence level (\code{"error.lower"}, \code{"error.upper"}).
}
\description{
Evaluate the (restricted) log-likelihood around Maximum Likelihood Estimate (MLE) of a linear mixed model.
The values of a given parameter are varied over a pre-defined grid and the corresponding (contrained) likelihood w.r.t. each value is evaluated.
The other parameters are either kept constant or set to maximize the contrained likelihood (profile likelihood).
In the latter case, confidence intervals consistent with a likelihood ratio test (LRT) can be output.
}
\details{
Each parameter defined by the argument \code{effets} is treated separately either to evaluate the constrained likelihood or compute a confidence interval (argument \code{profile.likelihood}).

Confidence intervals are evaluating such that the lower and upper bound correspond to the same likelihood value,
aiming at having intervals where lack of coverage is equaly likely due to low or to high bounds.
This is performed using a root finding algorithm (\code{\link{uniroot}}).

The constrained likelihood is evaluated as follow:\itemize{
\item the confidence interval of a parameter is discretized with \code{maxpt} points. Increasing the confidence level will lead to a larger range of parameter values.
\item this parameter is set to each discretization value.
\item the other parameters are either set to the (unconstrained) MLE (\code{profile.likelihood=FALSE})
or to constrained MLE (\code{profile.likelihood=TRUE}). The latter case is much more computer intensive as it implies re-running the estimation procedure.
\item the (restricted) log-likelihood is evaluated.
}
Since a locally quadratic log-likelihood with an Hessian equivariant in law implies normally distributed estimates (Geyer 2013)
it can help trusting confidence intervals and p-values in small samples with a non-normally distributed outcome.
}
\examples{
### Linear regression ####
data(gastricbypassW, package = "LMMstar")
e.lmm <- lmm(weight2 ~ weight1 + glucagonAUC1, data = gastricbypassW)

#### likelihood along a parameter axis (slice)
## no transformation
e.sliceNone <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "none")
plot(e.sliceNone)
## transformation
e.sliceLog <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "log")
plot(e.sliceLog)

#### profile likelihood (local maxima of the likelihood - crest line)
\dontrun{
e.pro <- profile(e.lmm, effects = "all", profile.likelihood = TRUE)
plot(e.pro)
}

#### confidence interval based on profile likelihood
\dontrun{
e.PLCI <- profile(e.lmm, effects = c("weight1","sigma"), profile.likelihood = "ci")
e.PLCI
}

### Random intercept model ####
## Data shown in Sahai and Ageel (2000) page 122.
## The analysis of variance: fixed, random, and mixed models. Springer
df <- rbind(data.frame(Y = c(7.2, 7.7, 8, 8.1, 8.3, 8.4, 8.4, 8.5, 8.6, 8.7,
                             9.1, 9.1, 9.1, 9.8, 10.1, 10.3), type = "Hb SS"),
           data.frame(Y = c(8.1, 9.2, 10, 10.4, 10.6, 10.9, 11.1, 11.9, 12, 12.1),
                      type = "Hb S/beta"),
           data.frame(Y = c(10.7, 11.3, 11.5, 11.6, 11.7, 11.8, 12, 12.1, 12.3,
                            12.6, 12.6, 13.3, 13.3, 13.8, 13.9), type = "HB SC")
)
df$type <- factor(df$type, levels = unique(df$type))

## retrive first column of table I in the original publication
## (doi:  10.1136/bmj.282.6260.283)
round(tapply(df$Y,df$type,mean),1)
round(tapply(df$Y,df$type,sd),1)

e.RI <- lmm(Y ~ (1|type), data = df)

#### delta method
confint(e.RI, effects = "correlation", df = FALSE) ## 0.76 [0.111; 0.955]
 
ranef(e.RI, effects = "variance", se = TRUE) ## 3.17 [0.429; 23.423]
ranef(e.RI, effects = "variance", se = TRUE, transform=FALSE) ## 3.17 [-3.167; 9.510]
## same as confint(e.RI, effects = "correlation", transform.rho = "cov", df = FALSE)

#### profile likelihood
\dontrun{
plot(profile(e.RI, effects = "correlation", profile.likelihood = TRUE, df = FALSE))
plot(profile(e.RI, effects = "correlation", profile.likelihood = TRUE, df = FALSE,
             transform.rho = "none", maxpts = seq(0.4,0.99,by=0.01)))

profile(e.RI, effects = "correlation", profile.likelihood = "ci")
## 0.760 [0.378; 0.983]
profile(e.RI, effects = "correlation", profile.likelihood = "ci", transform.rho = "cov")
}
}
\references{
Geyer, C. J. (2013). Asymptotics of maximum likelihood without the lln or clt or sample size going to infinity. In Advances in Modern Statistical Theory and Applications: A Festschrift in honor of Morris L. Eaton, pages 1–24. Institute of Mathematical Statistics.
}
\keyword{htest}
