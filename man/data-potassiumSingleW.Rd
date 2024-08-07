% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc-data.R
\docType{data}
\name{potassiumSingleW}
\alias{potassiumSingleW}
\title{Data From The Potassium Intake Study (Wide Format)}
\usage{
data(potassiumSingleW)
}
\description{
Data from the potassium intake study,
a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
This dataset is in the wide format (i.e. one line per patient).

\itemize{
\item \code{id}: patient identifier.
\item \code{sequence}: treatment group to which the patient has been randomized.
\item \code{treatment1}: treatment during the first time period.
\item \code{treatment2}: treatment during the second time period.
\item \code{auc1}: area under the curve of ?? during the first time period.
\item \code{auc2}: area under the curve of ?? during the second time period.
\item \code{bsauc1}: ??
\item \code{aldo1}: ??
\item \code{aldo2}: ??
}
}
\references{
Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
Nephrol Dial Transplant (1998) 110. doi: 10.1093/ndt/gfaa114
}
\keyword{datasets}
