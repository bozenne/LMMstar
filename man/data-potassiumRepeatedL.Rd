% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc-data.R
\docType{data}
\name{potassiumRepeatedL}
\alias{potassiumRepeatedL}
\title{Data From The Potassium Intake Study (Long Format with intermediate measurements)}
\usage{
data(potassiumRepeatedL)
}
\description{
Data from the potassium intake study,
a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
This dataset is in the long format (i.e. one line per measurement) and contains measurement over 6 timepoints for each time period.

\itemize{
\item \code{id}: patient identifier.
\item \code{sequence}: treatment group to which the patient has been randomized.
\item \code{period}: time period.
\item \code{treatment}: treatment during the time period.
\item \code{time}: time within each period.
\item \code{aldo}: ??
}
}
\references{
Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
Nephrol Dial Transplant (2020) 110. doi: 10.1093/ndt/gfaa114
}
\keyword{datasets}
