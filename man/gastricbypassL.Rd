% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc-data.R
\docType{data}
\name{gastricbypassL}
\alias{gastricbypassL}
\title{Data From The Gastric Bypass Study (Long Format)}
\usage{
data(gastricbypassL)
}
\description{
Data from the gastric bypass study
where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
This dataset is in the long format (i.e. one line per measurement).

\itemize{
\item id Patient identifier
\item visit The visit index.
\item time The time at which the visit took place.
\item weight Bodyweight (in kg) measured during the visit.
\item glucagon Glucagon measured during the visit.
}
}
\examples{
\dontrun{
data("gastricbypassW")
dtW <- data.table::as.data.table(gastricbypassW)
dtL <- data.table::melt(dtW, id.vars = "id",
                        measure.vars = patterns("weight","glucagonAUC"),
                        value.name = c("weight","glucagon"), variable.name = "time")
gastricbypassL <- as.data.frame(dtL)
gastricbypassL$visit <- gastricbypassL$time
gastricbypassL$time <- factor(gastricbypassL$visit, levels = 1:4,
                              labels = c("3 months before surgery","1 week before surgery",
                                         "1 week after surgery","3 months after surgery"))
gastricbypassL <- gastricbypassL[,c("id","visit","time","weight","glucagon")]
save(gastricbypassL, file = "data/gastricbypassL.rda")

str(gastricbypassL)
}
}
\references{
The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
}
\keyword{data}
