% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc-data.R
\docType{data}
\name{vitaminW}
\alias{vitaminW}
\title{Data From The Vitamin Study (Wide Format)}
\usage{
data(vitaminW)
}
\description{
Data from the vitamin Study,
a randomized study where the growth of guinea pigs was monitored before and after intake of vitamin E/placebo.
The weight of each guinea pig was recorded at the end of week 1, 3, 4, 5, 6, and 7. Vitamin E/placebo is given at the beginning of week 5.
This dataset is in the wide format (i.e. one line per patient).

\itemize{
\item group Treatment group: vitamin or placebo.
\item animal Identifier
\item weigth1 weight (in g) of the pig at the end of week 1 (before treatment).
\item weigth3 weight (in g) of the pig at the end of week 3 (before treatment).
\item weigth4 weight (in g) of the pig at the end of week 4 (before treatment).
\item weigth5 weight (in g) of the pig at the end of week 5 (after treatment).
\item weigth6 weight (in g) of the pig at the end of week 6 (after treatment).
\item weigth7 weight (in g) of the pig at the end of week 7 (after treatment).
}
}
\examples{
\dontrun{
vitaminW <- read.table("inst/dataTXT/vitamin.txt", header = TRUE, na.string = ".")
vitaminW$group <- as.factor(vitaminW$group)
vitaminW$animal <- as.factor(vitaminW$animal)
save(vitaminW, file = "data/vitaminW.rda")
}
}
\references{
TODO
}
\keyword{data}
