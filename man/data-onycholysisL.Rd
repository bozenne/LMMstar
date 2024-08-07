% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc-data.R
\docType{data}
\name{onycholysisL}
\alias{onycholysisL}
\title{Data From The toenail onycholysis Study (Long Format)}
\usage{
data(onycholysisL)
}
\description{
Data from the toenail onychomycosis study,
a randomized double-blind comparative study comparing Terbinafine (250 mg/day) against Itraconazole (200 mg/day) over 7 weeks.
This dataset is in the long format (i.e. one line per measurement).

\itemize{
\item \code{id}: patient identifier.
\item \code{group}: treatment arm to which the patient has been randomized (\code{"itraconazole"} or \code{"terbinafine"}).
\item \code{treatment}: treatment recieved by the patient at a given timepoint (\code{"none"}, \code{"itraconazole"}, or \code{"terbinafine"}).
\item \code{visit}: index of time at which the measurement was taken (\code{1} to \code{7}).
\item \code{time}: scheduled time, in months, for each measurement (\code{0}, \code{4}, \code{8}, \code{12}, \code{24}, \code{36}, \code{48}).
\item \code{obstime}: time, in months, at which the measurement was taken (numeric between \code{0} and \code{18.5}).
\item \code{response}: degree of onycholysis (separation of the nail plate from the nail-bed). Can be 0 for none or mild, 1 for moderate or severe.
}
}
\references{
De Backer et al. Twelve weeks of continuous oral therapy for toenail onychomycosis caused by dermatophytes: A double-blind comparative trial of terbinafine 250 mg/day versus itraconazole 200 mg/day,  Journal of the American Academy of Dermatology, (2020) 110. doi: 10.1016/s0190-9622(98)70486-4
}
\keyword{datasets}
