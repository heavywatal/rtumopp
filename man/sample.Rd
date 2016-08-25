% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\name{make_samples}
\alias{make_samples}
\title{Make ms-like matrix from genealogy column of population}
\usage{
make_samples(genealogy, nsam = length(genealogy), mu = NULL,
  segsites = NULL)
}
\arguments{
\item{genealogy}{list of int vectors}

\item{nsam}{number of samples}

\item{mu}{mutation rate per cell division (ignored if segsites is given)}

\item{segsites}{number of segregating sites}
}
\value{
int matrix [nsam x segsites]
}
\description{
Make ms-like matrix from genealogy column of population
}
