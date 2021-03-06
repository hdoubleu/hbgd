% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_holdout.R
\name{add_holdout_ind}
\alias{add_holdout_ind}
\title{Add indicator column for per-subject holdout}
\usage{
add_holdout_ind(dat, random = TRUE)
}
\arguments{
\item{dat}{data}

\item{random}{if TRUE, a random observation per subject will be designated as the holdout, if FALSE, the endpoint for each subject will be designated as the holdout}
}
\description{
Add indicator column for per-subject holdout
}
