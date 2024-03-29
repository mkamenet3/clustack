% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustack.R
\name{detectclusters}
\alias{detectclusters}
\title{detectclusters}
\usage{
detectclusters(
  sparsemat,
  Ex,
  Yx,
  numCenters,
  Time,
  maxclust,
  byloc = FALSE,
  model = "poisson",
  overdisp.est
)
}
\arguments{
\item{sparsemat}{Large sparse matrix of potential clusters.}

\item{Ex}{Vector of expected counts.}

\item{Yx}{Vector of observed counts.}

\item{numCenters}{Number of centroids.}

\item{Time}{Number of time periods.}

\item{maxclust}{Maximum number of clusters allowed. TODO - allow this to be unknown.}

\item{byloc}{If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{FALSE} (detection by potential cluster).}

\item{model}{A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.}

\item{overdisp.est}{Overdispersion estimate.}
}
\value{
Returns a large list.
}
\description{
Detect disease clusters either by location or by potential cluster.
}
\details{
The elements of the returned list: 1) the weighted relative risks for each cell; 2) the LRT statistic; 3) the selected number of clusters by (Q)BIC; 4) the selected number of clusters by (Q)AIC; 5) the selected number of clusters by (Q)AICc; 6) the selected number of clusters by (Q)BIC when forced to detect a cluster; 7) the selected number of clusters by (Q)AIC when forced to detect a cluster; 8) the selected number of clusters by (Q)AICc when forced to detect a cluster; 9) matrix of weights; 10) maximum potential cluster or location IDs identified; 11) large matrix of single-cluster relative risk estimates.
}
