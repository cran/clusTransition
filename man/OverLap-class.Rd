% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R
\docType{class}
\name{OverLap-class}
\alias{OverLap-class}
\alias{OverLap}
\title{Overlap between clusters}
\description{
Contains matrix of similarity indices between clusters, after clustering dynamic datasets at consecutive time points.
}
\section{Slots}{

\describe{
\item{\code{Overlap}}{A numeric matrix containing the similarity index between clusters extracted at time point \code{t_1} and \code{t_2}.
The rows of the matrix illustrate clusters extracted from first clustering \eqn{\xi_1(time point t_1)},whereas columns represent
clusters extracted from second clustering \eqn{\xi_2(time point t_2)}.}

\item{\code{rx}}{A numeric vector containg radius of each cluster from first clustering \eqn{\xi_1}.}

\item{\code{ry}}{A numeric vector containg radius of each cluster from second clustering \eqn{\xi_2}.}

\item{\code{Centersx}}{A numeric vector containing centers of clusters from first clustering \eqn{\xi_1}.}

\item{\code{Centersy}}{A numeric vector containing centers of clusters from second clustering \eqn{\xi_2}.}

\item{\code{avgDisx}}{A numeric vector containing average distance between points in a cluster from its center in first clustering \eqn{\xi_1}.}

\item{\code{avgDisy}}{A numeric vector containing average distance between points in a cluster from its center in second clustering \eqn{\xi_2}.}

\item{\code{clusterMem}}{A vector of integers containing cluster membership from second clustering \eqn{\xi_2}.}
}}

