% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{predict.avgSuperLearner}
\alias{predict.avgSuperLearner}
\title{Method to average predictions over multiple SuperLearners}
\usage{
\method{predict}{avgSuperLearner}(x, newdata, ...)
}
\arguments{
\item{x}{Object of class \code{avgSuperLearner}}

\item{newdata}{Prediction \code{data.frame}}

\item{...}{Other arguments (not used)}
}
\value{
Vector of predictions on \code{newdata}
}
\description{
Method to average predictions over multiple SuperLearners
}
