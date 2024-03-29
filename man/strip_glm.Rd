% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{strip_glm}
\alias{strip_glm}
\title{Helper function to remove unnecessary output from SuperLearner model to reduce output size}
\usage{
strip_glm(cm, earth = FALSE)
}
\arguments{
\item{cm}{SuperLearner GLM or Earth model}

\item{earth}{flag for if model is SuperLearner Earth model, defaults to GLM}
}
\value{
model with unnecessary features removed
}
\description{
Reference: https://www.r-bloggers.com/2019/12/need-to-save-rs-lm-or-glm-models-trim-the-fat/
}
