% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tva.R
\name{summary,stantvafit-method}
\alias{summary,stantvafit-method}
\title{Summary method for RStanTVA fits}
\usage{
\S4method{summary}{stantvafit}(object, pars, ...)
}
\arguments{
\item{object}{The RStanTVA fit.}

\item{pars}{(Optional) A character vector of variable names to extract.}

\item{...}{Additional arguments passed to \code{\link[rstan:summary,stanfit-method]{rstan::summary()}}, e.g. \code{probs} and \code{use_cache}.}
}
\value{
See \code{\link[rstan:summary,stanfit-method]{rstan::summary()}} for details.
}
\description{
Summarize the distributions of estimated parameters and derived quantities using the posterior draws.
}
\examples{
\donttest{
f <- read_stantva_fit("fit.rds")
summary(f, "C_Intercept", probs = c(.025, .975))
}
}
