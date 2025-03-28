% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tva.R
\name{write_stantva_model}
\alias{write_stantva_model}
\title{Write StanTVA model}
\usage{
write_stantva_model(model, file = stdout())
}
\arguments{
\item{model}{The StanTVA model object.}

\item{file}{The file name.}
}
\description{
Writes a StanTVA model to a file.
}
\examples{
\dontrun{write_stantva_model(model, "model.stan")}
}
