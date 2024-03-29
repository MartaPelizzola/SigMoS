% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corSig.R
\name{corSig}
\alias{corSig}
\title{Correlation of mutational signatures}
\usage{
corSig(H1, H2)
}
\arguments{
\item{H1}{Numeric matrix of mutational signatures. Each column should represent a signature.}

\item{H2}{Numeric matrix of mutational signatures. Each column should represent a signature.}
}
\value{
Vector of indexes to reorder the second matrix to match the first one.
Also returns the full correlation matrix between the signatures.
}
\description{
Correlation of mutational signatures
}
\examples{
# Find cosine similarity between signatures estimated using SigMoS
# with the Negative Binomial distribution and the Poisson distribution:
res_nb <- sigmos(BRCA21,k=3,patient_specific = TRUE)
res_pois <- sigmos(BRCA21,k=3,method="Poisson")
corSig(res_pois$Signatures, res_nb$Signatures)

}
