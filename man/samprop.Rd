\name{samprop}
\alias{samprop}
\title{Random Sampling: the Sampford Algorithm}
\description{
  \code{samprop} implements a modified version of the Sampford algorithm
  to obtain a sample without replacement and with unequal probabilities.
  The probability of the sample is proportional to the product of the
  probabilities of units contained in the sample.}

\usage{
 samprop(x, size, prob = NULL)
}
\arguments{
  \item{x}{Either a (numeric, complex, character or logical) vector of
    more than one element from which to choose, or a positive integer.}
  \item{size}{A positive integer giving the number of items to choose.}
  \item{prob}{A vector of weights for obtaining the elements of the
    vector being sampled.}
}
\details{
  If \code{x} has length 1, sampling takes place from \code{1:x}.

  The \code{prob} argument gives a vector of weights for obtaining the
  elements of the vector being sampled.  They need not sum to one, but
  they should be nonnegative.  The number of nonzero weights must be at
  least \code{size}.  Each weight cannot exceed a certain quantity, in
  order to satisfy Sampford algorithm requirements (see
  \code{\link[sampfling]{sampfle}}).

  This procedure draws a sample without replacement, which probability
  is proportional to the product of weights of units contained in the
  sample.
}
\value{
  A vector of size \code{size}, which contains the sampled units from
  \code{x}, or indices if \code{x} is an integer.
  }
\references{
  Cochran, W.G. (1977) \emph{Sampling techniques}, John Wiley and Sons.

  Sampford, M.R. (1967) On sampling without replacement with unequal
  probabilities of selection, \emph{Biometrika} 54:499-513.
}
\author{Carlos Enrique Carleos Artime \email{carleos@vodafone.es.}}
\seealso{ 
  \code{\link{sampfle}}
}
\examples{
library(sampfling)
samprop(10, 3, 1:10)
samprop(c("a","b","c"), 2, c(0.3,0.3,0.4))
}
\keyword{distribution}


