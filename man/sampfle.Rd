\name{sampfle}
\alias{sampfle}
\alias{sampfling}
\title{Random Sampling: the Sampford Algorithm}
\description{
  \code{sampfle} implements the Sampford algorithm to obtain a sample
  without replacement and with unequal probabilities.  The probability
  of the sample is known a priori.}

\usage{
 sampfle(x, size, prob = NULL)
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
  they should be nonnegative and, after normalizing, less than
  \code{1/size}.  The number of nonzero weights must be at least \code{size}.

  The Sampford algorithm draws the first unit \code{i} with probability
  proportional to \code{prob}, and subsequent units with probabilities
  proportional to \code{prob/(1-size*prob[i])} with replacement.  As
  soon as a unit is repeated, the entire sample is rejected and the
  process is restarted.  Thus a sample without replacement is obtained.

  The algorithm produces a sample \eqn{s}{s}, subset of \code{x} of size
  \eqn{n}{n}, with probability proportional to
  \deqn{\Pr[s]\propto \left(1-\sum z_i\right)\prod\frac{z_i}{1-n z_i}}
  {(1 - Sum z[i]) Prod(z[i]/(1-n*z[i]))}
  where \eqn{z_i}{z[i]} is the \eqn{i}{i}-th element of \code{prob},
  that is, \code{prob[i]}; and both the sum and the product extend
  over the drawn units, \eqn{i\in s}{i in s}.
}
\value{
  A vector, which contains the sampled units from \code{x}, or indices
  if \code{x} is an integer.
  }
\references{
  Cochran, W.G. (1977) \emph{Sampling techniques}, John Wiley and Sons.

  Sampford, M.R. (1967) On sampling without replacement with unequal
  probabilities of selection, \emph{Biometrika} 54:499-513.
}
\author{Carlos Enrique Carleos Artime \email{carleos@vodafone.es.}}
\seealso{ 
  \code{\link{samprop}}
}
\examples{
library(sampfling)
sampfle(10, 3, 1:10)
sampfle(c("a","b","c"), 2, c(0.3,0.3,0.4))
}
\keyword{distribution}


