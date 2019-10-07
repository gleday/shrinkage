\name{eightnum}
\alias{eightnum}

\title{
  eight-number summary
}

\description{
  This function computes an eight-number summary for a sample
}

\usage{
  eightnum(x)	
}

\arguments{
  \item{x}{vector of samples}
}

\examples{
\dontrun{
  set.seed(12345)
  xx <- rnorm(10000)
  eightnum(xx)
}
}
