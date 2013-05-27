\name{migration.weighted.gini.out}
\alias{migration.weighted.gini.out}
\title{Migration-weighted Out-migration Gini Index}
\usage{
  migration.weighted.gini.out(m, mgo)
}
\arguments{
  \item{m}{migration matrix}

  \item{mgo}{optionally passed (precomputed) Migration
  In-migration Gini Index}
}
\value{
  number
}
\description{
  Migration-weighted Out-migration Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.weighted.gini.out(migration.hyp)   #
migration.weighted.gini.out(migration.hyp2)  #
}
}
\references{
  \itemize{ \item M. Bell, M. Blake, P. Boyle, O.
  Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002)
  Cross-National Comparison of Internal Migration. Issues
  and Measures. \emph{Journal of the Royal Statistical
  Society. Series A (Statistics in Society)} \bold{165},
  435--464 }
}
\seealso{
  \code{\link{migration.weighted.gini.in}}
  \code{\link{migration.weighted.gini.mean}}

  \code{\link{migration.gini}}
  \code{\link{migration.gini.out}}
  \code{\link{migration.weighted.gini.in}}
  \code{\link{migration.weighted.gini.mean}}
}

