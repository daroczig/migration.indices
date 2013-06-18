\name{migration.weighted.gini.out}
\alias{migration.weighted.gini.out}
\title{Migration-weighted Out-migration Gini Index}
\usage{
  migration.weighted.gini.out(m,
    mgo = migration.gini.out(m))
}
\arguments{
  \item{m}{migration matrix}

  \item{mgo}{optionally passed (precomputed) Migration
  In-migration Gini Index}
}
\description{
  The Migration-weighted Out-migration Gini Index is a
  weighted version of the Out-migration Field Gini Index
  (\code{\link{migration.gini.out}}) "according to the zone
  of destination's share of total migration and the mean of
  the weighted values is computed as": \deqn{MWG^O = \frac{
  \sum_i G^O_i \frac{\sum_j M_{ij}}{\sum_{ij} M_{ij}}}{n}}
}
\examples{
data(migration.hyp)
migration.weighted.gini.out(migration.hyp)   # 0
migration.weighted.gini.out(migration.hyp2)  # 0.02083333
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

