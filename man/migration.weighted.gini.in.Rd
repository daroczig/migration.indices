\name{migration.weighted.gini.in}
\alias{migration.weighted.gini.in}
\title{Migration-weighted In-migration Gini Index}
\usage{
  migration.weighted.gini.in(m, mgi = migration.gini.in(m))
}
\arguments{
  \item{m}{migration matrix}

  \item{mgi}{optionally passed (precomputed) Migration
  In-migration Gini Index}
}
\description{
  The Migration-weighted In-migration Gini Index is a
  weighted version of the In-migration Field Gini Index
  (\code{\link{migration.gini.in}}) "according to the zone
  of destination's share of total migration and the mean of
  the weighted values is computed as": \deqn{MWG^I = \frac{
  \sum_j G^I_j \frac{\sum_j M_{ij}}{\sum_{ij} M_{ij}}}{n}}
}
\examples{
data(migration.hyp)
migration.weighted.gini.in(migration.hyp)   # 0.1222222
migration.weighted.gini.in(migration.hyp2)  # 0.05238095
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
  \code{\link{migration.gini}}
  \code{\link{migration.gini.in}}
  \code{\link{migration.weighted.gini.out}}
  \code{\link{migration.weighted.gini.mean}}
}

