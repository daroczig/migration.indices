\name{migration.gini.out}
\alias{migration.gini.out}
\title{Out-migration Field Gini Index}
\usage{
  migration.gini.out(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A numeric vector with the range of 0 to 1 where 0 means
  no spatial focusing and 1 shows maximum focusing.
}
\description{
  The Out-migration Field Gini Index is a decomposed
  version of the Rows Gini Index
  (\code{\link{migration.gini.row}}) representing "the
  contribution of each region's row to the total index" ()
  (\code{\link{migration.gini.total}}): \deqn{G^O_i =
  \frac{\sum_{j \neq i} \sum_{l \neq i,j} | M_{ij} - M_{il}
  | }{ 2(n-2) \sum_{j \neq k} M_{ij}}} These Gini indices
  facilitates the direct comparison of different
  territories without further standardization.
}
\examples{
data(migration.hyp)
migration.gini.out(migration.hyp)  # 0 0 0
migration.gini.out(migration.hyp2) # 0.000 0.25 0.000
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 \item M. Bell, M.
  Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell
  and G. Hugo (2002) Cross-National Comparison of Internal
  Migration. Issues and Measures. \emph{Journal of the
  Royal Statistical Society. Series A (Statistics in
  Society)} \bold{165}, 435--464 }
}
\seealso{
  \code{\link{migration.gini}}
  \code{\link{migration.gini.in}}
  \code{\link{migration.weighted.gini.out}}
}

