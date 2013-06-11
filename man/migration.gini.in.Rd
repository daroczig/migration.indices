\name{migration.gini.in}
\alias{migration.gini.in}
\title{In-migration Field Gini Index}
\usage{
  migration.gini.in(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A numeric vector with the range of 0 to 1 where 0 means
  no spatial focusing and 1 shows maximum focusing.
}
\description{
  The In-migration Field Gini Index is a decomposed version
  of the Columns Gini Index
  (\code{\link{migration.gini.col}}) representing "the
  contribution of each region's columns to the total index"
  () (\code{\link{migration.gini.total}}): \deqn{G^I_j =
  \frac{\sum_{i \neq j} \sum_{k \neq j,i} | M_{ij} - M_{kj}
  | }{ 2(n-2) \sum_{i \neq j} M_{ij}}} These Gini indices
  facilitates the direct comparison of different
  territories without further standardization.
}
\examples{
data(migration.hyp)
migration.gini.in(migration.hyp)   # 0.2000000 0.5000000 0.3333333
migration.gini.in(migration.hyp2)  # 0.2000000 0.0000000 0.4285714
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
  \code{\link{migration.gini.out}}
  \code{\link{migration.weighted.gini.in}}
}

