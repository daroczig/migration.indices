\name{migration.gini.total}
\alias{migration.gini.total}
\title{Total Flows Gini Index}
\usage{
  migration.gini.total(m, corrected = TRUE)
}
\arguments{
  \item{m}{migration matrix}

  \item{corrected}{Bell et al. (2002) updated the formula
  of Plane and Mulligan (1997) to have \eqn{2{n(n-1)-1}}
  instead of \eqn{2n(n-1)} in the denominator to "ensure
  that the index can assume the upper limit of 1".}
}
\value{
  A number between 0 and 1 where 0 means no spatial
  focusing and 1 shows that all migrants are found in one
  single flow.
}
\description{
  The Total Gini Index shows the overall concentration of
  migration with a simple number computed by comparing each
  cell of the migration matrix with every other cell except
  for the diagonal: \deqn{G^T = \frac{\sum_i \sum_{j \neq
  i} \sum_k \sum_{l \neq k} | M_{ij} - M_{kl} | }{
  (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}} This
  implementation solves the above formula by a simple loop
  for performance issues to compare all values to the
  others at one go, although smaller migration matrices
  could also be addressed by a much faster \code{dist}
  method. Please see the sources for more details.
}
\examples{
data(migration.hyp)
migration.gini.total(migration.hyp)           # 0.2666667
migration.gini.total(migration.hyp2)          # 0.225
migration.gini.total(migration.hyp, FALSE)    # 0.2222222
migration.gini.total(migration.hyp2, FALSE)   # 0.1875
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
  \code{\link{migration.gini.col}}
  \code{\link{migration.gini.row}}
  \code{\link{migration.gini.exchange}}
  \code{\link{migration.gini.in}}
  \code{\link{migration.gini.out}}
}

