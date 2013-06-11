\name{migration.gini.exchange}
\alias{migration.gini.exchange}
\title{Exchange Gini Index}
\usage{
  migration.gini.exchange(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A number between 0 and 1 where 0 means no spatial
  focusing and 1 shows maximum focusing.
}
\description{
  The Exchange Gini Index "indicates the contribution to
  spatial focusing represented by the \eqn{n(n-q)} net
  interchanges in the system": \deqn{G^T_{RC, CR} =
  \frac{\sum_i \sum_{j \neq i} | M_{ij} - M_{ji} | }{
  (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}} This
  implementation solves the above formula by simply
  substracting the transposed matrix's values from the
  original one at one go.
}
\examples{
data(migration.hyp)
migration.gini.exchange(migration.hyp)     # 0.05555556
migration.gini.exchange(migration.hyp2)    # 0.04166667
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 }
}
\seealso{
  \code{\link{migration.gini}}
  \code{\link{migration.gini.exchange.standardized}}
}

