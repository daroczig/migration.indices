\name{migration.gini.row}
\alias{migration.gini.row}
\title{Rows Gini Index}
\usage{
  migration.gini.row(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A number between 0 and 1 where 0 means no spatial
  focusing and 1 shows maximum focusing.
}
\description{
  The Rows Gini index concentrates on the "relative extent
  to which the destination selections of out-migrations are
  spatially focused": \deqn{G^T_R = \frac{\sum_i \sum_{j
  \neq i} \sum_{h \neq i,j} | M_{ij} - M_{ih} | }{
  (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}} This
  implementation solves the above formula by computing the
  \code{dist} matrix for each row.
}
\examples{
data(migration.hyp)
migration.gini.row(migration.hyp)  # 0
migration.gini.row(migration.hyp2) # 0.02083333
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 }
}
\seealso{
  \code{\link{migration.gini.col}}
  \code{\link{migration.gini.row.standardized}}
}

