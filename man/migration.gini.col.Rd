\name{migration.gini.col}
\alias{migration.gini.col}
\title{Columns Gini Index}
\usage{
  migration.gini.col(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A number between 0 and 1 where 0 means no spatial
  focusing and 1 shows maximum focusing.
}
\description{
  The Columns Gini index concentrates on the "relative
  extent to which the destination selections of
  in-migrations are spatially focused": \deqn{G^T_R =
  \frac{\sum_j \sum_{i \neq j} \sum_{g \neq i,j} | M_{ij} -
  M_{gj} | }{ (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}}
  This implementation solves the above formula by computing
  the \code{dist} matrix for each columns.
}
\examples{
data(migration.hyp)
migration.gini.col(migration.hyp)  # 0.05555556
migration.gini.col(migration.hyp2) # 0.04166667
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 }
}
\seealso{
  \code{\link{migration.gini.row}}
  \code{\link{migration.gini.col.standardized}}
}

