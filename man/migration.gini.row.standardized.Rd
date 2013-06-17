\name{migration.gini.row.standardized}
\alias{migration.gini.row.standardized}
\title{Standardized Rows Gini Index}
\usage{
  migration.gini.row.standardized(m,
    gini.total = migration.gini.total(m, FALSE))
}
\arguments{
  \item{m}{migration matrix}

  \item{gini.total}{optionally pass the pre-computed Total
  Flows Gini Index to save computational resources}
}
\value{
  A percentage range from 0\% to 100\% where 0\% means that
  the migration flows are uniform, while a higher value
  indicates spatial focusing.
}
\description{
  The standardized version of the Rows Gini Index
  (\code{\link{migration.gini.row}}) by dividing that with
  the Total Flows Gini Index
  (\code{\link{migration.gini.total}}): \deqn{G^{T*}_R =
  100\frac{G^T_R}{G^T}} As this index is standardized, it
  "facilitate comparisons from one period to the next of
  the rows" indices.
}
\examples{
data(migration.hyp)
migration.gini.row.standardized(migration.hyp)     # 0
migration.gini.row.standardized(migration.hyp2)    # 11.11111
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

