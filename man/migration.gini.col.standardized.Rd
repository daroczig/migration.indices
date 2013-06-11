\name{migration.gini.col.standardized}
\alias{migration.gini.col.standardized}
\title{Standardized Columns Gini Index}
\usage{
  migration.gini.col.standardized(m,
    gini.total = migration.gini.total(m))
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
  The standardized version of the Columns Gini Index
  (\code{\link{migration.gini.col}}) by dividing that with
  the Total Flows Gini Index
  (\code{\link{migration.gini.total}}): \deqn{G^{T*}_C =
  100\frac{G^T_C}{G^T}} As this index is standardized, it
  "facilitate comparisons from one period to the next" of
  the columns indices.
}
\examples{
data(migration.hyp)
migration.gini.col.standardized(migration.hyp)     # 25
migration.gini.col.standardized(migration.hyp2)    # 22.22222
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

