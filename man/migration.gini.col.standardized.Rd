\name{migration.gini.col.standardized}
\alias{migration.gini.col.standardized}
\title{Standardized Columns Gini Index}
\usage{
  migration.gini.col.standardized(m,
    migration.gini.total = migration.gini.total(m))
}
\arguments{
  \item{m}{migration matrix}

  \item{migration.gini.total}{optionally pass the
  pre-computed Total Flows Gini Index to save resources}
}
\value{
  percentage
}
\description{
  Standardized Columns Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.col.standardized(migration.hyp)     # 25
migration.gini.col.standardized(migration.hyp2)    # 22.22222
}
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

