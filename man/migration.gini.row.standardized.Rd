\name{migration.gini.row.standardized}
\alias{migration.gini.row.standardized}
\title{Standardized Rows Gini Index}
\usage{
  migration.gini.row.standardized(m,
    migration.gini.total = migration.gini.total(m))
}
\arguments{
  \item{m}{migration matrix}

  \item{migration.gini.total}{optionally pass the
  pre-computed Total Flows Gini Index to save resources}
}
\value{
  number
}
\description{
  Standardized Rows Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.row.standardized(migration.hyp)     # 0
migration.gini.row.standardized(migration.hyp2)    # 11.11111
}
}
\author{
  Gergely Daróczi
}
\references{
  David A. Plane and Gordon F. Mulligan (1997): Measuring
  Spatial Focusing in a Migration System. In. Demography,
  Vol. 34, No. 2 (May, 1997), pp. 251-262
}

