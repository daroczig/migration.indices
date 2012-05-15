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
  number
}
\description{
  Standardized Columns Gini Index
}
\examples{
\dontrun{
m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
migration.gini.col.standardized(m)     # 25
m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
migration.gini.col.standardized(m2)    # 22.22222
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

