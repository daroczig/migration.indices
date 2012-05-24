\name{migration.gini.exchange.standardized}
\alias{migration.gini.exchange.standardized}
\title{Standardized Exchange Gini Index}
\usage{
  migration.gini.exchange.standardized(m,
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
  Standardized Exchange Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.exchange.standardized(migration.hyp)  # 25
migration.gini.exchange.standardized(migration.hyp2) # 22.22222
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

