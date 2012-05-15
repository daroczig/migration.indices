\name{migration.gini.in}
\alias{migration.gini.in}
\title{In-migration Field Gini Index}
\usage{
  migration.gini.in(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  In-migration Field Gini Index
}
\examples{
\dontrun{
m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
migration.gini.in(m)   # 0.1000000 0.2500000 0.1666667
m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
migration.gini.in(m2)  # 0.1000000 0.0000000 0.2142857
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

