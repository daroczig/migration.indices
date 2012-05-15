\name{migration.gini.total}
\alias{migration.gini.total}
\title{Total Flows Gini Index}
\usage{
  migration.gini.total(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  Total Flows Gini Index
}
\examples{
\dontrun{
m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0), 3, 3)
migration.gini.total(m)        # 0.2222222
m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
migration.gini.total(m2)       # 0.1875
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

