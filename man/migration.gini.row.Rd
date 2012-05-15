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
  number
}
\description{
  Rows Gini Index
}
\examples{
\dontrun{
m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0), 3, 3)
migration.gini.row(m)  # 0
m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
migration.gini.row(m2) # 0.02083333
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

