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
  number
}
\description{
  Columns Gini Index
}
\examples{
\dontrun{
m <- matrix(c(0, 20, 30, 10, 0, 30, 10, 20, 0),3,3)
migration.gini.col(m)  # 0.05555556
m2 <- matrix(c(0, 30, 20, 20, 0, 20, 20, 50, 0), 3, 3 )
migration.gini.col(m2) # 0.04166667
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

