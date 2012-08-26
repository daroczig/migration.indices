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
data(migration.hyp)
migration.gini.col(migration.hyp)  # 0.05555556
migration.gini.col(migration.hyp2) # 0.04166667
}
}
\references{
  David A. Plane and Gordon F. Mulligan (1997) Measuring
  Spatial Focusing in a Migration System. \emph{Demography}
  \bold{34}, 251--262
}

