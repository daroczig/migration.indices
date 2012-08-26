\name{migration.gini.exchange}
\alias{migration.gini.exchange}
\title{Exchange Gini Index}
\usage{
  migration.gini.exchange(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  Exchange Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.exchange(migration.hyp)     # 0.05555556
migration.gini.exchange(migration.hyp2)    # 0.04166667
}
}
\references{
  David A. Plane and Gordon F. Mulligan (1997) Measuring
  Spatial Focusing in a Migration System. \emph{Demography}
  \bold{34}, 251--262
}

