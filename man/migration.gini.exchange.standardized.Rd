\name{migration.gini.exchange.standardized}
\alias{migration.gini.exchange.standardized}
\title{Standardized Exchange Gini Index}
\usage{
  migration.gini.exchange.standardized(m,
    gini.total = migration.gini.total(m, FALSE))
}
\arguments{
  \item{m}{migration matrix}

  \item{gini.total}{optionally pass the pre-computed Total
  Flows Gini Index to save resources}
}
\value{
  A percentage range from 0\% to 100\% where 0\% means that
  the migration flows are uniform, while a higher value
  indicates spatial focusing.
}
\description{
  The standardized version of the Exchange Gini Index
  (\code{\link{migration.gini.exchange}}) by dividing that
  with the Total Flows Gini Index
  (\code{\link{migration.gini.total}}): \deqn{G^{T*}_{RC,
  CR} = 100\frac{G^T_{RC, CR}}{G^T}} As this index is
  standardized, it "facilitate comparisons from one period
  to the next" of the exchange indices.
}
\examples{
data(migration.hyp)
migration.gini.exchange.standardized(migration.hyp)  # 25
migration.gini.exchange.standardized(migration.hyp2) # 22.22222
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 }
}
\seealso{
  \code{\link{migration.gini}}
  \code{\link{migration.gini.exchange}}
}

