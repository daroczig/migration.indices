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
data(migration.hyp)
migration.gini.row(migration.hyp)  # 0
migration.gini.row(migration.hyp2) # 0.02083333
}
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 }
}
\seealso{
  \code{\link{migration.gini.col}}
  \code{\link{migration.gini.row.standardized}}
}

