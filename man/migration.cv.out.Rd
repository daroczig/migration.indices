\name{migration.cv.out}
\alias{migration.cv.out}
\title{Out-migration Coefficient of Variation}
\usage{
  migration.cv.out(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A numeric vector of standardized values where a higher
  (\eqn{\neq 0}) shows more spatial focus.
}
\description{
  As "the coefficient of variation is defined as the
  standard deviation to mean ratio of a distribution", the
  Out-migration Coefficient of Variation is computed by
  dividing the standard deviation (with the nominator being
  \eqn{n} instead of \eqn{n-1}) of the out-migration flows
  by the mean.
}
\examples{
\dontrun{
data(migration.hyp)
migration.cv.out(migration.hyp)    # 0 0 0
migration.cv.out(migration.hyp2)   # 0.00 0.25 0.00
}
}
\references{
  \itemize{ \item Andrei Rogers and Stuart Sweeney (1998)
  Measuring the Spatial Focus of Migration Patterns.
  \emph{The Professional Geographer} \bold{50}, 232--242 }
}
\seealso{
  \code{\link{migration.cv.in}}
  \code{\link{migration.acv.in}}
  \code{\link{migration.acv.out}}
  \code{\link{migration.acv}}
}

