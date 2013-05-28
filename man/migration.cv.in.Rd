\name{migration.cv.in}
\alias{migration.cv.in}
\title{In-migration Coefficient of Variation}
\usage{
  migration.cv.in(m)
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
  In-migration Coefficient of Variation is computed by
  dividing the standard deviation (with the nominator being
  \eqn{n} instead of \eqn{n-1}) of the in-migration flows
  by the mean.
}
\examples{
\dontrun{
data(migration.hyp)
migration.cv.in(migration.hyp)    # 0.2000000 0.5000000 0.3333333
migration.cv.in(migration.hyp2)   # 0.2000000 0.0000000 0.4285714
}
}
\references{
  \itemize{ \item Andrei Rogers and Stuart Sweeney (1998)
  Measuring the Spatial Focus of Migration Patterns.
  \emph{The Professional Geographer} \bold{50}, 232--242 }
}
\seealso{
  \code{\link{migration.cv.out}}
  \code{\link{migration.acv.in}}
  \code{\link{migration.acv.out}}
  \code{\link{migration.acv}}
}

