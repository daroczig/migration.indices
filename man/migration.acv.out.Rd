\name{migration.acv.out}
\alias{migration.acv.out}
\title{Aggregated Out-migration Coefficient of Variation}
\usage{
  migration.acv.out(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A number where a higher (\eqn{\neq 0}) shows more spatial
  focus.
}
\description{
  The Aggregated Out-migration Coefficient of Variation is
  the weighted average of the Out-migration Coefficient of
  Variation (\code{\link{migration.cv.out}}).
}
\examples{
data(migration.hyp)
migration.acv.out(migration.hyp)    # 0
migration.acv.out(migration.hyp2)   # 0.125
}
\references{
  \itemize{ \item Andrei Rogers and Stuart Sweeney (1998)
  Measuring the Spatial Focus of Migration Patterns.
  \emph{The Professional Geographer} \bold{50}, 232--242 }
}
\seealso{
  \code{\link{migration.cv.in}}
  \code{\link{migration.cv.out}}
  \code{\link{migration.acv.in}}
  \code{\link{migration.acv}}
}

