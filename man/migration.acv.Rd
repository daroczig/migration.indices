\name{migration.acv}
\alias{migration.acv}
\title{Aggregated System-wide Coefficient of Variation}
\usage{
  migration.acv(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A number where a higher (\eqn{\neq 0}) shows more spatial
  focus.
}
\description{
  The Aggregated System-wide Coefficient of Variation is
  simply the sum of the Aggregated In-migration
  (\code{\link{migration.acv.in}}) and the Aggregated
  Out-migration Coefficient of Variation
  (\code{\link{migration.acv.out}}).
}
\examples{
data(migration.hyp)
migration.acv(migration.hyp)    # 0.3333333
migration.acv(migration.hyp2)   # 0.375
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
  \code{\link{migration.acv.out}}
}

