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
  number
}
\description{
  Aggregated Out-migration Coefficient of Variation
}
\examples{
\dontrun{
data(migration.hyp)
migration.acv.out(migration.hyp)    # 0
migration.acv.out(migration.hyp2)   # 0.125
}
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

