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
  number
}
\description{
  Out-migration Coefficient of Variation
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

