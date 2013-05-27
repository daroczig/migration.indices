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
  number
}
\description{
  In-migration Coefficient of Variation
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

