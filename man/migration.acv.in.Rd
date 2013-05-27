\name{migration.acv.in}
\alias{migration.acv.in}
\title{Aggregated In-migration Coefficient of Variation}
\usage{
  migration.acv.in(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  Aggregated In-migration Coefficient of Variation
}
\examples{
\dontrun{
data(migration.hyp)
migration.acv.in(migration.hyp)    # 0.3333333
migration.acv.in(migration.hyp2)   # 0.25
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
  \code{\link{migration.acv.out}}
  \code{\link{migration.acv}}
}

