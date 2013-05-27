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
migration.cv.out(migration.hyp)    # 0.2000000 0.5000000 0.3333333
migration.cv.out(migration.hyp2)   # 0.2000000 0.0000000 0.4285714
}
}
\references{
  \itemize{ \item Andrei Rogers and Stuart Sweeney (1998)
  Measuring the Spatial Focus of Migration Patterns.
  \emph{The Professional Geographer} \bold{50}, 232--242 }
}

