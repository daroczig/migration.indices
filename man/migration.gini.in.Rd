\name{migration.gini.in}
\alias{migration.gini.in}
\title{In-migration Field Gini Index}
\usage{
  migration.gini.in(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  In-migration Field Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.in(migration.hyp)   # 0.1000000 0.2500000 0.1666667
migration.gini.in(migration.hyp2)  # 0.1000000 0.0000000 0.2142857
}
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997) Measuring Spatial Focusing in a Migration System.
  \emph{Demography} \bold{34}, 251--262 \item M. Bell, M.
  Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell
  and G. Hugo (2002) Cross-National Comparison of Internal
  Migration. Issues and Measures. \emph{Journal of the
  Royal Statistical Society. Series A (Statistics in
  Society)} \bold{165}, 435--464 }
}

