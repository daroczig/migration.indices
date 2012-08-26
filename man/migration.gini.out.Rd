\name{migration.gini.out}
\alias{migration.gini.out}
\title{Out-migration Field Gini Index}
\usage{
  migration.gini.out(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  Out-migration Field Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.out(migration.hyp)  # 0 0 0
migration.gini.out(migration.hyp2) # 0.000 0.125 0.000
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

