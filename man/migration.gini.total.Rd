\name{migration.gini.total}
\alias{migration.gini.total}
\title{Total Flows Gini Index}
\usage{
  migration.gini.total(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  number
}
\description{
  Total Flows Gini Index
}
\examples{
\dontrun{
data(migration.hyp)
migration.gini.total(migration.hyp)        # 0.2222222
migration.gini.total(migration.hyp2)       # 0.1875
}
}
\author{
  Gergely Daróczi
}
\references{
  \itemize{ \item David A. Plane and Gordon F. Mulligan
  (1997): Measuring Spatial Focusing in a Migration System.
  In. Demography, Vol. 34, No. 2 (May, 1997), pp. 251-262
  \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P.
  Rees, J. Stillwell and G. Hugo (2002): Cross-National
  Comparison of Internal Migration. Issues and Measures.
  In. Journal of the Royal Statistical Society. Series A
  (Statistics in Society), Vol. 165, No. 3 (2002), pp.
  435-464 }
}

