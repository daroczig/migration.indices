\name{migration.connectivity}
\alias{migration.connectivity}
\title{Migration Connectivity Index}
\usage{
  migration.connectivity(m)
}
\arguments{
  \item{m}{migration matrix}
}
\value{
  A number between 0 and 1 where zero shows no connections
  between regions.
}
\description{
  The Migration Connectivity Index measures "the proportion
  of the total number of potential interregional flows
  which are not zero": \deqn{I_{MC} = \sum_i \sum_{j \neq
  i} \frac{MC_{ij}}{n(n-1)}} where \eqn{MC_{ij}} is 0 if
  the flow from \eqn{i} to \eqn{j} is zero and let it be 1
  otherwise.
}
\examples{
data(migration.hyp)
migration.connectivity(migration.hyp)
data(migration.world)
migration.connectivity(migration.world)
}
\references{
  \itemize{ \item M. Bell, M. Blake, P. Boyle, O.
  Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002)
  Cross-National Comparison of Internal Migration. Issues
  and Measures. \emph{Journal of the Royal Statistical
  Society. Series A (Statistics in Society)} \bold{165},
  435--464 }
}

