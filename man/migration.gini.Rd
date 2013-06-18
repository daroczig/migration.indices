\name{migration.gini}
\alias{migration.gini}
\title{Spatial Gini Indexes}
\usage{
  migration.gini(m, corrected = TRUE)
}
\arguments{
  \item{m}{migration matrix}

  \item{corrected}{to use Bell et al. (2002) updated
  formulas instead of Plane and Mulligan (1997)}
}
\value{
  List of all Gini indices.
}
\description{
  This is a wrapper function computing all the following
  Gini indices: \itemize{ \item Total Flows Gini Index
  (\code{\link{migration.gini.total}}) \item Rows Gini
  Index (\code{\link{migration.gini.row}}) \item
  Standardized Rows Gini Index
  (\code{\link{migration.gini.row.standardized}}) \item
  Columns Gini Index (\code{\link{migration.gini.col}})
  \item Standardized Columns Gini Index
  (\code{\link{migration.gini.col.standardized}}) \item
  Exchange Gini Index
  (\code{\link{migration.gini.exchange}}) \item
  Standardized Exchange Gini Index
  (\code{\link{migration.gini.exchange.standardized}})
  \item Out-migration Field Gini Index
  (\code{\link{migration.gini.out}}) \item
  Migration-weighted Out-migration Gini Index
  (\code{\link{migration.weighted.gini.out}}) \item
  In-migration Field Gini Index
  (\code{\link{migration.gini.in}}) \item
  Migration-weighted In-migration Gini Index
  (\code{\link{migration.weighted.gini.in}}) \item
  Migration-weighted Mean Gini Index
  (\code{\link{migration.weighted.gini.mean}}) }
}
\examples{
data(migration.hyp)
migration.gini(migration.hyp)
migration.gini(migration.hyp2)
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
\seealso{
  \code{\link{migration.gini.col}}
  \code{\link{migration.gini.row}}
  \code{\link{migration.gini.exchange}}
  \code{\link{migration.gini.in}}
  \code{\link{migration.gini.out}}
}

