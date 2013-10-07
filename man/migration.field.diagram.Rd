\name{migration.field.diagram}
\alias{migration.field.diagram}
\title{Joint plot for in and out-migration fields}
\usage{
  migration.field.diagram(m, method = c("gini", "acv"),
    title = "Migration field diagram",
    xlab = "Out-migration", ylab = "In-migration")
}
\arguments{
  \item{m}{migration matrix}

  \item{method}{measurement of in and out-migration}

  \item{title}{plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}
}
\description{
  This migration field diagram makes easy to visualize both
  direction of migration. E.g. points above the diagonal
  "are outward redistributors, while those below that line
  are inward redistributors."
}
\examples{
\dontrun{
data(migration.world)
par(mfrow = c(2, 1))
migration.field.diagram(migration.world)
migration.field.diagram(migration.world, method = 'acv')
}
}
\references{
  \itemize{ \item Source code was adopted from Michael Ward
  and Kristian Skrede Gleditsch (2008) \emph{Spatial
  Regression Models}. Thousand Oaks, CA: Sage.
  \url{http://privatewww.essex.ac.uk/~ksg/code/srm_enhanced_code_v5.R}
  with the permission of the authors.  \item Case study and
  use case: Andrei Rogers and Stuart Sweeney (1998)
  Measuring the Spatial Focus of Migration Patterns.
  \emph{The Professional Geographer} \bold{50}, 232--242 }
}

