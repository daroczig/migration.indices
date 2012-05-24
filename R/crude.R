#' Crude Migration Rate
#'
#'
#' @param m migration matrix
#' @param PAR population at risk (estimated average population size)
#' @param k scaling constant (set to \code{100} by default to result in percentage)
#' @return percentage (when \code{k=100})
#' @author Gergely Daróczi
#' @references Philip Rees, Martin Bell, Oliver Duke-Williams and Marcus Blake (2000): Problems and Solutions in the Measurement of Migration Intensities: Australia and Britain Compared. In. Population Studies, Vol. 54, No. 2 (Jul., 2000), pp. 207-222
#' @examples \dontrun{
#' data(migration.world)
#' migration.cmr(migration.world, 6e+9)
#' }
#' @export
migration.cmr <- function(m, PAR, k = 100) {

    check.migration.matrix(m)
    if (missing(PAR))
        stop('Estimated population size was not provided!')

    k*sum(m)/PAR

}
