#' In-migration Coefficient of Variation
#'
#' @param m migration matrix
#' @return number
#' @references \itemize{
#' \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.cv.in(migration.hyp)    # 0.2000000 0.5000000 0.3333333
#' migration.cv.in(migration.hyp2)   # 0.2000000 0.0000000 0.4285714
#' }
#' @export
migration.cv.in <- function(m) {
    diag(m) <- NA
    n       <- ncol(m) - 1
    apply(m, 2, function(m.row) (sd(m.row, na.rm = TRUE) * sqrt((n - 1) / n)) / mean(m.row, na.rm = TRUE))
}


#' Out-migration Coefficient of Variation
#'
#' @param m migration matrix
#' @return number
#' @references \itemize{
#' \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.cv.out(migration.hyp)    # 0 0 0
#' migration.cv.out(migration.hyp2)   # 0.00 0.25 0.00
#' }
#' @export
migration.cv.out <- function(m) {
    diag(m) <- NA
    n       <- ncol(m) - 1
    apply(m, 1, function(m.row) (sd(m.row, na.rm = TRUE) * sqrt((n - 1) / n)) / mean(m.row, na.rm = TRUE))
}
