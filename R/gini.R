#' Check Migration Matrix
#'
#' Checks if provided R object looks like a migration matrix.
#'
#' A migration matrix is a symmetric matrix with the exact same row and column names. The diagonal equals to zero. The upper triangle shows the in- and the lower triangle shows the out-migration.
#' @param m R object to check
#' @return (invisibly) TRUE
#' @keywords internal
check.migration.matrix <- function(m) {

    ## dummy checks on provided matrix
    if (missing(m))
        stop('No data provided!')
    if (!is.matrix(m))
        stop('Wrong data type (!matrix) provided!')
    if (nrow(m) != ncol(m))
        stop('Wrong data tpye (!symmetrical matrix) provided!')
    if (!is.numeric(m))
        stop('Wrong data tpye (!numeric) provided!')
    if (length(which(is.na(m[xor(upper.tri(m), lower.tri(m))]))) > 0)
        stop('Missing values (outside of diagonal) found in provided matrix!')
    if (!any(all(is.na(diag(m))), all(diag(m) == 0)))
        stop('Diagonal should be zero or missing.')

    return(invisible(TRUE))

}


#' Total Flows Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references \itemize{
#' \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.total(migration.hyp)        # 0.2222222
#' migration.gini.total(migration.hyp2)       # 0.1875
#' }
#' @export
migration.gini.total <- function(m) {

    check.migration.matrix(m)

    n           <- nrow(m)
    m.val       <- m[xor(upper.tri(m), lower.tri(m))]
    return(sum(apply(as.data.frame(m.val), 1, function(x) sum(abs(m.val-x))), na.rm = TRUE)/(2*n*(n-1)*sum(m, na.rm = TRUE)))

    ## faster method (fails with "low memory")
    diag(m)     <- NA
    sum(dist(as.vector(m)), na.rm = TRUE)*2/((2*n*(n-1)-1)*sum(m, na.rm = TRUE))

}


#' Rows Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.row(migration.hyp)  # 0
#' migration.gini.row(migration.hyp2) # 0.02083333
#' }
#' @export
migration.gini.row <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    sum(apply(m, 1, function(m.row) sum(dist(m.row), na.rm = TRUE)*2))/(2*n*(n-1)*sum(m, na.rm = TRUE))

}


#' Standardized Rows Gini Index
#'
#' @param m migration matrix
#' @param migration.gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return percentage
#' @references David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.row.standardized(migration.hyp)     # 0
#' migration.gini.row.standardized(migration.hyp2)    # 11.11111
#' }
#' @export
migration.gini.row.standardized <- function(m, migration.gini.total = migration.gini.total(m)) {

    100 * migration.gini.row(m) / migration.gini.total
}


#' Columns Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.col(migration.hyp)  # 0.05555556
#' migration.gini.col(migration.hyp2) # 0.04166667
#' }
#' @export
migration.gini.col <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    sum(apply(m, 2, function(m.row) sum(dist(m.row), na.rm = TRUE)*2))/(2*n*(n-1)*sum(m, na.rm = TRUE))

}


#' Standardized Columns Gini Index
#'
#' @param m migration matrix
#' @param migration.gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return percentage
#' @references David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.col.standardized(migration.hyp)     # 25
#' migration.gini.col.standardized(migration.hyp2)    # 22.22222
#' }
#' @export
migration.gini.col.standardized <- function(m, migration.gini.total = migration.gini.total(m)) {

    100 * migration.gini.col(m) / migration.gini.total
}


#' Exchange Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.exchange(migration.hyp)     # 0.05555556
#' migration.gini.exchange(migration.hyp2)    # 0.04166667
#' }
migration.gini.exchange <- function(m) {

    check.migration.matrix(m)

    n           <- nrow(m)
    m.t         <- t(m)
    m.val       <- m[xor(upper.tri(m), lower.tri(m))]
    m.t.val     <- m.t[xor(upper.tri(m.t), lower.tri(m.t))]

    sum(abs(m.val - m.t.val))/(2*n*(n-1)*sum(m))

}


#' Standardized Exchange Gini Index
#'
#' @param m migration matrix
#' @param migration.gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return number
#' @references David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.exchange.standardized(migration.hyp)  # 25
#' migration.gini.exchange.standardized(migration.hyp2) # 22.22222
#' }
#' @export
migration.gini.exchange.standardized <- function(m, migration.gini.total = migration.gini.total(m)) {

    100 * migration.gini.exchange(m) / migration.gini.total

}


#' Out-migration Field Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references \itemize{
#' \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.out(migration.hyp)  # 0 0 0
#' migration.gini.out(migration.hyp2) # 0.000 0.125 0.000
#' }
migration.gini.out <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    apply(m, 1, function(m.row) sum(dist(m.row), na.rm = TRUE) * 2) / (2 * (n - 2) * rowSums(m, na.rm = TRUE))

}


#' Migration-weighted Out-migration Gini Index
#'
#' @param m migration matrix
#' @param mgi optionally passed (precomputed) Migration In-migration Gini Index
#' @return number
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.weighted.gini.out(migration.hyp)   #
#' migration.weighted.gini.out(migration.hyp2)  #
#' }
#' @seealso \code{\link{migration.weighted.gini.in}} \code{\link{migration.weighted.gini.mean}}
#' @export
migration.weighted.gini.out <- function(m, mgo) {

    diag(m)     <- NA
    n           <- nrow(m)
    m.sum       <- sum(m, na.rm = TRUE)

    if (missing(mgo))
        mgo <- migration.gini.out(m)

    sum(mgo * colSums(m, na.rm = TRUE) / m.sum) / n

}


#' In-migration Field Gini Index
#'
#' @param m migration matrix
#' @return number
#' @references \itemize{
#' \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini.in(migration.hyp)   # 0.1000000 0.2500000 0.1666667
#' migration.gini.in(migration.hyp2)  # 0.1000000 0.0000000 0.2142857
#' }
#' @export
migration.gini.in <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    apply(m, 2, function(m.row) sum(dist(m.row), na.rm = TRUE) * 2) / (2 * (n - 2) * colSums(m, na.rm = TRUE))

}


#' Migration-weighted In-migration Gini Index
#'
#' @param m migration matrix
#' @param mgi optionally passed (precomputed) Migration In-migration Gini Index
#' @return number
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.weighted.gini.in(migration.hyp)   # 0.1222222
#' migration.weighted.gini.in(migration.hyp2)  # 0.05238095
#' }
#' @seealso \code{\link{migration.weighted.gini.out}} \code{\link{migration.weighted.gini.mean}}
#' @export
migration.weighted.gini.in <- function(m, mgi) {

    diag(m)     <- NA
    n           <- nrow(m)
    m.sum       <- sum(m, na.rm = TRUE)

    if (missing(mgi))
        mgi <- migration.gini.in(m)

    sum(mgi * rowSums(m, na.rm = TRUE) / m.sum) / n

}

#' Migration-weighted Mean Gini Index
#'
#' @param m migration matrix
#' @param mwgi optionally passed (precomputed) Migration-weighted In-migration Gini Index
#' @param mwgo optionally passed (precomputed) Migration-weighted Out-migration Gini Index
#' @return number
#' @references M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.weighted.gini.mean(migration.hyp)  # 0.06111111
#' migration.weighted.gini.mean(migration.hyp2) # 0.03660714
#' }
#' @seealso \code{\link{migration.weighted.gini.in}} \code{\link{migration.weighted.gini.out}}
#' @export
migration.weighted.gini.mean <- function(m, mwgi, mwgo) {

    if (missing(mwgi))
        mwgi <- migration.weighted.gini.in(m)
    if (missing(mwgo))
        mwgo <- migration.weighted.gini.out(m)

    (mwgi + mwgo) / 2

}


#' Spatial Gini Indexes
#'
#' @param m migration matrix
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.gini(migration.hyp)
#' migration.gini(migration.hyp2)
#' }
#' @export
migration.gini <- function(m) {

    res <- list(
             migration.gini.total         = migration.gini.total(m),
             migration.gini.exchange      = migration.gini.exchange(m),
             migration.gini.row           = migration.gini.row(m),
             migration.gini.col           = migration.gini.col(m),
             migration.gini.in            = migration.gini.in(m),
             migration.gini.out           = migration.gini.out(m)
      )

    res$migration.gini.row.standardized           <- migration.gini.row.standardized(m, res$migration.gini.total)
    res$migration.gini.col.standardized           <- migration.gini.col.standardized(m, res$migration.gini.total)
    res$migration.gini.exchange.standardized      <- migration.gini.exchange.standardized(m, res$migration.gini.total)
    res$migration.gini.in.weighted                <- migration.weighted.gini.in(m, res$migration.gini.in)
    res$migration.gini.out.weighted               <- migration.weighted.gini.out(m, res$migration.gini.out)
    res$migration.gini.mean.weighted              <- migration.weighted.gini.mean(m, res$migration.gini.in.weighted, res$migration.gini.out.weighted)

    class(res) <- 'migration.gini'
    return(res)

}


#' @method print migration.gini
#' @S3method print migration.gini
print.migration.gini <- function(x) {

    cat('\n')
    cat('Total Flows Gini Index:             ', x$migration.gini.total, '\n')
    cat('Rows Gini Index:                    ', x$migration.gini.row, '\n')
    cat('Standardized Rows Gini Index:       ', x$migration.gini.row.standardized, '\n')
    cat('Columns Gini Index:                 ', x$migration.gini.col, '\n')
    cat('Standardized Columns Gini Index:    ', x$migration.gini.col.standardized, '\n')
    cat('Exchange Gini Index:                ', x$migration.gini.exchange, '\n')
    cat('Standardized Exchange Gini Index:   ', x$migration.gini.exchange.standardized, '\n')
    cat('In-migration Field Gini Index:      ', 'vector', '\n')
    cat('Weighted In-migration Gini Index:   ', x$migration.gini.in.weighted, '\n')
    cat('Out-migration Field Gini Index:     ', 'vector', '\n')
    cat('Weighted Out-migration Gini Index:  ', x$migration.gini.out.weighted, '\n')
    cat('Migration-weighted Mean Gini Index: ', x$migration.gini.mean.weighted, '\n')
    cat('\n')

}
