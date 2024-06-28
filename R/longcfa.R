#' Wrapper function for longitudinal CFA (Experimental)
#'
#' This function estimates a longitudinal measurement model by first
#' generating model syntax (using the [longcfa_syntax()] function), and
#' then fit the model using [lavaan::cfa()]. An input matrix is needed
#' where the rows are the indicators and the columns are the time points.
#'
#' Currently only supports model with one latent variable at each time
#' point.
#'
#' @param ind_matrix A $p \times T$ character matrix specifying the
#'   names of the indicator variables across time points. Each column
#'   corresponds to a time point.
#' @param lv_names A vector of names of $T$ latent variables.
#' @param model A character string showing additional syntax to be
#'   added to the model. Defaults to `NULL`.
#' @param lag_cov Logical; whether the same indicator is allowed to
#'   correlate over time.
#' @param long_equal A character vector indicating types of parameters
#'   to be constrained equal across time points. This is similar to
#'   the `group.equal` argument in `lavaan::cfa()`. Currently, only
#'   `"loadings"`, `"intercepts"`, and `"residuals"` are supported.
#' @param long_partial A named list of matrices specifying specific
#'   indicators to have different parameter values for a specific time
#'   point. The list should have names "loadings", "intercepts", and
#'   "residuals", and each element is a 2-column matrix where each row
#'   specifies the indicator (column 1) and the time point (column 2)
#'   for the parameter to be free. See example below.
#' @param ... Other arguments passed to `lavaan::cfa()`, such as `data`.
#'
#' @return A fit object as returned by `lavaan::cfa()`.
#'
#' @examples
#' library(lavaan)
#' # Indicator matrix
#' spec <- matrix(c(
#'     "y1", "y2", "y3", "y4",
#'     "y5", "y6", "y7", "y8"
#' ), ncol = 2)
#' # Scalar invariance
#' fit <- longcfa(spec,
#'                lv_names = c("dem60", "dem65"),
#'                data = PoliticalDemocracy,
#'                long_equal = c("loadings", "intercepts"))
#' summary(fit)
#' # Partial invariance
#' fit2 <- longcfa(spec,
#'                 lv_names = c("dem60", "dem65"),
#'                 data = PoliticalDemocracy,
#'                 long_equal = c("loadings", "intercepts"),
#'                 long_partial = list(
#'                     loadings = matrix(c(1, 2), ncol = 2),
#'                     intercepts = matrix(c(1, 3, 2, 2), ncol = 2)
#'                 ))
#' summary(fit2)
#' @export
longcfa <- function(ind_matrix, lv_names, model = NULL,
                    lag_cov = FALSE,
                    long_equal = NULL,
                    long_partial = NULL, ...) {
    syn <- longcfa_syntax(ind_matrix,
        lv_names = lv_names,
        lag_cov = lag_cov,
        long_equal = long_equal,
        long_partial = long_partial
    )
    lavaan::cfa(syn, ...,
        auto.fix.first = FALSE,
        int.lv.free = TRUE
    )
}

# TODO: Accomodate input pattern matrix for multiple latent variables
# TODO: Need to handle NA

longcfa_syntax <- function(ind_matrix, lv_names, lag_cov = FALSE,
                           long_equal = NULL, long_partial = NULL) {
    if ("loadings" %in% long_equal) {
        load_labels <- gen_labels(dim(ind_matrix), ".l",
                                  partial = long_partial$loadings)
        latvar_syntax <- latent_var_syntax(lv_names[1])
    } else {
        load_labels <- NULL
        latvar_syntax <- latent_var_syntax(lv_names)
    }
    if ("intercepts" %in% long_equal) {
        int_labels <- gen_labels(dim(ind_matrix), ".i",
                                 partial = long_partial$intercepts)
        latmean_syntax <- latent_mean_syntax(lv_names[1])
    } else {
        int_labels <- NULL
        latmean_syntax <- latent_mean_syntax(lv_names)
    }
    if ("residuals" %in% long_equal) {
        uniq_labels <- gen_labels(dim(ind_matrix), ".u",
                                  partial = long_partial$residuals)
    } else {
        uniq_labels <- NULL
    }
    syn <- lapply(seq_len(ncol(ind_matrix)), function(t) {
        valid_pos <- which(!is.na(ind_matrix[, t]))
        paste0(
            "# Time ", t, "\n",
            one_factor_syntax(ind_matrix[valid_pos, t, drop = FALSE],
                lv_name = lv_names[t],
                load_lab = load_labels[valid_pos, t, drop = FALSE],
                int_lab = int_labels[valid_pos, t, drop = FALSE],
                uniq_lab = uniq_labels[valid_pos, t, drop = FALSE]
            ),
            "\n"
        )
    })
    if (lag_cov) {
        lagcov_syntax <- lag_cov_syntax(ind_matrix)
    } else {
        lagcov_syntax <- NULL
    }
    paste(c(syn,
            paste0(latvar_syntax, "\n"),
            paste0(latmean_syntax, "\n"),
            lagcov_syntax),
          collapse = "\n")
}

one_factor_syntax <- function(ind_names, lv_name = ".eta",
                              load_lab = NULL, int_lab = NULL,
                              uniq_lab = NULL) {
    if (is.null(load_lab)) {
        load_ind <- ind_names
    } else {
        load_ind <- paste(load_lab, ind_names, sep = " * ")
    }
    out <- paste0(lv_name, " =~ ", paste0(load_ind, collapse = " + "))
    if (!is.null(int_lab)) {
        out <- c(out, paste(ind_names, "~", int_lab, "* 1"))
    }
    if (!is.null(uniq_lab)) {
        out <- c(out, paste(ind_names, "~~", uniq_lab, "*", ind_names))
    }
    paste(out, collapse = "\n")
}

lag_cov_syntax <- function(ind_matrix) {
    out <- "# Lag Covariances"
    for (i in seq_len(nrow(ind_matrix))) {
        ind_names_i <- ind_matrix[i, ]
        ind_names_i <- ind_names_i[!is.na(ind_names_i)]
        p_i <- length(ind_names_i)
        for (j in seq_len(p_i - 1)) {
            out <- c(out, paste(ind_names_i[j], "~~",
                                paste(ind_names_i[(j + 1):p_i],
                                      collapse = " + ")))
        }
    }
    paste(out, collapse = "\n")
}

gen_labels <- function(size, prefix, equal = TRUE,
                       partial = NULL) {
    if (equal) {
        out <- matrix(seq_len(size[1]), nrow = size[1], ncol = size[2])
        if (!is.null(partial)) {
            partial <- as.matrix(partial)
            out[partial] <- apply(partial, MARGIN = 1,
                                  FUN = paste0, collapse = "")
        }
    } else {
        out <- outer(seq_len(size[1]),
            Y = seq_len(size[2]),
            FUN = paste0
        )
    }
    out[] <- paste0(prefix, out)
    out
}

latent_var_syntax <- function(lv_names) {
    paste(
        c(
            "# Latent variances",
            paste(lv_names, "~~ 1 *", lv_names)
        ),
        collapse = "\n"
    )
}

latent_mean_syntax <- function(lv_names) {
    paste(
        c(
            "# Latent means",
            paste(lv_names, "~ 0 * 1")
        ),
        collapse = "\n"
    )
}

sub_vars <- function(vnames, replacement, x) {
    for (i in seq_along(vnames)) {
        x <- gsub(vnames[i], replacement = replacement[i], x = x)
    }
    x
}
