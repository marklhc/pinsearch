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
    syn <- lapply(seq_len(nrow(ind_matrix)), function(t) {
        paste0(
            "# Time ", t, "\n",
            one_factor_syntax(ind_matrix[, t],
                lv_name = lv_names[t],
                load_lab = load_labels[, t],
                int_lab = int_labels[, t],
                uniq_lab = uniq_labels[, t]
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
