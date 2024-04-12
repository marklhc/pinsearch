type2op <- function(type) {
  switch(type,
         loadings = "=~",
         intercepts = "~1",
         thresholds = "|",
         residuals = "~~",
         residual.covariances = "~~")
}
# type2col <- function(type) {
#   switch(type,
#          loadings = "rhs",
#          intercepts = "lhs",
#          residuals = "lhs")
# }
op2col <- function(op) {
    switch(op,
           `=~` = "rhs",
           `~1` = "lhs",
           `~~` = "lhs")
}

initialize_partable <- function(pt0, ninv_items) {
    ninv_lds <- ninv_items[ninv_items$type == "loadings", ]
    # Augmented ninv_items (intercepts are freed if loadings are free)
    if (nrow(ninv_lds) > 0) {
        ninv_items <- rbind(
            ninv_items,
            data.frame(
                lhs = ninv_lds$rhs,
                rhs = "",
                group = ninv_lds$group,
                type = "intercepts"
            )
        )
    }
    for (j in seq_len(nrow(ninv_items))) {
        pt0 <- remove_cons(pt0,
            lhs = ninv_items$lhs[j], rhs = ninv_items$rhs[j],
            group = ninv_items$group[j], op = type2op(ninv_items$type[j])
        )
    }
    ninv_res <- ninv_items[ninv_items$type == "residuals", ]
    for (j in seq_len(nrow(ninv_res))) {
        # Get items covaried with item j
        cov_items1 <- pt0[pt0$op == "~~" &
            pt0$lhs == ninv_res$lhs[j] &
            pt0$group == ninv_res$group[j] &
            pt0$lhs != pt0$rhs, ]
        cov_items2 <- pt0[pt0$op == "~~" &
            pt0$rhs == ninv_res$lhs[j] &
            pt0$group == ninv_res$group[j] &
            pt0$lhs != pt0$rhs, ]
        cov_items <- rbind(cov_items1, cov_items2)
        for (k in seq_len(nrow(cov_items))) {
            pt0 <- remove_cons(pt0,
                lhs = cov_items$lhs[k],
                rhs = cov_items$rhs[k],
                group = cov_items$group[k], op = "~~"
            )
        }
    }
    pt0$id <- seq_len(nrow(pt0))
    pt0
}

fdr_alpha <- function(i, m, q = .05) {
    # Based on Benjamini & Gavrilov (2009)
    i * q / (m + 1 - i * (1 - q))
}

cfa2 <- function(...) {
    # Needed for executing `do.call(cfa2, ...)`
    lavaan::cfa(...)
}

#' Search for noninvariant parameters across groups.
#'
#' The function implements the sequential selection method similar to
#' that discussed in
#' [Yoon and Millsap (2007)](https://doi.org/10.1080/10705510701301677).
#' The function proceeds in the order of metric, scalar (threshold),
#' and strict invariance. In each stage, invariance constraints in all
#' items, and the constraint associated with the biggest test
#' statistic above a predefined threshold is freed, before recomputing
#' the test statistic for the next constraint to free.
#'
#' @details Note that when an item has a non-invariant loading, the
#'   corresponding intercept constraint will automatically be freed,
#'   as intercept difference across groups is sensitive to the location
#'   of the zero point for the latent variable and the item.
#'
#' For a particular stage of invariance constraints, the Benjamini &
#'   Gavrilov method uses an adjusted \eqn{alpha} level of
#'   \deqn{iq / [m + 1 - i(1 - q)]}
#'   where \eqn{i} is the step index in the search, \eqn{m} is the
#'   maximum number of constraints that can be freed, and \eqn{q} is
#'   the desirable significance level.
#'
#' @param config_mod Syntax of a configural invariance model to be passed to
#'   [lavaan::cfa()].
#' @param ... Additonal arguments passed to [lavaan::cfa()].
#' @param type Character variable indicating the stage of invariance to be
#'   searched. Currently supported options are (a) for continuous indicators,
#'   "loadings", "intercepts", "residuals", and "residual.covariances", and
#'   (b) "loadings", "thresholds", "residual.covariances", in an increasingly
#'   strict order. A stricter model (e.g., "residual.covariances") will have
#'   constraints of all previous stages.
#' @param inv_test Character variable indicating the statistical test to be
#'   used for specification search. Currently supported options are `"mod"` for
#'   modification index using `lavaan::modindices()`, `"score"` for score
#'   test statistic using `lavaan::lavTestScore()`, and (experimental) `"lrt"`
#'   for likelihood ratio test statistic using `lavaan::lavTestLRT()`.
#' @param sig_level Significance level used to determine whether the parameter
#'   associated with the highest modification index should be removed. Default
#'   is .05.
#' @param control_fdr Logical; whether to use adjust for false discovery rate
#'   for multiple testing. If `TRUE`, the method by Benjamini & Gavrilov (2009)
#'   will be used.
#' @param effect_size Logical; whether to compute dmacs (two groups) or
#'   fmacs (> two groups) effect size or not (default).
#'   This is an experimental feature.
#' @param progress Logical; an experimental feature of showing a progress bar
#'   if `TRUE`. Because the number of steps is unknown until the stopping
#'   criteria are reached, the progress bar may be inaccurate.
#'
#' @return A list of three elements:
#' \itemize{
#'   \item{`Partial Invariance Fit`}{A [`lavaan::lavaan-class`]
#'     object containing the final partial invariance model.}
#'   \item{`Non-Invariant Items`}{A data frame of non-invariant
#'     parameters.}
#'   \item{`effect_size`}{Effect size statistics obtained from
#'     [pin_es()].}
#' }
#' @references Yoon, M., & Millsap, R. E. (2007). Detecting violations of
#'   factorial invariance using data-based specification searches: A
#'   Monte Carlo study. Structural Equation Modeling: A Multidisciplinary
#'   Journal, 14(3), 435-463.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' library(lavaan)
#' library(MASS)
#' # Simulate random data
#' set.seed(14) # for reproducible results
#' mod_ninv <- "f =~ c(1, 0.8, .6)*y1 + c(0.8, 1.2, 1.0)*y2"
#' mod_inv <- paste0("1*y", 3:6, collapse = " + ")
#' sim_mod <- paste(
#'     paste(mod_ninv, "+", mod_inv),
#'     "f ~~ c(1, 1.3, 1.5)*f
#'      f ~ c(0, 0.5, 1.0)*1
#'      y1 ~ c(0, 0.3, 0)*1
#'      y3 ~ c(0.3, 0, -0.3)*1
#'      y1 ~~ c(1, .5, 1)*y1",
#'     sep = "\n"
#' )
#' # The uniqueness of each item is assumed to be 1.0
#' dat_sim <- simulateData(sim_mod, sample.nobs = c(100, 100, 100))
#' # Fit configural model:
#' sam_config <- paste(
#'     paste0("f =~ ", paste0("y", 1:6, collapse = " + "))
#' )
#' pinSearch(sam_config,
#'     data = dat_sim, group = "group",
#'     type = "intercepts"
#' )
#' @export
pinSearch <- function(config_mod,
                       ...,
                       type = c(
                           "loadings", "intercepts", "thresholds",
                           "residuals", "residual.covariances"
                       ),
                       inv_test = c("mod", "score", "lrt"),
                       sig_level = .05,
                       control_fdr = FALSE,
                       effect_size = FALSE,
                       progress = FALSE) {
    type <- match.arg(type)
    inv_test <- match.arg(inv_test)
    dots <- list(...)
    base_fit <- do.call(cfa2,
        args = c(list(model = config_mod, do.fit = FALSE), dots)
    )
    base_call <- stats::getCall(base_fit)
    # base_call <- stats::getCall(base_fit)
    base_opt <- base_fit@Options
    # stored_fit <- list()
    ind_names <- get_ovnames(base_fit)
    if (is.null(base_call$std.lv)) {
        dots$std.lv <- TRUE
        if (interactive()) message("`std.lv` is set to TRUE by default")
    }
    if (base_opt$categorical) {
        types <- c("loadings", "thresholds", "residual.covariances")
        if (is.null(base_call$parameterization)) {
            dots$parameterization <- "theta"
        }
        if (dots$parameterization == "theta") {
            # add constraints on unique variances
            config_mod <- c(config_mod,
                            paste(ind_names, "~~ 1 *", ind_names))
            if (interactive()) {
                message("Unique variances are constrained to 1 ",
                        "for identification")
            }
        } else if (dots$parameterization == "delta") {
            # add constraints on unique variances
            config_mod <- c(config_mod,
                            paste(ind_names, "~*~ 1 *", ind_names))
            warning(
                "Delta parameterization has not been tested ",
                "and likely results in untrustworthy results."
            )
        }
        if (!type %in% types) {
            stop("`type = ", type, "` cannot be used with ordered items")
        }
    } else {
        types <- c(
            "loadings", "intercepts", "residuals",
            "residual.covariances"
        )
        if (!type %in% types) {
            stop("`type = ", type, "` cannot be used with continuous items")
        }
    }
    base_fit <- do.call(cfa2,
        args = c(list(model = config_mod), dots)
    )
    n_type <- which(types == type) # number of stages
    fn_get_inv <- switch(inv_test,
        mod = get_invmod,
        score = get_invscore,
        lrt = get_invlrt
    )
    ninv_items <- data.frame(
        lhs = character(),
        rhs = character(),
        group = numeric(),
        type = character(),
        stringsAsFactors = FALSE
    )
    for (i in seq_len(n_type)) {
        typei <- types[i]
        if (progress) {
            message(
                "\n[", i, "/", n_type, "] Searching for ",
                typei, " noninvariance\n"
            )
        }
        # new_fit <- lavaan::update(base_fit, ..., config_mod,
        #                           group.equal = types[seq_len(i)],
        #                           do.fit = i <= 1)
        new_fit <- do.call(
            cfa2,
            c(list(
                model = config_mod, group.equal = types[seq_len(i)],
                do.fit = i <= 1
            ), dots)
        )
        pt0 <- lavaan::parTable(new_fit)
        if (i >= 2) {
            pt0 <- initialize_partable(pt0, ninv_items = ninv_items)
            # new_fit <- lavaan::update(new_fit, ..., pt0,
            #     group.equal = types[seq_len(i)],
            #     do.fit = TRUE
            # )
            new_fit <- do.call(
                cfa2,
                c(list(model = pt0, group.equal = types[seq_len(i)]), dots)
            )
        }
        if (typei == "residual.covariances" &&
            base_fit@test$standard$df >= new_fit@test$standard$df) {
            next
        }
        lrt_base_new <- lavaan::lavTestLRT(base_fit, new_fit)
        df_diff <- lrt_base_new[2, "Df diff"]
        if (isTRUE(all.equal(lrt_base_new[2, "Chisq diff"], 0)) &&
            df_diff == 0) {
                message("All free ", typei, " are noninvariant!")
                next
            }
        if (lrt_base_new[2, "Pr(>Chisq)"] >= sig_level) {
            base_fit <- new_fit
            next # skip to next stage
        } else {
            op <- type2op(typei)
            if (!control_fdr) {
                p_enter <- sig_level
            } else {
                num_free <- 1
                p_enter <- fdr_alpha(num_free, m = df_diff, q = sig_level)
            }
            row_to_free <- do.call(
                fn_get_inv,
                c(list(
                    object = new_fit, type = typei, alpha = p_enter,
                    ind_names = ind_names, group.equal = types[seq_len(i)]
                ), dots)
            )
            if (progress) {
                total_mod <- remain_mod <- attr(row_to_free, which = "size")
                pb <- txtProgressBar(min = 0, max = total_mod, style = 3)
                pb_count <- 0
            }
            while (!is.null(row_to_free)) {
                pt0 <- do.call(remove_cons,
                    args = c(list(pt0), row_to_free[c("group", "lhs", "rhs")],
                        op = op
                    )
                )
                ninv_items <- rbind(
                    ninv_items,
                    c(row_to_free[c("group", "lhs", "rhs")], type = typei)
                )
                if (progress) {
                    pb_count <- max(total_mod - remain_mod + 1, pb_count + 1)
                    setTxtProgressBar(pb, pb_count)
                }
                new_fit <- do.call(
                    cfa2,
                    c(list(model = pt0, group.equal = types[seq_len(i)]), dots)
                )
                # new_fit <- fit_cfa(pt0, eq = types[seq_len(i)], ...)
                if (control_fdr) {
                    num_free <- num_free + 1
                    p_enter <- fdr_alpha(num_free, m = df_diff, q = sig_level)
                }
                row_to_free <- do.call(
                    fn_get_inv,
                    c(list(
                        object = new_fit, type = typei, alpha = p_enter,
                        ind_names = ind_names, group.equal = types[seq_len(i)]
                    ), dots)
                )
                remain_mod <- attr(row_to_free, which = "size")
            }
            base_fit <- new_fit
        }
    }
    if (progress) close(pb)
    out <- list(
        `Partial Invariance Fit` = new_fit,
        `Non-Invariant Items` = ninv_items
    )
    if (effect_size) {
        out$effect_size <- es_lavaan(new_fit)
    }
    out
}
