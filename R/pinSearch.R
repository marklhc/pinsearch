lavaanify_cfa <- function(config_mod, ngroups, group.equal, varTable,
                          parameterization = "theta") {
    lavaan::lavaanify(config_mod, meanstructure = TRUE, ngroups = ngroups,
                      parameterization = parameterization,
                      varTable = varTable,
                      group.equal = group.equal,
                      std.lv = TRUE, int.ov.free = TRUE, int.lv.free = FALSE,
                      auto.var = TRUE, auto.cov.lv.x = TRUE,
                      auto.cov.y = TRUE, auto.th = TRUE, auto.delta = TRUE,
                      auto.efa = TRUE)
}
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
remove_cons <- function(object, lhs, rhs, group, op) {
    # Remove constraints
    # col <- op2col(op)
    # matched <- if (op == "~~")
    #     object$lhs == object$rhs
    # else
    #     TRUE
    # Identify row with the specified parameter
    row_i <- which(object$lhs == lhs & object$rhs == rhs &
                       object$op == op & object$group == group)
    mi_plab <- object$plabel[row_i]
    object$label[row_i] <- ""
    if (length(mi_plab) == 0)
        return(object)
    ind_plabs <-
        object$plabel[object$op == op & object$lhs == lhs & object$rhs == rhs]
    pt_eq <- object[object$op == "==" &
                        object$lhs %in% ind_plabs &
                        object$rhs %in% ind_plabs,]
    # If it appears on rhs, remove it
    pt_new <-
        object[!object$id %in% pt_eq[pt_eq$rhs == mi_plab, "id"],]
    # Replace lhs with the next group (make sure it's not on the table)
    if (mi_plab %in% pt_eq$lhs) {
        pt_new[pt_new$op == "==" &
                   pt_new$lhs == mi_plab, "lhs"] <- ind_plabs[group + 1]
        pt_new <-
            pt_new[!(pt_new$op == "==" & pt_new$lhs == pt_new$rhs),]
    }
    pt_new
}
get_invmod <- function(object, type, alpha = .05, ind_names, pt0) {
    if (missing(pt0)) pt0 <- lavaan::parTable(object)
    op <- type2op(type)
    pt0_eq <- pt0[pt0$label != "" & pt0$free >= 0, ]
    mis <- lavaan::modindices(
        object,
        op = op,
        minimum.value = stats::qchisq(alpha, 1, lower.tail = FALSE),
        sort. = TRUE,
        free.remove = FALSE
    )
    # Only parameters in the original model
    if (type == "loadings") {
        mis_ids <- rownames(mis)[mis$rhs %in% ind_names &
                                     rownames(mis) %in% pt0_eq$id]
    } else if (type %in% c("intercepts", "thresholds")) {
        mis_ids <- rownames(mis)[mis$lhs %in% ind_names &
                                     rownames(mis) %in% pt0_eq$id]
    } else if (type == "residuals") {
        mis_ids <- rownames(mis)[mis$lhs %in% ind_names &
                                     rownames(mis) %in% pt0_eq$id &
                                     mis$lhs == mis$rhs]
    } else if (type == "residual.covariances") {
        mis_ids <- rownames(mis)[mis$lhs %in% ind_names &
                                     mis$rhs %in% ind_names &
                                     rownames(mis) %in% pt0_eq$id &
                                     mis$lhs != mis$rhs]
    }
    if (length(mis_ids) == 0) return(NULL)
    pt_mis <- pt0_eq[match(mis_ids, pt0_eq$id), ]
    # Only parameters that are free
    pt_mis <- pt_mis[pt_mis$free >= 0, ]
    out <- mis[rownames(pt_mis)[1], ]
    attr(out, "size") <- nrow(pt_mis)
    out
}

get_invscore <- function(object, type, alpha = .05, ind_names, pt0, ...) {
    if (missing(pt0)) pt0 <- lavaan::parTable(object)
    op <- type2op(type)
    pt0_eq <- pt0[pt0$label != "" & pt0$free >= 0 & pt0$op == op, ]
    # Find relevant constraints
    pt0_cons <- pt0[pt0$op == "==", ]
    cons_to_test <- which(pt0_cons$lhs %in% pt0_eq$plabel |
        pt0_cons$rhs %in% pt0_eq$plabel)
    mis <- lavaan::lavTestScore(object, release = cons_to_test, ...)
    cons_to_free <- mis$uni[which.max(mis$uni$X2), ]
    if (cons_to_free$p.value > alpha) return(NULL)
    out <- pt0_eq[which(pt0_eq$plabel == cons_to_free$rhs), ]
    attr(out, "size") <- sum(mis$uni$p.value < alpha)
    out
}

initialize_partable <- function(mod, ngp, ninv_items, group.equal,
                                varTable,
                                parameterization = "theta") {
  pt0 <- lavaanify_cfa(mod,
                       ngroups = ngp,
                       group.equal = group.equal,
                       parameterization = parameterization,
                       varTable = varTable)
  ninv_lds <- ninv_items[ninv_items$type == "loadings", ]
  for (j in seq_len(nrow(ninv_lds))) {
      pt0 <- remove_cons(pt0, lhs = ninv_lds$lhs[j], rhs = ninv_lds$rhs[j],
                         group = ninv_lds$group[j], op = "=~")
      pt0 <-
          remove_cons(pt0, lhs = ninv_lds$rhs[j], rhs = "",
                      group = ninv_lds$group[j], op = "~1")
  }
  ninv_ints <- ninv_items[ninv_items$type == "intercepts", ]
  for (j in seq_len(nrow(ninv_ints))) {
      pt0 <- remove_cons(pt0, lhs = ninv_ints$lhs[j], rhs = ninv_ints$rhs[j],
                         group = ninv_ints$group[j], op = "~1")
  }
  ninv_thres <- ninv_items[ninv_items$type == "thresholds", ]
  for (j in seq_len(nrow(ninv_thres))) {
      pt0 <- remove_cons(pt0, lhs = ninv_thres$lhs[j], rhs = ninv_thres$rhs[j],
                         group = ninv_thres$group[j], op = "|")
  }
  ninv_res <- ninv_items[ninv_items$type == "residuals", ]
  for (j in seq_len(nrow(ninv_res))) {
      pt0 <- remove_cons(pt0, lhs = ninv_res$lhs[j], rhs = ninv_res$rhs[j],
                         group = ninv_res$group[j], op = "~~")
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
          pt0 <- remove_cons(pt0, lhs = cov_items$lhs[k],
                             rhs = cov_items$rhs[k],
                             group = cov_items$group[k], op = "~~")
      }
  }
  pt0$id <- seq_len(nrow(pt0))
  pt0
}

fdr_alpha <- function(i, m, q = .05) {
    # Based on Benjamini & Gavrilov (2009)
    i * q / (m + 1 - i * (1 - q))
}

#' Search for invariant loadings/intercepts.
#'
#' @param config_mod Syntax of a configural invariance model to be passed to
#'   [lavaan::cfa()].
#' @param data A data frame to be passed to [lavaan::cfa()].
#' @param group Character indicating the variable name in `data` that
#'   defines the grouping variable in multiple-group CFA.
#' @param ordered Character vector indicating names of variables to be treated
#'   as binary or ordinal. IF `NULL`, all items are treated as continuous.
#' @param parameterization Character, either "theta" or "delta". "theta" should
#'   be used for invariance testing, and is the only method tested.
#' @param ... Additonal arguments passed to [lavaan::cfa()].
#' @param type Character variable indicating the stage of invariance to be
#'   searched. Currently supported options are (a) for continuous indicators,
#'   "loadings", "intercepts", "residuals", and "residual.covariances", and
#'   (b) "loadings", "thresholds", "residual.covariances", in an increasingly
#'   strict order. A stricter model (e.g., "residual.covariances") will have
#'   constraints of all previous stages.
#' @param test Character variable indicating the statistical test to be used
#'   for specification search. Currently supported options are `"mod"` for
#'   modification index using `lavaan::modindices()` and `"score"` for score
#'   test statistic using `lavaan::lavTestScore()`.
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
#' 
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' library(lavaan)
#' library(MASS)
#' # Simulate random data
#' set.seed(14)  # for reproducible results
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
#' dat_sim = simulateData(sim_mod, sample.nobs = c(100, 100, 100))
#' # Fit configural model:
#' sam_config <- paste(
#'     paste0('f =~ ', paste0("y", 1:6, collapse = " + "))
#' )
#' pinSearch(sam_config, data = dat_sim, group = "group",
#'           type = "intercepts")
#' @export
pinSearch <- function(config_mod,
                      data = NULL,
                      group = NULL,
                      ordered = NULL,
                      parameterization = "theta",
                      ...,
                      type = c("loadings", "intercepts", "thresholds",
                               "residuals", "residual.covariances"),
                      test = c("mod", "score"),
                      sig_level = .05,
                      control_fdr = FALSE,
                      effect_size = FALSE,
                      progress = FALSE) {
    type <- match.arg(type)
    test <- match.arg(test)
    base_fit <- lavaan::cfa(config_mod, group = group, data = data,
                            ordered = ordered,
                            parameterization = parameterization,
                            std.lv = TRUE, ...)
    # stored_fit <- list()
    ngp <- lavaan::lavInspect(base_fit, "ngroups")
    ninv_items <- data.frame(
        lhs = character(),
        rhs = character(),
        group = numeric(),
        type = character(),
        stringsAsFactors = FALSE
    )
    ind_names <- get_ovnames(base_fit)
    if (!is.null(ordered)) {
        types <- c("loadings", "thresholds", "residual.covariances")
        if (parameterization == "theta") {
            # add constraints on unique variances
            config_mod <- paste(c(config_mod,
                                  paste(ind_names, "~~ 1 *", ind_names)),
                                collapse = "\n")
            message("Unique variances are constrained to 1 for identification")
        }
        if (parameterization == "delta") {
            # add constraints on unique variances
            config_mod <- paste(c(config_mod,
                                  paste(ind_names, "~*~ 1 *", ind_names)),
                                collapse = "\n")
            warning("Delta parameterization has not been tested ",
                    "and likely results in untrustworthy results.")
        }
        if (!type %in% types) {
            stop("`type = ", type, "` cannot be used with ordered items")
        }
    } else {
        types <- c("loadings", "intercepts", "residuals",
                   "residual.covariances")
        if (!type %in% types) {
            stop("`type = ", type, "` cannot be used with continuous items")
        }
    }
    n_type <- which(types == type)  # number of stages
    fit_cfa <- function(mod, eq, ...) {
        lavaan::cfa(mod,
            group = group, data = data,
            ordered = ordered,
            parameterization = parameterization,
            std.lv = TRUE,
            group.equal = eq, ...
        )
    }
    fn_get_inv <- switch(test, mod = get_invmod, score = get_invscore)
    for (i in seq_len(n_type)) {
        typei <- types[i]
        if (progress) message("\n[", i, "/", n_type, "] Searching for ",
                              typei, " noninvariance\n")
        if (i == 1) {
            new_fit <- fit_cfa(config_mod, eq = "loadings", ...)
            pt0 <- lavaan::parTable(new_fit)
        } else if (i >= 2) {
            pt0 <- initialize_partable(config_mod, ngp = ngp,
                                       ninv_items = ninv_items,
                                       group.equal = types[seq_len(i)],
                                       varTable = base_fit@Data@ov,
                                       parameterization = parameterization)
            new_fit <- fit_cfa(pt0, eq = types[seq_len(i)], ...)
        }
        if (types[i] == "residual.covariances" &&
            base_fit@test$standard$df >= new_fit@test$standard$df) {
            next
        }
        lrt_base_new <- lavaan::lavTestLRT(base_fit, new_fit)
        df_diff <- lrt_base_new[2, "Df diff"]
        if ((lrt_base_new[2, "Chisq diff"] == 0 &&
             df_diff == 0) ||
            lrt_base_new[2, "Pr(>Chisq)"] >= sig_level) {
            base_fit <- new_fit
            next  # skip to next stage
        } else {
            # if (typei == "loadings") {
            #     mi_op <- "=~"
            # } else if (typei == "intercepts") {
            #     mi_op <- "~1"
            # } else if (typei %in% c("residuals", "residual.covariances")) {
            #     mi_op <- "~~"
            # }
            op <- type2op(typei)
            if (!control_fdr) {
                p_enter <- sig_level
            } else {
                num_free <- 1
                p_enter <- fdr_alpha(num_free, m = df_diff, q = sig_level)
            }
            row_to_free <- fn_get_inv(new_fit, type = typei,
                                      alpha = p_enter,
                                      ind_names = ind_names)
            if (progress) {
                total_mod <- remain_mod <- attr(row_to_free, which = "size")
                pb <- txtProgressBar(min = 0, max = total_mod, style = 3)
                pb_count <- 0
            }
            while (!is.null(row_to_free)) {
                to_free_gp <- row_to_free$group
                to_free_lhs <- row_to_free$lhs
                to_free_rhs <- row_to_free$rhs
                pt0 <- remove_cons(pt0, lhs = to_free_lhs, rhs = to_free_rhs,
                                   group = to_free_gp, op = op)
                ninv_items <- rbind(ninv_items,
                                    data.frame(
                                        lhs = to_free_lhs,
                                        rhs = to_free_rhs,
                                        group = to_free_gp,
                                        type = typei
                                    ))
                if (progress) {
                    pb_count <- max(total_mod - remain_mod + 1, pb_count + 1)
                    setTxtProgressBar(pb, pb_count)
                }
                # new_fit <- lavaan::update(new_fit, model = pt0,
                #                           group.equal = types[seq_len(i)])
                new_fit <- fit_cfa(pt0, eq = types[seq_len(i)], ...)
                if (control_fdr) {
                    num_free <- num_free + 1
                    p_enter <- fdr_alpha(num_free, m = df_diff, q = sig_level)
                }
                row_to_free <- fn_get_inv(new_fit, type = typei,
                                          alpha = p_enter,
                                          ind_names = ind_names)
                remain_mod <- attr(row_to_free, which = "size")
            }
            base_fit <- new_fit
        }
    }
    if (progress) close(pb)
    out <- list(`Partial Invariance Fit` = new_fit,
                `Non-Invariant Items` = ninv_items)
    if (effect_size) {
        out$effect_size <- es_lavaan(new_fit)
    }
    out
}
