lavaanify_cfa <- function(config_mod, ngroups, group.equal) {
    lavaan::lavaanify(config_mod, meanstructure = TRUE, ngroups = ngroups,
                      group.equal = group.equal,
                      std.lv = TRUE, int.ov.free = TRUE, int.lv.free = FALSE,
                      auto.var = TRUE, auto.cov.lv.x = TRUE,
                      auto.cov.y = TRUE, auto.th = TRUE, auto.delta = TRUE,
                      auto.efa = TRUE)
}
op2col <- function(op) {
    switch(op,
           `=~` = "rhs",
           `~1` = "lhs",
           `~~` = "lhs")
}
remove_cons <- function(object, item, group, op) {
    # Remove constraints
    col <- op2col(op)
    matched <- if (op == "~~")
        object$lhs == object$rhs
    else
        TRUE
    # Identify row with the specified parameter
    row_i <- which(object[[col]] == item & matched &
                       object$op == op & object$group == group)
    mi_plab <- object$plabel[row_i]
    object$label[row_i] <- ""
    if (length(mi_plab) == 0)
        return(object)
    ind_plabs <-
        object$plabel[object$op == op & object[[col]] == item &
                          matched]
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
get_invmod <- function(object, op, mod_min, ind_names, pt0) {
    if (missing(pt0)) pt0 <- lavaan::parTable(object)
    col <- op2col(op)
    pt0_eq <- pt0[pt0$label != "" & pt0$free >= 0, ]
    mis <- lavaan::modindices(
        object,
        op = op,
        minimum.value = mod_min,
        sort. = TRUE,
        free.remove = FALSE
    )
    # Only parameters in the original model
    if (op == "~~") {
        mis_ids <- rownames(mis)[mis[[col]] %in% ind_names &
                                     rownames(mis) %in% pt0_eq$id &
                                     mis$lhs == mis$rhs]
    } else {
        mis_ids <- rownames(mis)[mis[[col]] %in% ind_names &
                                     rownames(mis) %in% pt0_eq$id]
    }
    if (length(mis_ids) == 0) return(NULL)
    pt_mis <- pt0_eq[match(mis_ids, pt0_eq$id), ]
    # Only parameters that are free
    pt_mis <- pt_mis[pt_mis$free >= 0, ]
    mis[rownames(pt_mis)[1], ]
}
initialize_partable <- function(mod, ngp, ninv_items, group.equal) {
  pt0 <- lavaanify_cfa(mod,
                       ngroups = ngp,
                       group.equal = group.equal)
  ninv_lds <- ninv_items[ninv_items$type == "loadings", ]
  for (j in seq_len(nrow(ninv_lds))) {
      pt0 <- remove_cons(pt0, ninv_lds$item[j], ninv_lds$group[j], "=~")
      pt0 <-
          remove_cons(pt0, ninv_lds$item[j], ninv_lds$group[j], "~1")
  }
  ninv_ints <- ninv_items[ninv_items$type == "intercepts", ]
  for (j in seq_len(nrow(ninv_ints))) {
      pt0 <- remove_cons(pt0, ninv_ints$item[j], ninv_ints$group[j], "~1")
  }
  pt0$id <- seq_len(nrow(pt0))
  pt0
}

#' Search for invariant loadings/intercepts.
#'
#' @param config_mod Syntax of a configural invariance model to be passed to
#'   \code{\link[lavaan]{cfa}}.
#' @param data A data frame to be passed to \code{\link[lavaan]{cfa}}.
#' @param group Character indicating the variable name in \code{data} that
#'   defines the grouping variable in multiple-group CFA.
#' @param ... Additonal arguments passed to \code{\link[lavaan]{cfa}}.
#' @param type Character variable indicating the stage of invariance to be
#'   searched.
#' @param sig_level Significance level used to determine whether the parameter
#'   associated with the highest modification index should be removed. Default
#'   is .05.
#'
#' @return The sum of \code{x} and \code{y}.
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
                      ...,
                      type = c("loadings", "intercepts", "residuals"),
                      sig_level = .05) {
    type <- match.arg(type)
    base_fit <- lavaan::cfa(config_mod, group = group, data = data,
                            std.lv = TRUE, ...)
    # stored_fit <- list()
    ngp <- lavaan::lavInspect(base_fit, "ngroups")
    ninv_items <- data.frame(
        items = character(),
        group = numeric(),
        type = character(),
        stringsAsFactors = FALSE
    )
    chisq_cv <- stats::qchisq(sig_level, 1, lower.tail = FALSE)
    ind_names <- get_ovnames(base_fit)
    types <- c("loadings", "intercepts", "residuals")
    n_type <- which(types == type)  # number of stages
    for (i in seq_len(n_type)) {
        typei <- types[i]
        if (i == 1) {
            new_fit <- lavaan::cfa(config_mod, group = group, data = data,
                                   std.lv = TRUE,
                                   group.equal = "loadings", ...)
            pt0 <- lavaan::parTable(new_fit)
        } else if (i >= 2) {
            pt0 <- initialize_partable(config_mod, ngp = ngp,
                                       ninv_items = ninv_items,
                                       group.equal = types[seq_len(i)])
            new_fit <- lavaan::cfa(
                pt0,
                group = group,
                data = data,
                group.equal = types[seq_len(i)],
                std.lv = TRUE,
                ...
            )
        }
        if (lavaan::lavTestLRT(base_fit, new_fit)[2, "Pr(>Chisq)"] >= sig_level) {
            base_fit <- new_fit
            next  # skip to next stage
        } else {
            if (typei == "loadings") {
                mi_op <- "=~"
                mi_col <- "rhs"
            } else if (typei == "intercepts") {
                mi_op <- "~1"
                mi_col <- "lhs"
            } else if (typei == "residuals") {
                mi_op <- "~~"
                mi_col <- "lhs"
            }
            largest_mi_row <- get_invmod(new_fit, op = mi_op,
                                         mod_min = chisq_cv,
                                         ind_names = ind_names)
            while (!is.null(largest_mi_row)) {
                mi_gp <- largest_mi_row$group
                mi_ind <- largest_mi_row[[mi_col]]
                pt0 <- remove_cons(pt0, item = mi_ind, group = mi_gp,
                                   op = mi_op)
                ninv_items <- rbind(ninv_items,
                                    data.frame(
                                        item = mi_ind,
                                        group = mi_gp,
                                        type = typei
                                    ))
                # new_fit <- lavaan::update(new_fit, model = pt0,
                #                           group.equal = types[seq_len(i)])
                new_fit <- lavaan::cfa(
                    pt0,
                    group = group,
                    data = data,
                    group.equal = types[seq_len(i)],
                    std.lv = TRUE,
                    ...
                )
                largest_mi_row <- get_invmod(new_fit, op = mi_op,
                                             mod_min = chisq_cv,
                                             ind_names = ind_names)
            }
            base_fit <- new_fit
        }
    }
    list(`Partial Invariance Fit` = new_fit,
         `Non-Invariant Items` = ninv_items)
}
