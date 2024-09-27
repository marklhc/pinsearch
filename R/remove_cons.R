remove_cons <- function(pt, lhs, rhs, group, op, check_min2 = FALSE) {
    # Identify row with the specified parameter
    row_i <- which(pt$lhs == lhs & pt$rhs == rhs &
                       pt$op == op & pt$group == group)
    # Find other items of the same type
    eq_plabel <- unique(unlist(pt[pt$op == "==", c("lhs", "rhs")]))
    row_labelled <- which(pt$op == op & pt$plabel %in% eq_plabel)
    ind_labelled <- unique(pt[row_labelled, op2col(op)])
    if (length(ind_labelled) <= 2 && check_min2) {
        return(pt)
    }
    mi_plab <- pt$plabel[row_i]
    pt$label[row_i] <- ""
    if (length(mi_plab) == 0)
        return(pt)
    ind_plabs <-
        pt$plabel[pt$op == op & pt$lhs == lhs & pt$rhs == rhs]
    pt_eq <- pt[pt$op == "==" &
                        pt$lhs %in% ind_plabs &
                        pt$rhs %in% ind_plabs,]
    # If it appears on rhs, remove it
    pt_new <-
        pt[!pt$id %in% pt_eq[pt_eq$rhs == mi_plab, "id"],]
    # Replace lhs with the next group (make sure it's not on the table)
    if (mi_plab %in% pt_eq$lhs) {
        pt_new[pt_new$op == "==" &
                   pt_new$lhs == mi_plab, "lhs"] <- ind_plabs[group + 1]
        pt_new <-
            pt_new[!(pt_new$op == "==" & pt_new$lhs == pt_new$rhs),]
    }
    pt_new
}

get_invmod <- function(object, type, alpha = .05, pt0, ...) {
    if (missing(pt0)) pt0 <- lavaan::parTable(object)
    op <- type2op(type)
    pt0_eq <- pt0[pt0$label != "" & pt0$free >= 0, ]
    suppressWarnings({
        mis <- lavaan::modindices(
            object,
            op = op,
            minimum.value = stats::qchisq(alpha, 1, lower.tail = FALSE),
            sort. = TRUE,
            free.remove = FALSE
        )
    })
    # TODO: See if ind_names is needed
    # Only parameters in the original model
    # if (type == "loadings") {
    #     mis_ids <- rownames(mis)[mis$rhs %in% ind_names &
    #                                  rownames(mis) %in% pt0_eq$id]
    # } else if (type %in% c("intercepts", "thresholds")) {
    #     mis_ids <- rownames(mis)[mis$lhs %in% ind_names &
    #                                  rownames(mis) %in% pt0_eq$id]
    # } else if (type == "residuals") {
    #     mis_ids <- rownames(mis)[mis$lhs %in% ind_names &
    #                                  rownames(mis) %in% pt0_eq$id &
    #                                  mis$lhs == mis$rhs]
    # } else if (type == "residual.covariances") {
    #     mis_ids <- rownames(mis)[mis$lhs %in% ind_names &
    #                                  mis$rhs %in% ind_names &
    #                                  rownames(mis) %in% pt0_eq$id &
    #                                  mis$lhs != mis$rhs]
    # }
    mis_ids <- rownames(mis)[rownames(mis) %in% pt0_eq$id]
    if (length(mis_ids) == 0) return(NULL)
    pt_mis <- pt0_eq[match(mis_ids, pt0_eq$id), ]
    # Only parameters that are free
    pt_mis <- pt_mis[pt_mis$free >= 0, ]
    out <- mis[rownames(pt_mis)[1], ]
    attr(out, "size") <- nrow(pt_mis)
    out
}

get_invscore <- function(object, type, alpha = .05, pt0, ...) {
    if (missing(pt0)) pt0 <- lavaan::parTable(object)
    op <- type2op(type)
    pt0_eq <- pt0[pt0$label != "" & pt0$free >= 0 & pt0$op == op, ]
    # Find relevant constraints
    pt0_cons <- pt0[pt0$op == "==", ]
    cons_to_test <- which(pt0_cons$lhs %in% pt0_eq$plabel |
        pt0_cons$rhs %in% pt0_eq$plabel)
    mis <- lavaan::lavTestScore(object, release = cons_to_test)
    cons_to_free <- mis$uni[which.max(mis$uni$X2), ]
    if (cons_to_free$p.value > alpha) return(NULL)
    out <- pt0_eq[which(pt0_eq$plabel == cons_to_free$rhs), ]
    attr(out, "size") <- sum(mis$uni$p.value < alpha)
    out
}

get_invlrt <- function(object, type, alpha = .05, pt0, ...) {
    if (missing(pt0)) pt0 <- lavaan::parTable(object)
    op <- type2op(type)
    pt0_eq <- pt0[pt0$label != "" & pt0$free >= 0 & pt0$op == op, ]
    # Find relevant constraints
    # pt0_cons <- pt0[pt0$op == "==", ]
    to_extract <- c("Chisq diff", "Df diff", "Pr(>Chisq)")
    if (FALSE) {
        # by item
        par_sets <- split(pt0_eq, f = pt0_eq[c("lhs", "op", "rhs")])
        lrt_mat <- matrix(
            nrow = length(par_sets), ncol = length(to_extract),
            dimnames = list(NULL, to_extract)
        )
        for (j in seq_along(par_sets)) {
            set <- par_sets[[j]]
            pt_new <- pt0
            for (i in seq_len(nrow(set))) {
                pt_new <- remove_cons(
                    pt_new, set$lhs[i], set$rhs[i],
                    set$group[i], op
                )
            }
            if (nrow(pt0) == nrow(pt_new)) {
                lrt_mat[j, ] <- c(0, 0, NA)
            } else {
                lrt_j <- lavaan::lavTestLRT(
                    object,
                    lavaan::cfa(pt_new, ...)
                )
                lrt_mat[j, ] <- as.numeric(lrt_j[2, to_extract])
            }
        }
        if (min(lrt_mat[, "Pr(>Chisq)"]) > alpha) return(NULL)
        pt0_eq <- par_sets[[which.min(lrt_mat[, "Pr(>Chisq)"])]]
    }
    lrt_mat <- matrix(
        nrow = nrow(pt0_eq), ncol = length(to_extract),
        dimnames = list(NULL, to_extract)
    )
    for (i in seq_len(nrow(pt0_eq))) {
        pt_new <- remove_cons(
            pt0, pt0_eq$lhs[i], pt0_eq$rhs[i], pt0_eq$group[i], op
        )
        if (nrow(pt0) == nrow(pt_new)) {
            lrt_mat[i, ] <- c(0, 0, NA)
        } else {
            lrt_i <- lavaan::lavTestLRT(
                object,
                lavaan::cfa(pt_new, ...)
            )
            lrt_mat[i, ] <- as.numeric(lrt_i[2, to_extract])
        }
    }
    mis <- cbind(pt0_eq, lrt_mat)

    out <- mis[which.min(mis$"Pr(>Chisq)"), ]
    if (nrow(out) == 0 || out$"Pr(>Chisq)" > alpha) return(NULL)
    attr(out, "size") <- sum(mis$"Pr(>Chisq)" < alpha)
    out
}

# TODO:
# - Separate functions for (a) extracting a list of statistics and
#   (b) extracting the smallest one
# - Clean up get_invscore()
# - Make functions consistent