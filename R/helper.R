get_ovnames <- function(object) {
    object@pta$vnames$ov
}

get_lvnames <- function(object) {
    object@pta$vnames$lv
}

scale_constraints <- function(inds,
                              parameterization = c("theta", "delta")) {
    parameterization <- match.arg(parameterization)
    if (parameterization == "theta") {
        paste(inds, "~~ 1 *", inds)
    } else if (parameterization == "delta") {
        paste(inds, "~*~ 1 *", inds)
    }
}

add_scale_constraints_group <- function(mod, inds_list, ...) {
    mod_vec <- strsplit(mod, split = "\n")[[1]]
    loc_group <- grep("group:", mod_vec)
    loc_add <- c(loc_group[-1] - 1, length(mod_vec))
    for (i in rev(seq_along(loc_add))) {
        to_add <- scale_constraints(inds = inds_list[[i]], ...)
        mod_vec <- append(mod_vec, values = to_add,
                          after = loc_add[i])
    }
    paste(mod_vec, collapse = "\n")
}

latent_var_syntax <- function(lv_names, val) {
    paste(
        c(
            "# Latent variances",
            paste(lv_names, "~~", val, "*", lv_names)
        ),
        collapse = "\n"
    )
}

latent_mean_syntax <- function(lv_names, val) {
    paste(
        c(
            "# Latent means",
            paste(lv_names, "~", val, "* 1")
        ),
        collapse = "\n"
    )
}

add_iden_constraints_group <- function(
    mod, lvs,
    type = c("config", "loadings", "intercepts", "thresholds",
             "residuals", "residual.covariances")) {
    type <- match.arg(type)
    mod_vec <- strsplit(mod, split = "\n")[[1]]
    loc_group <- grep("group:", mod_vec)
    loc_add <- c(loc_group[-1] - 1, length(mod_vec))
    num_groups <- length(loc_group)
    if (type == "config") {
        lat_var <- mapply(latent_var_syntax, lv_names = lvs,
                          val = rep(1, num_groups))
        lat_mean <- mapply(latent_mean_syntax, lv_names = lvs,
                           val = rep(0, num_groups))
    } else if (type == "loadings") {
        lat_var <- mapply(latent_var_syntax, lv_names = lvs,
                          val = c(1, rep(NA, num_groups - 1)))
        lat_mean <- mapply(latent_mean_syntax, lv_names = lvs,
                           val = rep(0, num_groups))
    } else {
        lat_var <- mapply(latent_var_syntax, lv_names = lvs,
                          val = c(1, rep(NA, num_groups - 1)))
        lat_mean <- mapply(latent_mean_syntax, lv_names = lvs,
                           val = c(0, rep(NA, num_groups - 1)))
    }
    for (i in rev(seq_along(loc_add))) {
        to_add <- paste0(lat_var[[i]], "\n", lat_mean[[i]])
        mod_vec <- append(mod_vec,
            values = to_add,
            after = loc_add[i]
        )
    }
    paste(mod_vec, collapse = "\n")
}
