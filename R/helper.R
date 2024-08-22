get_ovnames <- function(object) {
  object@pta$vnames$ov
}

get_lvnames <- function(object) {
    object@pta$vnames$lv[[1]]
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
    lines_added <- 0
    for (i in seq_along(loc_add)) {
        to_add <- scale_constraints(inds = inds_list[[i]], ...)
        mod_vec <- append(mod_vec, values = to_add,
                          after = loc_add[i] + lines_added)
        lines_added <- lines_added + length(to_add)
    }
    paste(mod_vec, collapse = "\n")
}
