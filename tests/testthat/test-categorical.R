test_that("pinSearch() works properly for binary data", {
    ps1 <- pinSearch(" f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ",
        data = df1, group = "group", type = "thresholds",
        ordered = paste0("yy", 1:7)
    )
    expect_identical(ps1[[2]]$lhs, c("yy7", "yy2"))
})

test_that("pinSearch() works properly for 1 noninvariant item", {
    ps1b <- cfa(" f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ",
        data = df1, group = "group",
        group.equal = c("loadings", "thresholds"),
        group.partial = c("yy7|t1"),
        ordered = TRUE,
    )
    expect_no_warning(pin_effsize(ps1b))
})

test_that("pinSearch() works properly for noninvariant uniqueness", {
    # Unique variances should shift to loadings and intercepts
    ps2 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                     data = df2, group = "group", type = "thresholds",
                     ordered = paste0("yy", 1:7))
    expect_true("yy7" %in% ps2[[2]]$lhs | "yy7" %in% ps2[[2]]$rhs)
})

test_that("pinSearch() works properly for noninvariant unique covariances", {
    ps3 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                     data = df3, group = "group", type = "residual.covariances",
                     ordered = paste0("yy", 1:7))
    expect_setequal(ps3[[2]]$lhs, c("yy2", "yy7"))
    ps4 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7
                       yy2 ~~ yy3
                       yy4 ~~ yy5 ',
                     data = df3, group = "group", type = "residual.covariances",
                     ordered = paste0("yy", 1:7))
    expect_identical(ps4[[2]][3, "lhs"], "yy2")
    expect_identical(ps4[[2]][3, "rhs"], "yy3")
})

test_that("pinSearch() gives error when type = 'residuals' for ordered items", {
    expect_error(
        pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
              data = df1, group = "group", type = "residuals",
              ordered = paste0("yy", 1:7))
    )
})

test_that("type = 'thresholds' gives error for continuous items", {
    expect_error(
        pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                  data = df1, group = "group", type = "thresholds")
    )
})

test_that("Works for three groups and ordinal items", {
    ps5 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                     data = dfo, group = "group", type = "thresholds",
                     ordered = paste0("yy", 1:7))
    expect_equal(sort(with(ps5[[2]], paste(lhs, rhs))),
                 c("yy2 t1", "yy4 t3", "yy5 t3", "yy7 t1"))
    ps5_re <- cfa(' f =~ NA * yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7
                    f ~~ c(0.5, NA, NA) * f
                    f ~ c(1, NA, NA) * 1
                    yy1 ~~ 1 * yy1
                    yy2 ~~ 1 * yy2
                    yy3 ~~ 1 * yy3
                    yy4 ~~ 1 * yy4
                    yy5 ~~ 1 * yy5
                    yy6 ~~ 1 * yy6
                    yy7 ~~ 1 * yy7
                    yy2 | c(t2, t2, t23) * t1
                    yy4 | c(t4, t4, t43) * t3
                    yy5 | c(t5, t5, t53) * t3
                    yy7 | c(t7, t7, t73) * t1 ',
                  data = dfo, group = "group", ordered = TRUE,
                  group.equal = c("loadings", "thresholds", "residuals"),
                  group.partial = c("yy7|t1", "yy2|t1", "yy4|t3", "yy5|t3"),
                  parameterization = "theta")
    fmacs_ps5 <- pin_effsize(ps5[[1]])
    expect_equal(c(fmacs_ps5),
                 c(.1227653, .07881984, .06000593, .07972769),
                 tolerance = 0.00001)
    expect_equal(fmacs_ps5,
                 pin_effsize(ps5_re),
                 tolerance = 0.00001)
})
