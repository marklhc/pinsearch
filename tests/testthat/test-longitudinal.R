library(lavaan)
library(semTools)

test_that(
    "Same configural model as semTools",
    code = {
        mod_cat <- "
    FU1 =~ u1 + u2 + u3 + u4
    FU2 =~ u5 + u6 + u7 + u8
"
        ## the 2 factors are actually the same factor (FU) measured twice
        long_fac <- list(FU = c("FU1", "FU2"))

        ## CONFIGURAL model: no constraints across groups or repeated measures
        mod_config <- measEq.syntax(
            configural.model = mod_cat,
            data = datCat,
            ordered = paste0("u", 1:8),
            parameterization = "theta",
            ID.fac = "std.lv",
            ID.cat = "Wu.Estabrook.2016",
            longFacNames = long_fac,
            longIndNames = list(
                i1 = c("u1", "u5"),
                i2 = c("u2", "u6"),
                i3 = c("u3", "u7"),
                i4 = c("u4", "u8")
            ),
            long.equal = "thresholds",
            long.partial = "i1|t1",
            auto = TRUE
        )

        fit1 <- cfa(
            as.character(mod_config),
            data = datCat,
            ordered = paste0("u", 1:8),
            parameterization = "theta"
        )

        fit2 <- longcfa(
            matrix(
                c(
                    "u1", "u2", "u3", "u4",
                    "u5", "u6", "u7", "u8"
                ),
                ncol = 2
            ),
            lv_names = c("FU1", "FU2"),
            lag_cov = TRUE,
            data = datCat, ordered = TRUE,
            parameterization = "theta"
        )
        expect_equal(fitmeasures(fit1),
                     fitmeasures(fit2),
                     ignore_attr = TRUE)
    }
)

# Test missing items

# Test conflicting syntax

