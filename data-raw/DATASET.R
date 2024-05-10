## code to prepare `lui_sim` dataset goes here
library(haven)
library(dplyr)
orig_dat <- haven::read_sav("https://osf.io/download/wxjsg/")
dat <- orig_dat |>
    filter(eth %in% c(1, 2, 4), gender %in% 1:2) |>
    mutate( # across(class1:class15, .fns = as.numeric),
        across(c(eth, gender), .fns = as_factor),
        group = interaction(eth, gender),
        # Make binary
        across(c(class1:class15),
            .fns = \(x) {
                factor(x >= 3, labels = c("disagree/neutral", "agree"))
            }
        )
    ) |>
    arrange(group)
config_mod <- "
f1 =~ class1 + class2 + class3 + class4 + class5 + class6 + 
      class7 + class8 + class9 + class10 + class11 + class12 + 
      class13 + class14 + class15
"
config_fit <- cfa(c(config_mod, "class1 ~~ class2"),
    data = dat, group = "group",
    ordered = TRUE
)
# Simulate data
lui_sim <- simulateData(
    parTable(config_fit),
    sample.nobs = c(160, 80, 70, 350, 140, 110)
)

usethis::use_data(lui_sim, overwrite = TRUE)
