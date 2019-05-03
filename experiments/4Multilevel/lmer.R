install.packages("brms")

library(brms)
library(lme4)
library(dplyr)

lmm <- lmer(Compensation ~ (1 | Company),
            data=directors, REML=FALSE)
summary(lmm)

brmm <- brm(Compensation ~ Age + Male + (1 | Company),
            data=directors)
