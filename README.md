# MNLFAic

Anchor-free estimation and inference for differential item functioning (DIF) in moderated nonlinear factor analysis (MNLFA) with item covariates.

- **Model:** a two-parameter logistic MNLFA in which any number of continuous or categorical person covariates affect the latent mean, the latent variance and every item's difficulty. Each item's DIF is split into a part explained by item covariates (e.g., a word's part of speech or age of acquisition) and item-specific residual DIF.
- **Identification:** no anchor items. The model is identified by assuming that residual DIF is sparse, which generalizes the minimal-L1 condition of Chen et al. (2023). Without item covariates, the analysis is Chen et al.'s method for an arbitrary number of continuous or categorical covariates.
- **Inference:** a threshold-centered parametric bootstrap: the model estimates are redrawn from their estimated sampling distribution and the identification step is repeated on each draw, with the draws centered at estimates in which residual DIF smaller than a threshold (in standard errors) is set to zero. It gives intervals and tests for impact, item-covariate effects and item-specific DIF, with Benjamini-Hochberg flags for residual DIF and diagnostics for the identifying assumption.

**Status:** under development (version 0.0.0.9000). The interface may still change.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sanfordstudent/MNLFAic", build_vignettes = TRUE)
```

## Example

```r
library(MNLFAic)

# Simulated vocabulary test: 1,000 children (gender, centered age), 30 words
# (age of acquisition, part of speech). Verbs are harder for girls, later-acquired
# words get easier with age, and word 4 has item-specific gender DIF.
set.seed(1)
persons <- data.frame(female = rbinom(1000, 1, 0.5), age = runif(1000, -1.5, 1.5))
words <- data.frame(aoa = rnorm(30, 5, 1.5), pos = rep(c("noun", "verb", "adjective"), 10))
effects <- rbind(aoa = c(0, -0.25), pos_adjective = c(0, 0), pos_verb = c(0.3, 0))
residual <- matrix(0, 30, 2); residual[4, 1] <- 0.6
sim <- simulate_mnlfa_ic(persons, ~ female + age, words, ~ aoa + pos,
                         slopes = runif(30, 0.8, 2), intercepts = rnorm(30), impact = c(0.2, 0.5),
                         effects = effects, residual_dif = residual,
                         item_omit = list(pos = "noun"), seed = 1)

# responses: N x J matrix of 0/1 (NA allowed); one row of person_data per person,
# one row of item_data per item
fit <- mnlfa_ic(sim$responses, person_data = persons, person = ~ female + age,
                item_data = words, item = ~ aoa + pos,
                control = mnlfa_control(item_omit = list(pos = "noun"), seed = 1))

fit                                 # impact and the number of flagged items
summary(fit)                        # impact, item-covariate effects, flagged residual DIF
dif_table(fit, flagged_only = TRUE) # item-level results
diagnose(fit)                       # checks of the identifying assumption
```

See `vignette("MNLFAic")` for a worked example on simulated data, and `vignette("method")` for the model and the procedure.

## Citation

`citation("MNLFAic")`. The method without item covariates is from:

Chen, Y., Li, C., Ouyang, J., & Xu, G. (2023). DIF statistical inference without knowing anchoring items. *Psychometrika, 88*, 1097-1122. <https://doi.org/10.1007/s11336-023-09930-9>

## License

GPL (>= 3).
