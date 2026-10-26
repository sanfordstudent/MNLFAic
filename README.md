# MNLFAic

Anchor-free estimation and inference for differential item functioning (DIF) in moderated nonlinear factor analysis with item covariates (MNLFA-IC).

- **Model:** a two-parameter logistic MNLFA in which any number of continuous or categorical person covariates affect the latent mean, the latent variance and every item's difficulty. Each item's DIF is split into a part explained by item covariates (e.g., a word's part of speech or age of acquisition) and item-specific residual DIF.
- **Identification:** no anchor items. The model is identified by assuming that residual DIF is sparse, which generalizes the minimal-L1 condition of Chen et al. (2023). Without item covariates, the analysis is Chen et al.'s method extended to an arbitrary number of continuous or categorical covariates.
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

# Simulated mathematics test: 1,000 students (English-learner status, grade 3-5
# centered at 4), 30 items (reading load, content area). Reading-heavy items are
# harder for English learners, geometry items get easier with grade, the spread
# of ability grows with grade, and item 4 has item-specific English-learner DIF.
set.seed(1)
persons <- data.frame(el = rbinom(1000, 1, 0.3), grade = sample(3:5, 1000, replace = TRUE) - 4)
items <- data.frame(reading = round(runif(30, 10, 80)),
                    domain = rep(c("number", "geometry", "data"), 10))
effects <- rbind(reading = c(0.3, 0), domain_data = c(0, 0), domain_geometry = c(0, -0.25))
residual <- matrix(0, 30, 2); residual[4, 1] <- 0.6
sim <- simulate_mnlfa_ic(persons, ~ el + grade, items, ~ reading + domain,
                         slopes = runif(30, 0.8, 2), intercepts = rnorm(30),
                         impact = c(-0.4, 0.5), variance_impact = c(0, 0.2),
                         effects = effects, residual_dif = residual,
                         item_omit = list(domain = "number"), seed = 1)

# responses: N x J matrix of 0/1 (NA allowed); one row of person_data per person,
# one row of item_data per item
fit <- mnlfa_ic(sim$responses, person_data = persons, person = ~ el + grade,
                item_data = items, item = ~ reading + domain,
                control = mnlfa_control(item_omit = list(domain = "number"), seed = 1))

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
