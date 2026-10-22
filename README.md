# aDIFic

*Working name ("automated DIF with item covariates"); the final package name is still to be decided.*

Anchor-free estimation and inference for differential item functioning (DIF) in moderated nonlinear factor analysis (MNLFA) with item covariates.

- **Model:** a two-parameter logistic MNLFA in which any number of continuous or categorical person covariates affect the latent mean, the latent variance and every item's difficulty. Each item's DIF is split into a part explained by item covariates (e.g., a word's part of speech or age of acquisition) and item-specific residual DIF.
- **Identification:** no anchor items. The model is identified by assuming that residual DIF is sparse, which generalizes the minimal-L1 condition of Chen et al. (2023). Without item covariates, the analysis is Chen et al.'s method for several person covariates at once.
- **Inference:** a threshold-centered parametric bootstrap gives intervals and tests for impact, item-covariate effects and item-specific DIF, with Benjamini-Hochberg flags for residual DIF and diagnostics for the identifying assumption.

**Status:** under development (version 0.0.0.9000). The interface may still change.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sanfordstudent/aDIFic", build_vignettes = TRUE)
```

## Example

```r
library(aDIFic)

# responses: N x J matrix of 0/1 (NA allowed)
# persons:   data frame with one row per person (female, age centered)
# words:     data frame with one row per item (aoa, pos)
fit <- mnlfa_ic(responses, person_data = persons, person = ~ female + age,
                item_data = words, item = ~ aoa + pos,
                control = mnlfa_control(item_omit = list(pos = "noun"), seed = 1))

fit                                 # impact and the number of flagged items
summary(fit)                        # impact, item-covariate effects, flagged residual DIF
dif_table(fit, flagged_only = TRUE) # item-level results
diagnose(fit)                       # checks of the identifying assumption
```

See `vignette("aDIFic")` for a worked example on simulated data, and `vignette("method")` for the model and the procedure.

## Citation

`citation("aDIFic")`. The method without item covariates is from:

Chen, Y., Li, C., Ouyang, J., & Xu, G. (2023). DIF statistical inference without knowing anchoring items. *Psychometrika, 88*, 1097-1122.

## Licence

GPL (>= 3).
