# aDIFic

*Working name ("automated DIF with item covariates"); the final package name is still to be decided.*

Anchor-free estimation and inference for differential item functioning (DIF) in moderated nonlinear factor analysis (MNLFA) with item covariates.

- **Model:** item difficulty DIF with respect to any number of continuous or categorical person covariates, decomposed into a part explained by item covariates (e.g., word class, age of acquisition) and item-specific residual DIF.
- **Identification:** no anchor items. The model is identified by a sparsity assumption on residual DIF, generalizing the minimal-L1 condition of Chen et al. (2023).
- **Inference:** a threshold-centered parametric bootstrap gives intervals and tests for impact, item-covariate effects, and item-specific DIF.

**Status:** under development (version 0.0.0.9000), not yet ready for use. The interface is not final.

## Licence

GPL (>= 3).
