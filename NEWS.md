# Version 0.2.1

* Dependency on package `ggparty` was removed due to maintenance issues. Visualization via `ggplot2()` is therefore not avialable any more. Plotting is still possible, however.
  
# Version 0.2.0

* First CRAN release of the package, accompanying the publication of a working
  paper on arXiv: Karapetyan S, Zeileis A, Henriksen A, Hapfelmeier A (2023).
  "Tree Models for Assessing Covariate-Dependent Method Agreement."
  arXiv 2306.04456, _arXiv.org E-Print Archive_.
  [doi:10.48550/arXiv.2306.04456](https://doi.org/10.48550/arXiv.2306.04456)


# Version 0.1.0

* Internal first development version of `coat`, a new method and R package
  for tree-based modeling of conditional method agreement. This leverages
  the `ctree()` and `mob()` functions from R package
  [partykit](https://CRAN.R-project.org/package=partykit) for the tree-based
  modeling in combination with the well-established Bland-Altman analysis for
  method agreement.
