# mixedcirc <img src="/assets/hexSticker.png" width="150" align="right" />

## Overview 

`mixedcirc` uses mixed models to perform differential circadian rhythm analysis.

## Installation

Install requirements

```r
BiocManager::install( c("roxygen2","multcomp", "doFuture", "foreach", "future", "future.apply", "funtimes", "lme4", "lmerTest", "limma", "ggplot2",
                        "ggsci", "ggpubr", "multtest", "mixOmics", "cowplot", "dplyr","emmeans","MuMIn","GenomicRanges","genomation","performance",
                        "rmcorr","RhpcBLASctl","aod"))
```

```r
devtools::install_github("mixedcirc/mixedcirc")
library(mixedcirc)
```

To install the vignettes, use:

```r
devtools::install_github("mixedcirc/mixedcirc", build_vignettes = TRUE)
library(mixedcirc)
```

## License

[Apache License Version 2.0](LICENSE.md)

## Maintainer

[Payam Emami](https://github.com/PayamEmami)



