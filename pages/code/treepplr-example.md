# treepplr reform


- [Diversification models (CRBD)](#diversification-models-crbd)
- [Tree inference (pruning)](#tree-inference-pruning)
- [Host-repertoire (subroot)](#host-repertoire-subroot)
- [Qualitative trait evolution](#qualitative-trait-evolution)

I’ve been reworking treepplr to make it faster, easier to use, and less
redundant. Along the way I reorganized the code so it’s easier to
maintain, and fixed a few bugs. This vignette shows how to convert R
objects (lists, data frames, phylo objects, etc.) into TreePPL input
files, and how to convert TreePPL output back into R. I’ll explain the
new and updated functions as they come up.

The core workflow in treepplr is:

**R objects** → `tp_data()` → `tp_sampler()` & `tp_run()` → **R
objects**

This vignette covers:

1.  [Diversification models (CRBD)](#diversification-models-crbd)

2.  [Tree inference (pruning)](#tree-inference-pruning)

3.  [Host-repertoire (subroot)](#host-repertoire-subroot)

4.  [Qualitative trait evolution](#qualitative-trait-evolution)

> Make sure you’re using treepplr from my branch `improvements`.

``` r
library(devtools)
library(ape)

load_all()
```

## Diversification models (CRBD)

## Tree inference (pruning)

## Host-repertoire (subroot)

## Qualitative trait evolution
