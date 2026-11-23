# treepplr reform


Source:
[treepplr-example.qmd](https://github.com/foersterst/foersterst.github.io/blob/master/pages/code/treepplr-example.qmd)

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

CRBD, ClaDS2, and BDD2 models are all fitted the same way, so I’ll use
CRBD as the example here. You can check the correct spelling of each
model name with `tp_model_library()`, which has a new look:

``` r
tp_model_library()
```

    ── diversification ──

    • bdd2

    • bdd2_delayed

    • clads2

    • clads2_delayed

    • crbd

    • crbd_analytical

    • crbd_delayed

    ── host-repertoire-evolution ──

      • flat-root-prior-HRM

      • subroot-HRM

    ── lang ──

      • anonymous

      • blub

      • coin

      • distributions

      • externals

      • hello

      • matrix-tests

      • mini-mat-test

      • poly

      • ser-list

      • tensors

    ── metabarcoding ──

      • optmodel

    ── pheno-mol ──

      • qt

    ── tree-inference ──

      • substmodel_belief_propagation

      • tree_inference

      • tree_inference_pruning

      • tree_inference_pruning_gtr

      • tree_inference_pruning_gtr_gamma

      • tree_inference_pruning_scaled

For this example, I’m gonna use a toy tree with four tips and height =
10:

``` r
phy <- read.tree(text = "((A:4,B:4):6,(C:7,D:7):3);")

plot.phylo(phy, direction = "upwards")
```

![](treepplr-example_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
# tree height
max(node.depth.edgelength(phy))
```

    [1] 10

``` r
class(phy)
```

    [1] "phylo"

\

## Tree inference (pruning)

## Host-repertoire (subroot)

## Qualitative trait evolution
