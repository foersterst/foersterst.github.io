

---
title: "treepplr reform"
---

* TOC
{:toc}

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

![](treepplr-example_files/figure-commonmark/unnamed-chunk-4-1.png)

``` r
# tree height
max(node.depth.edgelength(phy))
```

    [1] 10

``` r
class(phy)
```

    [1] "phylo"

We need to create a named list with the input items required by the CRBD
model. Each element of the list must be named exactly like an argument
of the model; the order doesn’t matter. To see the exact names, we use
the new function `tp_expected_input()`:

``` r
tp_expected_input(tppl = "crbd")
```

    ── Expected Input ──

    Model file:
    '/tmp/treeppl-0.4/6yjjihgyj569rg7bw7jkqxy3z1-treeppl-0.4/lib/mcore/treeppl/models/diversification/crbd.tppl'

    Model name: crbd

    Expected arguments:

    • tree (tpPhylo)

    • rho (numeric)

The CRBD expects two arguments: `tree` and `rho`, and these are the
names we have provide to the items in our list. `rho` is a numeric
object and `tree` is a tpPhylo object that is specific to treepplr. We
create the tpPhylo object from an `ape::phylo` object using
`tp_phylo_to_tpPhylo()`:

``` r
phy_2 <- tp_phylo_to_tpPhylo(phy = phy)

# tpPhylo objects have a printing method similar to ape::phylo:
phy_2
```

    Phylogenetic tree of class <tpPhylo>, with 4 tips and 3 internal nodes.

    Rooted; includes branch lengths; root age = 10

> We’ll see later that some models require trees with extra information
> at the leaves and internal nodes. Function `tp_phylo_to_tpPhylo()`
> allows such trees to be created by using the arguments `leaf_extra`
> and `node_extra`. If needed, it is also possible to change the labels
> for node and leaf constructors using the arguments `node_ctor` and
> `leaf_ctor`, respectively.

Now, we create our named list with the data required by our model:

``` r
L <- list(tree = phy_2, rho = 0.35)
```

Next, we call `tp_data()` to assemble the JSON file required by TreePPL
from the R objects we provided:

``` r
data_crbd <- tp_data(model = "crbd", data = L, check = TRUE, print = FALSE)
```

    ── Input Check ──

    JSON file:
    '/private/var/folders/vh/0rc0zn655lsd4yxfvkgdk1hr0000gp/T/RtmpoJNGgJ//input-crbd-4d708a5.json'

    Model name: crbd

    ✔ All 2 arguments found: `tree` and `rho`.

`tp_data()` has two new arguments, `check` and `print`. When
`check = TRUE`, the function calls the new function `tp_input_check()`,
which reads the `model function` signature of the specified model and
compares the argument names with the top-level field names in the JSON
input file. It reports missing and extra fields. It only checks names,
not types.

\

## Tree inference (pruning)

## Host-repertoire (subroot)

## Qualitative trait evolution
