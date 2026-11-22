---
title: "gbscan pipeline: Getting COI sequences from GenBank"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get COI sequences for the scorpion genus Tityus}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

* TOC
{:toc}

This vignette shows the full gbscan pipeline. The example gets the COI
sequences of the scorpion genus *Tityus* (Buthidae) from GenBank. The
pipeline has five steps:

1. [Examine the marker dictionary](#step-1-examine-the-marker-dictionary).
2. [Count the records](#step-2-count-the-records).
3. [Download the records](#step-3-download-the-records).
4. [Read the metadata and filter the records](#step-4-read-the-metadata-and-filter-the-records).
5. [Write the sequences to FASTA files](#step-5-write-the-sequences-to-fasta-files).

The code in this vignette needs an internet connection. The vignette
does not run the code when the package is built.

## Installation

Install the development version of gbscan from GitHub:

``` r
pak::pak("https://github.com/foersterst/gbscan.git")
```

Then, load the required packages:

``` r
library(gbscan)
library(dplyr)
```

## Before you start: get an NCBI API key

NCBI accepts 3 requests per second without a key. With a key, NCBI
accepts 10 requests per second. Large downloads are faster with a key.

1.  Make an NCBI account at <https://www.ncbi.nlm.nih.gov/account/>.
2.  Go to the account settings and make an API key.
3.  Open your `.Renviron` file with `usethis::edit_r_environ()`.
4.  Add the line `ENTREZ_KEY=your_key`.
5.  Save the file and restart R.

Do not write the key in a script. A script can go to other persons or to
a public repository.

## Step 1: examine the marker dictionary

The dictionary tells gbscan how to search for each marker and how to
identify the marker in a record.

``` r
dict <- gb_dictionary()
dict
#> # A tibble: 29 × 7
#>    marker genome        description                       query                                                         regex min_len max_len
#>    <chr>  <chr>         <chr>                             <chr>                                                         <chr>   <dbl>   <dbl>
#>  1 COI    mitochondrial cytochrome c oxidase subunit I    "COI[Gene Name] OR COX1[Gene Name] OR CO1[Gene Name] OR COXI… "^(c…     100    5000
#>  2 COII   mitochondrial cytochrome c oxidase subunit II   "COII[Gene Name] OR COX2[Gene Name] OR CO2[Gene Name] OR COX… "^(c…     100    5000
#>  3 COIII  mitochondrial cytochrome c oxidase subunit III  "COIII[Gene Name] OR COX3[Gene Name] OR CO3[Gene Name] OR CO… "^(c…     100    5000
#>  4 CytB   mitochondrial cytochrome b                      "CYTB[Gene Name] OR COB[Gene Name] OR \"cytochrome b\"[Title… "^(c…     100    5000
#>  5 ND1    mitochondrial NADH dehydrogenase subunit 1      "ND1[Gene Name] OR NAD1[Gene Name] OR \"NADH dehydrogenase s… "^(n…     100    5000
#>  6 ND4    mitochondrial NADH dehydrogenase subunit 4      "ND4[Gene Name] OR NAD4[Gene Name] OR \"NADH dehydrogenase s… "^(n…     100    5000
#>  7 ND5    mitochondrial NADH dehydrogenase subunit 5      "ND5[Gene Name] OR NAD5[Gene Name] OR \"NADH dehydrogenase s… "^(n…     100    5000
#>  8 ATP6   mitochondrial ATP synthase F0 subunit 6         "ATP6[Gene Name] OR ATPASE6[Gene Name] OR \"ATP synthase F0 … "^(a…     100    5000
#>  9 12S    mitochondrial 12S ribosomal RNA (small subunit) "12S[Gene Name] OR rrnS[Gene Name] OR s-rRNA[Gene Name] OR \… "^(1…     100    5000
#> 10 16S    mitochondrial 16S ribosomal RNA (large subunit) "16S[Gene Name] OR rrnL[Gene Name] OR l-rRNA[Gene Name] OR \… "^(1…     100    5000
#> # ℹ 19 more rows

dict |>
  filter(marker == "COI") |>
  pull(query)
#> [1] "COI[Gene Name] OR COX1[Gene Name] OR CO1[Gene Name] OR COXI[Gene Name] OR \"cytochrome oxidase subunit I\"[Title] OR \"cytochrome c oxidase subunit I\"[Title] OR \"cytochrome oxidase subunit 1\"[Title] OR \"cytochrome c oxidase subunit 1\"[Title]"
```

The `query` column has the GenBank search terms. The `regex` column has
the pattern that identifies COI in the features of a record.

To change the dictionary, write a copy to a CSV file, edit the file, and
read it again:

``` r
gb_dictionary("my_markers.csv", write = TRUE)
# Edit my_markers.csv in a text editor or spreadsheet program.
dict <- gb_dictionary("my_markers.csv")
```

## Step 2: count the records

Count the records before you download them. Then you know the size of
the download.

``` r
counts <- gb_count(taxa = "Tityus", markers = c("COI", "mitogenome"))
#> # A tibble: 2 × 4
#>   taxon  marker     query                         count
#>   <chr>  <chr>      <chr>                         <int>
#> 1 Tityus COI        "\"Tityus\"[Organism] AND..."    89
#> 2 Tityus mitogenome "\"Tityus\"[Organism] AND..."     2
```

The function prints the counts with shortened queries. The saved result
has the full queries. You can paste a query into the GenBank web page to
examine the results:

``` r
counts$query
#> [1] "\"Tityus\"[Organism] AND (COI[Gene Name] OR COX1[Gene Name] OR CO1[Gene Name] OR COXI[Gene Name] OR \"cytochrome oxidase subunit I\"[Title] OR \"cytochrome c oxidase subunit I\"[Title] OR \"cytochrome oxidase subunit 1\"[Title] OR \"cytochrome c oxidase subunit 1\"[Title]) AND 100:5000[SLEN]"
#> [2] "\"Tityus\"[Organism] AND (mitochondrion[Filter] AND \"complete genome\"[Title]) AND 10000:30000[SLEN]"
```

The `mitogenome` marker finds complete mitochondrial genomes. Many
mitogenomes do not have the COI terms in their title. Include
`mitogenome` to get the COI sequences from these genomes too.

## Step 3: download the records

``` r
gb_download(
  taxa = "Tityus",
  markers = c("COI", "mitogenome"),
  dir = "tityus_data",
  batch = 500
)
#> 
#> ── Tityus | COI ──
#> 
#> ℹ Searching GenBank✔ Searching GenBank [2.3s]
#> ℹ Found 89 records.✔ Saved the records in 'tityus_data/genbank/COI/Tityus'.
#> ── Tityus | mitogenome ──
#> 
#> ℹ Searching GenBank✔ Searching GenBank [962ms]
#> ℹ Found 2 records.✔ Saved the records in 'tityus_data/genbank/mitogenome/Tityus'.
#> ── NCBI Taxonomy ──
#> 
#> ✔ Downloaded taxonomy data for 49 taxa.
```

The function saves the files in `tityus_data/genbank/`. It also
downloads the NCBI Taxonomy data for the records.

If the download stops, run the same command again. The function
continues from the last complete batch. This is important for large
downloads (for example, tens of thousands of records).

To get records that were added to GenBank after your download, remove
the file `ids.txt` in the folder of the query. Then run the command
again. The function downloads only the new records.

## Step 4: read the metadata and filter the records

``` r
meta <- gb_parse_meta("tityus_data")
#> ✔ Read 91 records and 140 features.
#> ℹ 31 features have no marker.
records <- meta$records
features <- meta$features
```

`records` has one row for each record. It has the taxonomic ranks
(species, genus, tribe, subfamily, family, order, class, phylum), the
specimen data, and the reference data.

`features` has one row for each gene-level feature. The `marker` column
shows the marker that gbscan found for the feature.

NCBI Taxonomy often has no tribe or subfamily for invertebrates. Use the
`ranks` argument to fill the empty ranks from a second database, if
those are available there.

``` r
meta <- gb_parse_meta("tityus_data", ranks = "gbif")
#> ℹ Saved the gbif lookups in 'tityus_data/ranks_gbif.csv'.
#> ✔ Filled ranks of 0 records from "gbif".
#> ✔ Read 91 records and 140 features.
#> ℹ 31 features have no marker.
records <- meta$records

# not available on GBIF
records |>
  count(subfamily, filled_ranks)
#> # A tibble: 1 × 3
#>   subfamily filled_ranks     n
#>   <chr>     <chr>        <int>
#> 1 <NA>      <NA>            91
```

The function looks up each genus one time and saves the results in
`tityus_data/ranks_gbif.csv`. It fills only empty ranks. The column
`filled_ranks` shows the ranks that came from GBIF.

Examine the features without a marker:

``` r
features |>
  filter(is.na(marker)) |>
  count(gene, product, sort = TRUE)
#> # A tibble: 7 × 3
#>   gene  product                           n
#>   <chr> <chr>                         <int>
#> 1 <NA>  <NA>                             21
#> 2 ND2   NADH dehydrogenase subunit 2      2
#> 3 ND3   NADH dehydrogenase subunit 3      2
#> 4 ND4L  NADH dehydrogenase subunit 4L     2
#> 5 ND6   NADH dehydrogenase subunit 6      2
#> 6 ATP8  ATP synthase F0 subunit 8         1
#> 7 ATP8  ATP synthase subunit 8            1
```

If a feature is COI but has no marker, change the COI regular expression
in the dictionary. Then run `gb_parse_meta()` again. You do not have to
download again.

Join the two tables and select the COI sequences:

``` r
coi <- features |>
  filter(marker == "COI") |>
  left_join(records, by = "accession")
```

Keep only identified species, sequences of 500 bp or more, and the
longest sequence of each species:

``` r
best <- coi |>
  filter(
    !is.na(species),
    !grepl(" (sp|cf|aff|nr)\\.", species),
    feature_length >= 500
  ) |>
  arrange(desc(feature_length)) |>
  distinct(species, .keep_all = TRUE)

nrow(best)
#> [1] 42
```

Save the tables for later use:

``` r
saveRDS(meta, "tityus_meta.rds")
```

## Step 5: write the sequences to FASTA files

Write the selected sequences. The `keep` argument gives the features to
write. The `header` argument gives the columns for the FASTA headers.

``` r
gb_parse_fasta(
  dir = "tityus_data",
  out_dir = "fasta",
  keep = best,
  header = c("accession", "species", "country")
)
```

The function writes the file `fasta/COI/COI.fasta`. Each marker gets its
own folder.

The function cuts the COI gene from the mitogenomes. To write the full
record sequences, use `region = "record"`.

Align the sequences with an alignment program, for example MAFFT:

``` sh
mafft --auto fasta/COI/COI.fasta > fasta/COI/COI_aligned.fasta
```

## More taxa and markers

All functions accept more than one taxon and more than one marker. For
example, this command downloads the common Sanger markers for two
genera:

``` r
gb_download(
  taxa = c("Tityus", "Centruroides"),
  markers = c("COI", "12S", "16S", "18S", "28S", "mitogenome"),
  dir = "buthidae_data"
)
```
