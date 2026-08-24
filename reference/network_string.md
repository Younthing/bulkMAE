# Retrieve STRING interactions for selected assay features

Retrieve STRING interactions for selected assay features

## Usage

``` r
network_string(
  x,
  experiment,
  assay,
  genes = NULL,
  species = 9606,
  version = "12.0",
  score_threshold = 400,
  input_directory = "",
  remove_unmapped = TRUE
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- genes:

  Optional feature identifiers; defaults to assay row names.

- species:

  NCBI taxonomy identifier.

- version:

  STRING database version.

- score_threshold:

  Minimum STRING score.

- input_directory:

  Optional STRING download/cache directory.

- remove_unmapped:

  Remove identifiers not mapped by STRING.

## Value

A data frame of STRING interactions.
