# CHANGELOG

### version 1.2.0
+ Rename *example* to *example-flu-seqneut-2025* in the configuration, results, data, etc to allow multiple examples to be added.
+ Eliminate `ruff.toml` in favor of `pyproject.toml`
+ Update `conda` environment to use latest software versions.

### version 1.1.0
+ Build the trees faster with `iqtree`
  - The tree building uses multiple threads if specified when running snakemake.
  - Use the `-fast` option to `iqtree`

+ Handle if `treetime` cannot infer a date

## version 1.0.0
Initial version.
