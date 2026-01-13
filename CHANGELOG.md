# CHANGELOG

### version 1.3.0
+ The package is intended to allow multiple identical protein sequences (as long as they have identical strain names), but in practice that throw an error during parsing of branches with no mutations in `scripts/process_treetime_output.py`. Fixed this so now no error is thrown.

### version 1.2.0
+ Rename *example* to *example-flu-seqneut-2025* in the configuration, results, data, etc to allow multiple examples to be added; then added the second example *example-H5NX-seqneutVSVdG*.
+ Eliminate `ruff.toml` in favor of `pyproject.toml`
+ Update `conda` environment to use latest software versions.
+ Allow rooting of trees by date by allowing *outgroup* to be *null*.
+ Add check that no spaces in strain names
+ Allow dates of format `YYYY-MM-DD` where the `MM` and `DD` are optional and can also be `XX`; this is in addition to allowing numeric dates (eg, 2022.5) as before.
+ Added *metadata_columns* to config to show node tooltips not in *color_by_metadata*

### version 1.1.0
+ Build the trees faster with `iqtree`
  - The tree building uses multiple threads if specified when running snakemake.
  - Use the `-fast` option to `iqtree`

+ Handle if `treetime` cannot infer a date

## version 1.0.0
Initial version.
