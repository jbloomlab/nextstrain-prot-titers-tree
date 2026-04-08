# CHANGELOG

### version 1.6.0
+ Previously gaps (deletions) were labeled and counted as mutations on tip nodes, but not inferred on ancestral nodes. Change from using `--gtr jtt92` for `treetime` to a custom model that includes gaps so that they are inferred and counted like amino-acid mutations.

#### version 1.5.1
+ Avoid bug that causes crashes for deep trees via a recursion limit for `copy.deepcopy` in a script; instead just make copy by writing and reading tree as new object.
+ Make sure snakemake pipeline correctly tracks as input any `description` in `addtl_export_args`.

### version 1.5.0
+ Update `snakemake` to version 9.19 and `augur` to version 33.0.1.
+ Set `black` target version in `pyproject.toml`

### version 1.4.0
+ Set branch lengths to be exactly equal to number of amino-acid substitutions by rescaling IQ-TREE output to match `timetree` inferred number of mutation to (addresses [this issue](https://github.com/jbloomlab/nextstrain-prot-titers-tree/issues/12)). This behavior is implemented if the `set_branch_lengths_to_n_mutations` flag is either set to true in the configuration YAML or just excluded altogether. To get the old behavior, you must add to config: `set_branch_lengths_to_n_mutations: false`.

#### version 1.3.1
+ Make error messages in Python scripts more informative.
+ Better handle computing min / max for color scales (can avoid some errors otherwise triggered).

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
