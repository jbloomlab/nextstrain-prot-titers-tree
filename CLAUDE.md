# Instructions for Claude Code

## Bloom lab coding standards

@bloomlab-coding-standards/CLAUDE.md

The [standards](https://github.com/jbloomlab/bloomlab-coding-standards) are pinned at a
commit of that submodule; update periodically with
`git submodule update --remote bloomlab-coding-standards`.

## Additional guidelines for this project

- This is a pipeline that project repositories include as a git submodule and drive with
  their own configuration, so the lab layout looks different here: configuration is the root
  `config_example-*.yaml` files rather than a single `config.yml`, `data/` holds the example
  inputs and the substitution model rather than a project's own data, and `results/` is
  written per configuration under `results_subdir`. Rules live in `Snakefile`; raise
  splitting them into `rules/*.smk` if it grows appreciably longer than it is now.
- [README](README.md) and the comments in the example configurations are the canonical
  documentation of what the pipeline does and what every configuration key means. Read them
  and do not duplicate them here.
- Including projects pin a tag, so a configuration that worked at the previous tag must keep
  working. Add a new option as an optional key with a default; discuss before breaking one.
- Verify a change with `ruff check .`, `black --check .`, `snakefmt --check .`, and
  `snakemake --configfile config_example-flu-seqneut-2025.yaml --lint`, plus
  `snakemake -n --configfile <config>` for each example configuration, which catches
  configuration and DAG errors cheaply. One rule can be rebuilt on its own with
  `snakemake -j1 <output path> --configfile <config>`. Note that
  `snakemake --forcerun <rule> <target>` does not do what it looks like: the target is read
  as a second argument to `--forcerun`, so `rule all` is built instead.
- Running an example to completion rewrites the git-tracked JSONs in `auspice/`, and a rerun
  differs in labels that are not stable from run to run: the names assigned to internal
  nodes, and the order of `meta.filters`, which `augur` builds from a set. Regenerate them
  deliberately and in their own commit, and check the diff holds nothing but that churn.
- `pyproject.toml` holds the version and is the single source of truth for it. A version
  accumulates changes over several pull requests, and the git tag is what marks it released.

