"""Merge the inline and file specifications of `color_by_metadata` and validate them.

The merged specification is written for the rules that build the auspice configuration,
along with the corresponding `augur export v2` arguments for the `export` rule. Both are
resolved here rather than when the pipeline is parsed, since `color_by_metadata_file`
can be built by a rule and so may not exist at that point.

"""

import shlex
import sys

import yaml

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

valid_scale_types = [
    "viridis_linear",
    "viridis_log",
    "viridis_r_linear",
    "viridis_r_log",
]
valid_options = {
    "scale_type",
    "fixed_min",
    "fixed_max",
    "exclude_auto_scale",
    "title",
}
# these are always in the tree, and `augur export v2` handles them itself
reserved_cols = {"strain", "date"}

inline = snakemake.params.color_by_metadata
print(f"Read {len(inline)} inline `color_by_metadata` column(s): {sorted(inline)}")

if snakemake.params.has_file:
    color_by_metadata_file = snakemake.input.color_by_metadata_file
    with open(color_by_metadata_file) as f:
        from_file = yaml.safe_load(f) or {}
    if not isinstance(from_file, dict):
        raise ValueError(
            f"`color_by_metadata_file` {color_by_metadata_file} must contain a mapping "
            f"of column name to options, but contains a {type(from_file).__name__}"
        )
    print(
        f"Read {len(from_file)} column(s) from `color_by_metadata_file` "
        f"{color_by_metadata_file}: {sorted(from_file)}"
    )
else:
    from_file = {}
    print("No `color_by_metadata_file` specified")

# a column specified in both places is an error rather than one silently winning
in_both = sorted(set(inline) & set(from_file))
if in_both:
    raise ValueError(
        "The following columns are specified both inline in `color_by_metadata` and in "
        f"`color_by_metadata_file`, so it is ambiguous which options apply: {in_both}"
    )

# inline columns are ordered first, as they are the hand-written ones
color_by_metadata = {}
for col, options in list(inline.items()) + list(from_file.items()):
    col = str(col)
    # a column with no options can be written as `col:` rather than `col: {}`
    if options is None:
        options = {}
    if not isinstance(options, dict):
        raise ValueError(
            f"Options for `color_by_metadata` column {col!r} must be a mapping, but "
            f"are a {type(options).__name__}: {options}"
        )
    invalid_options = sorted(set(options) - valid_options)
    if invalid_options:
        raise ValueError(
            f"Invalid option(s) {invalid_options} for `color_by_metadata` column "
            f"{col!r}. Valid options are {sorted(valid_options)}."
        )
    if ("scale_type" in options) and (options["scale_type"] not in valid_scale_types):
        raise ValueError(
            f"For column {col!r}, scale_type {options['scale_type']!r} is not valid. "
            f"Valid options: {valid_scale_types}"
        )
    if ("title" in options) and not isinstance(options["title"], str):
        raise ValueError(
            f"For column {col!r}, 'title' must be a string but is "
            f"{options['title']!r}"
        )
    color_by_metadata[col] = options

# columns are passed to `augur export v2` as shell words, and auspice keys the colorings
# on them, so whitespace in a column name cannot be handled
metadata_columns = [str(c) for c in snakemake.params.metadata_columns]
for col_set, cols in [
    ("color_by_metadata", list(color_by_metadata)),
    ("metadata_columns", metadata_columns),
]:
    cols_w_space = [c for c in cols if len(c.split()) != 1]
    if cols_w_space:
        raise ValueError(
            f"`{col_set}` columns cannot contain whitespace: {cols_w_space}"
        )

# fail here rather than leaving a reserved column to `augur export v2`, whose error
# would not point at the configuration
reserved_specified = sorted(
    reserved_cols.intersection(set(color_by_metadata) | set(metadata_columns))
)
if reserved_specified:
    raise ValueError(
        f"Do not specify {reserved_specified} in `color_by_metadata` or "
        "`metadata_columns`; they are always included in the tree"
    )

# a column that is colored by is already included in the tree, so it does not need to be
# specified as an additional column shown in tooltips
metadata_columns = [c for c in metadata_columns if c not in color_by_metadata]

resolved = {
    "color_by_metadata": color_by_metadata,
    "metadata_columns": metadata_columns,
}
with open(snakemake.output.resolved, "w") as f:
    yaml.safe_dump(resolved, f, sort_keys=False)
print(
    f"\nWrote {len(color_by_metadata)} resolved `color_by_metadata` column(s) and "
    f"{len(metadata_columns)} `metadata_columns` to {snakemake.output.resolved}"
)

export_args = []
for arg, cols in [
    ("--color-by-metadata", list(color_by_metadata)),
    ("--metadata-columns", metadata_columns),
]:
    if cols:
        export_args.append(" ".join([arg, *(shlex.quote(c) for c in cols)]))
export_args = " ".join(export_args)
with open(snakemake.output.export_args, "w") as f:
    f.write(export_args + "\n")
print(f"Wrote these `augur export v2` arguments:\n{export_args}")
