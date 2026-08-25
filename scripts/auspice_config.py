"""Build the auspice configuration JSON, including the color scale for each coloring."""

import json
import sys

import matplotlib

import numpy

import pandas as pd

import yaml

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# get discrete colors spaced across color map
n_scale_points = 8  # use this many points to define the color scale
hex_colors = {
    cmap: [
        matplotlib.colors.to_hex(
            matplotlib.colormaps[cmap].resampled(n_scale_points)(i)
        )
        for i in range(n_scale_points)
    ]
    for cmap in ["viridis", "viridis_r"]
}


metadata = pd.read_csv(snakemake.input.metadata, sep="\t", dtype={"strain": str})
with open(snakemake.input.resolved_color_by_metadata) as f:
    color_by_metadata = yaml.safe_load(f)["color_by_metadata"]
missing_cols = set(color_by_metadata) - set(metadata.columns)
if missing_cols:
    raise ValueError(
        f"The following columns from 'color_by_metadata' are not present in metadata: "
        f"{sorted(missing_cols)}. "
        f"Available metadata columns: {sorted(metadata.columns.tolist())}"
    )

colorings = []
for col, col_d in color_by_metadata.items():
    if "scale_type" in col_d:
        # Compute min/max for continuous color scales
        if "exclude_auto_scale" in col_d:
            if not isinstance(col_d["exclude_auto_scale"], list):
                raise ValueError(
                    f"For column '{col}', 'exclude_auto_scale' must be a list, "
                    f"but got type {type(col_d['exclude_auto_scale']).__name__}: "
                    f"{col_d['exclude_auto_scale']}"
                )
            # strain names here have had special characters replaced by
            # `prep_and_sanitize_data.py`, so a name that does not match is either a
            # typo or written in its original form
            unmatched = sorted(
                set(col_d["exclude_auto_scale"]) - set(metadata["strain"])
            )
            if unmatched:
                raise ValueError(
                    f"For column '{col}', these 'exclude_auto_scale' strains are not "
                    f"in the metadata: {unmatched}. Names must match the metadata "
                    "strain names with special characters replaced by '_'."
                )
            for_lims = metadata[~metadata["strain"].isin(col_d["exclude_auto_scale"])]
        else:
            for_lims = metadata
        try:
            minval = (
                col_d["fixed_min"] if ("fixed_min" in col_d) else for_lims[col].min()
            )
            maxval = (
                col_d["fixed_max"] if ("fixed_max" in col_d) else for_lims[col].max()
            )
        except TypeError as e:
            # Check for mixed types
            value_types = metadata[col].dropna().apply(type).unique()
            sample_values = metadata[col].dropna().head(10).tolist()
            raise ValueError(
                f"Failed to compute min/max for column '{col}' with scale_type '{col_d['scale_type']}'. "
                f"The column contains mixed data types: {[t.__name__ for t in value_types]}. "
                f"Columns with continuous scale_type must contain only numeric values. "
                f"Sample values from column: {sample_values}. "
                f"Check your data to ensure '{col}' contains only numeric values, not a mix of numbers and strings."
            ) from e
        except ValueError as e:
            # Other value errors (e.g., empty column)
            raise ValueError(
                f"Failed to compute min/max for column '{col}' with scale_type '{col_d['scale_type']}'. "
                f"Original error: {e}"
            ) from e
        if minval > metadata[col].min():
            minprefix = "<="
        else:
            minprefix = ""
        if maxval < metadata[col].max():
            maxprefix = ">="
        else:
            maxprefix = ""
        if maxval < minval:
            raise ValueError(
                f"For column '{col}', maximum value ({maxval}) is less than minimum value ({minval}). "
                f"Check 'fixed_min' and 'fixed_max' settings or the data values in this column."
            )

        scale_type = col_d["scale_type"]
        if scale_type.endswith(("_linear", "_log")):
            scale = scale_type.split("_")[-1]
            cmap = "_".join(scale_type.split("_")[:-1])
        else:
            raise ValueError(f"invalid {scale_type=}")
        if cmap in hex_colors:
            hexcols = hex_colors[cmap]
        else:
            raise ValueError(f"invalid {scale_type=}")
        if scale == "linear":
            scalevals = numpy.linspace(minval, maxval, num=len(hexcols))
        elif scale == "log":
            if minval <= 0:
                positive_vals = for_lims[col][for_lims[col] > 0]
                if positive_vals.empty:
                    raise ValueError(
                        f"For column '{col}' with log scale, all values are <= 0. "
                        f"Cannot compute a log scale with no positive values."
                    )
                minval = positive_vals.min()
                minprefix = "<="
            scalevals = numpy.logspace(
                numpy.log(minval) / numpy.log(2),
                numpy.log(maxval) / numpy.log(2),
                num=len(hexcols),
                base=2,
            )
        else:
            raise ValueError(f"{scale=}")
        legendlabels = [f"{v:.3g}" for v in scalevals]
        legendlabels[0] = minprefix + legendlabels[0]
        legendlabels[-1] = maxprefix + legendlabels[-1]
        color_scale = {
            "type": "continuous",
            "scale": [[float(v), c] for v, c in zip(scalevals, hexcols)],
            "legend": [
                {"value": v, "display": d} for v, d in zip(scalevals, legendlabels)
            ],
        }
    else:
        color_scale = {}
    # without a 'title', auspice labels the coloring with the column name itself
    title = {"title": col_d["title"]} if "title" in col_d else {}
    colorings.append({"key": col} | title | color_scale)

json_d = {
    "display_defaults": snakemake.params.display_defaults,
    "colorings": colorings,
    "panels": (
        ["tree", "entropy"] + (["measurements"] if snakemake.params.has_titers else [])
    ),
}

with open(snakemake.output.auspice_config, "w") as f:
    json.dump(json_d, f, indent=2)
