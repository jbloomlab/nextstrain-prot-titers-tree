"""Check input metadata and titers and remove special characters from strain names.

The alignment is only read here, to check it has the same strains as the metadata; it is
written by `prep_alignment.py`, so that the tree is not rebuilt when only the metadata or
titers change.

"""

import re
import sys

import Bio.SeqIO

import pandas as pd

import strain_names

import yaml

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input.resolved_color_by_metadata) as f:
    resolved_color_by_metadata = yaml.safe_load(f)

# repeated names, sequence lengths, and an empty alignment are checked by
# `prep_alignment.py`, so just the names are needed here
alignment_strains = {
    seq.id for seq in Bio.SeqIO.parse(snakemake.input.alignment, format="fasta")
}

req_metadata_cols = {"strain", "date"}
# strain names are identifiers, so read them as strings rather than letting an
# all-numeric set of names be read as numbers that no longer match the alignment
metadata = pd.read_csv(snakemake.input.metadata, sep="\t", dtype={"strain": str})

if not req_metadata_cols.issubset(metadata.columns):
    raise ValueError(f"{metadata.columns=} lacks {req_metadata_cols=}")

for col_set in ["color_by_metadata", "metadata_columns"]:
    cols = list(resolved_color_by_metadata[col_set])
    if len(cols) != len(set(cols)):
        raise ValueError(f"Duplicate entries in `{col_set}`: {cols}")
    if req_metadata_cols.intersection(cols):
        raise ValueError(f"Do not include {req_metadata_cols} in `{col_set}` of {cols}")
    if not set(cols).issubset(metadata.columns):
        raise ValueError(f"{metadata.columns=} lacks `{col_set}` {cols}")

# make sure dates all valid
for strain, date in metadata[["strain", "date"]].itertuples(index=False):
    if not (
        isinstance(date, (float, int))
        or (
            isinstance(date, str)
            and re.fullmatch(r"\d\.?\d*|\d{4}(?:\-(?:\d{2}|XX)){0,2}", date)
        )
    ):
        raise ValueError(f"Invalid {date=} for {strain=}")

# Check metadata has unique strains
n_metadata_strains = len(metadata)
n_unique_strains = metadata["strain"].nunique()
if n_metadata_strains != n_unique_strains:
    duplicated = (
        metadata[metadata["strain"].duplicated(keep=False)].groupby("strain").size()
    )
    dup_str = "\n".join(
        f"  {strain}: {count} occurrences" for strain, count in duplicated.items()
    )
    raise ValueError(
        f"Metadata has {n_metadata_strains} rows but only {n_unique_strains} unique strains.\n"
        f"Duplicated strains:\n{dup_str}"
    )

# Check strain names do not contain whitespace
strains_w_space = metadata[metadata["strain"].str.contains(r"\s", regex=True)][
    "strain"
].tolist()
if strains_w_space:
    raise ValueError(
        f"Found {len(strains_w_space)} strain names containing whitespace"
        + strain_names.listed(strains_w_space)
    )

# Check metadata and alignment cover the same strains
metadata_strains = set(metadata["strain"])
if metadata_strains != alignment_strains:
    err_parts = [
        (
            f"Metadata has {len(metadata_strains)} strains and the alignment has "
            f"{len(alignment_strains)}, and their names differ."
        )
    ]
    for descr, missing in [
        ("Strains in metadata but not alignment", metadata_strains - alignment_strains),
        ("Strains in alignment but not metadata", alignment_strains - metadata_strains),
    ]:
        if missing:
            err_parts.append(descr + strain_names.listed(missing))
    raise ValueError("\n".join(err_parts))

strain_renames = strain_names.renames(metadata["strain"].unique())

metadata = metadata.assign(strain=lambda x: x["strain"].map(strain_renames))
# written at full precision, unlike the titers below: these are the values given in the
# input metadata, and they pass through to auspice as the user wrote them
metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)

for i, collection in enumerate(snakemake.params.titer_collection_keys):
    print(f"\nProcessing titer collection {collection}")
    titers = pd.read_csv(
        snakemake.input.titers[i], sep="\t", dtype={"strain": str, "serum": str}
    )
    titer_cols = snakemake.params.titer_cols[i]
    if not set(titer_cols).issubset(titers.columns):
        raise ValueError(f"{titer_cols=} not all in {titers.columns=}")
    extra_strains = set(titers["strain"]) - set(strain_renames)
    if extra_strains:
        raise ValueError(
            f"strains w titers but not specified as part of tree: {extra_strains}"
        )
    titers = titers[titer_cols].assign(strain=lambda x: x["strain"].map(strain_renames))
    titers_not_in_metadata = set(titers["strain"]) - set(metadata["strain"])
    if titers_not_in_metadata:
        raise ValueError(
            f"Found {len(titers_not_in_metadata)} strains with titers that are not in "
            "metadata" + strain_names.listed(titers_not_in_metadata)
        )
    titers_per_strain_serum = (
        titers.groupby(["strain", "serum"])
        .aggregate(
            n_titers=pd.NamedAgg(titer_cols[2], "size"),
            titers=pd.NamedAgg(titer_cols[2], "unique"),
        )
        .sort_values("n_titers", ascending=False)
        .query("n_titers > 1")
    )
    if len(titers_per_strain_serum):
        raise ValueError(
            f"multiple titers for some strains/sera:\n{titers_per_strain_serum}"
        )
    titers.to_csv(
        snakemake.output.titers[i], sep="\t", float_format="%.5g", index=False
    )
