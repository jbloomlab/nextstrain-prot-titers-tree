"""Check input metadata and alignment and remove special charactes from strain names."""

import collections
import re
import sys

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

alignment = list(Bio.SeqIO.parse(snakemake.input.alignment, format="fasta"))
if not len(alignment):
    raise ValueError(f"Empty alignment in {snakemake.input.alignment}")

req_metadata_cols = {"strain", "date"}
metadata = pd.read_csv(snakemake.input.metadata, sep="\t")

if not req_metadata_cols.issubset(metadata.columns):
    raise ValueError(f"{metadata.columns=} lacks {req_metadata_cols=}")

for col_set in ["color_by_metadata", "metadata_columns"]:
    cols = snakemake.params[col_set]
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

# Check metadata and alignment have same length and same strain names
alignment_strains = {seq.id for seq in alignment}
metadata_strains = set(metadata["strain"])

if len(metadata) != len(alignment):
    only_metadata = metadata_strains - alignment_strains
    only_alignment = alignment_strains - metadata_strains
    err_parts = [
        f"Metadata has {len(metadata)} strains but alignment has {len(alignment)} sequences."
    ]
    if only_metadata:
        n_show = 10
        shown = list(only_metadata)[:n_show]
        more_msg = (
            f" (showing first {n_show} of {len(only_metadata)})"
            if len(only_metadata) > n_show
            else ""
        )
        err_parts.append(
            f"Strains in metadata but not alignment{more_msg}:\n  " + "\n  ".join(shown)
        )
    if only_alignment:
        n_show = 10
        shown = list(only_alignment)[:n_show]
        more_msg = (
            f" (showing first {n_show} of {len(only_alignment)})"
            if len(only_alignment) > n_show
            else ""
        )
        err_parts.append(
            f"Strains in alignment but not metadata{more_msg}:\n  " + "\n  ".join(shown)
        )
    raise ValueError("\n".join(err_parts))

if metadata_strains != alignment_strains:
    only_metadata = metadata_strains - alignment_strains
    only_alignment = alignment_strains - metadata_strains
    err_parts = [
        f"Metadata and alignment both have {len(metadata)} entries, but strain names differ."
    ]
    if only_metadata:
        n_show = 10
        shown = list(only_metadata)[:n_show]
        more_msg = (
            f" (showing first {n_show} of {len(only_metadata)})"
            if len(only_metadata) > n_show
            else ""
        )
        err_parts.append(
            f"Strains in metadata but not alignment{more_msg}:\n  " + "\n  ".join(shown)
        )
    if only_alignment:
        n_show = 10
        shown = list(only_alignment)[:n_show]
        more_msg = (
            f" (showing first {n_show} of {len(only_alignment)})"
            if len(only_alignment) > n_show
            else ""
        )
        err_parts.append(
            f"Strains in alignment but not metadata{more_msg}:\n  " + "\n  ".join(shown)
        )
    raise ValueError("\n".join(err_parts))

strains_w_space = metadata[metadata["strain"].str.contains(r"\s", regex=True)]["strain"]
if len(strains_w_space):
    raise ValueError(f"following strain names contain whitespace:\n{strains_w_space}")

strain_renames = {
    orig: re.sub(r"[^A-Za-z0-9._\-/]", "_", orig.replace("'", ""))
    for orig in metadata["strain"].unique()
}
if len(strain_renames) != len(set(strain_renames.values())):
    raise ValueError(f"re-named strains not unique:\n{strain_renames=}")

metadata = metadata.assign(strain=lambda x: x["strain"].map(strain_renames))
metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)

no_outgroup = snakemake.params.no_outgroup
if not no_outgroup:
    outgroup = Bio.SeqIO.read(snakemake.input.outgroup, "fasta")
    if "outgroup" in set(strain_renames.values()):
        raise ValueError(
            'Strain name "outgroup" is reserved for the tree root. '
            "One of your strains is being renamed to 'outgroup', which conflicts with this reserved name. "
            f"Original strain names that map to 'outgroup': "
            f"{[k for k, v in strain_renames.items() if v == 'outgroup']}"
        )

seqlengths = collections.defaultdict(int)
with open(snakemake.output.alignment, "w") as f:
    for seq in alignment:
        strain = seq.id
        if strain not in strain_renames:
            raise ValueError(
                f"Strain '{strain}' found in alignment but not in metadata. "
                "This should not happen as earlier validation checked metadata and alignment have matching strains."
            )
        strain = strain_renames[strain]
        seq = str(seq.seq)
        seqlengths[len(seq)] += 1
        f.write(f">{strain}\n{seq}\n")

    if len(seqlengths) != 1:
        raise ValueError(f"Not all sequences same length in alignment:\n{seqlengths=}")

    if not no_outgroup:
        alignment_length = list(seqlengths)[0]
        outgroup_length = len(outgroup)
        if outgroup_length != alignment_length:
            raise ValueError(
                f"Outgroup sequence length ({outgroup_length}) does not match "
                f"alignment sequence length ({alignment_length}). "
                f"All sequences must be the same length."
            )
        f.write(f">outgroup\n{str(outgroup.seq)}\n")

if snakemake.params.have_titers:
    titers = pd.read_csv(snakemake.input.titers, sep="\t")
    titer_cols = snakemake.params.titer_cols
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
        n_show = 10
        shown = list(titers_not_in_metadata)[:n_show]
        more_msg = (
            f" (showing first {n_show} of {len(titers_not_in_metadata)})"
            if len(titers_not_in_metadata) > n_show
            else ""
        )
        raise ValueError(
            f"Found {len(titers_not_in_metadata)} strains with titers that are not in metadata{more_msg}:\n  "
            + "\n  ".join(shown)
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
    titers.to_csv(snakemake.output.titers, sep="\t", float_format="%.5g", index=False)
