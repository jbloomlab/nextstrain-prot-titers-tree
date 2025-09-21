"""Check input metadata and alignment and remove special charactes from strain names."""

import collections
import re
import sys

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

alignment = list(Bio.SeqIO.parse(snakemake.input.alignment, format="fasta"))
assert len(alignment), "empty alignment in {snakemake.input.alignment=}"

req_metadata_cols = {"strain", "date"}
metadata = pd.read_csv(snakemake.input.metadata, sep="\t")

if not req_metadata_cols.issubset(metadata.columns):
    raise ValueError(f"{metadata.columns=} lacks {req_metadata_cols=}")

color_by_metadata = snakemake.params.color_by_metadata
assert len(color_by_metadata) == len(set(color_by_metadata))
assert not req_metadata_cols.intersection(color_by_metadata)
if not set(color_by_metadata).issubset(metadata.columns):
    raise ValueError(f"{metadata.columns=} lacks {color_by_metadata=}")

# make sure dates all numeric
dates = pd.to_numeric(metadata["date"], errors="raise")
assert dates.notnull().all(), "date not all non-null"

assert len(metadata) == len(alignment)
assert len(metadata) == metadata["strain"].nunique()

strain_renames = {
    orig: re.sub(r"[^A-Za-z0-9._\-/]", "_", orig.replace("'", ""))
    for orig in metadata["strain"].unique()
}
if len(strain_renames) != len(set(strain_renames.values())):
    raise ValueError(f"re-named strains not unique:\n{strain_renames=}")

metadata.assign(strain=lambda x: x["strain"].map(strain_renames)).to_csv(
    snakemake.output.metadata, sep="\t"
)

outgroup = Bio.SeqIO.read(snakemake.input.outgroup, "fasta")
assert "outgroup" not in set(strain_renames.values())

seqlengths = collections.defaultdict(int)
with open(snakemake.output.alignment, "w") as f:
    for seq in alignment:
        strain = seq.id
        assert strain in strain_renames, strain
        strain = strain_renames[strain]
        seq = str(seq.seq)
        seqlengths[len(seq)] += 1
        f.write(f">{strain}\n{seq}\n")

    if len(seqlengths) != 1:
        raise ValueError(f"Not all sequences same length in alignment:\n{seqlengths=}")
    assert len(outgroup) == list(seqlengths)[0], f"{seqlengths=}, {len(outgroup)=}"
    f.write(f">outgroup\n{str(outgroup.seq)}\n")
