"""Check the input alignment, sanitize its strain names, and add the outgroup.

This is kept apart from `prep_and_sanitize_data.py` so that the alignment, and so the
tree built from it, is rewritten only when the alignment or outgroup changes rather than
whenever the metadata or titers do.

"""

import collections
import sys

import Bio.SeqIO

import strain_names

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

alignment = list(Bio.SeqIO.parse(snakemake.input.alignment, format="fasta"))
if not len(alignment):
    raise ValueError(f"Empty alignment in {snakemake.input.alignment}")

# Check the alignment has no repeated sequence name, which would otherwise be
# collapsed when the names are mapped to their sanitized forms
alignment_strains = [seq.id for seq in alignment]
if len(set(alignment_strains)) != len(alignment):
    repeated = [s for s, n in collections.Counter(alignment_strains).items() if n > 1]
    raise ValueError(
        f"Alignment has {len(alignment)} sequences but only "
        f"{len(set(alignment_strains))} distinct names, so some are repeated"
        + strain_names.listed(repeated)
    )

strain_renames = strain_names.renames(alignment_strains)

seqlengths = collections.Counter(len(seq) for seq in alignment)
if len(seqlengths) != 1:
    raise ValueError(f"Not all sequences same length in alignment:\n{seqlengths=}")
alignment_length = next(iter(seqlengths))

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
    if len(outgroup) != alignment_length:
        raise ValueError(
            f"Outgroup sequence length ({len(outgroup)}) does not match "
            f"alignment sequence length ({alignment_length}). "
            f"All sequences must be the same length."
        )

with open(snakemake.output.alignment, "w") as f:
    f.writelines(f">{strain_renames[seq.id]}\n{seq.seq!s}\n" for seq in alignment)
    if not no_outgroup:
        f.write(f">outgroup\n{outgroup.seq!s}\n")
print(f"Wrote {len(alignment)} sequences to {snakemake.output.alignment}")
