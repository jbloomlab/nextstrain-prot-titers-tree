"""Process results of ``treetime ancestral`` to JSON of amino-acid mutations."""

import json
import re
import sys

import Bio.Phylo
import Bio.SeqIO


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

tree = Bio.Phylo.read(snakemake.input.annotated_tree, "nexus")
ancestral_sequences = Bio.SeqIO.parse(snakemake.input.ancestral_sequences, "fasta")
gene = snakemake.params.gene
aa_muts = snakemake.output.aa_muts

mut_pat = re.compile(r'mutations="(?P<mutations>[^"]*)"')

nodes = {}  # JSON nodes for mutations

for clade in tree.find_clades(order="preorder"):
    if (clade.name is None) and (clade.confidence is not None):
        name = str(clade.confidence)  # Nexus wrapper can parse names as confidence
    else:
        name = clade.name
    if not name:
        raise ValueError(
            f"Tree contains a clade without a name. All clades (nodes and tips) in the tree "
            f"must have names. Clade object: {clade}"
        )
    muts = []
    # Biopython puts comment text on .comment for Nexus
    comment = getattr(clade, "comment", None)
    if comment:
        m = mut_pat.search(comment)
        if m:
            mutation_annotations = mut_pat.findall(comment)
            if len(mutation_annotations) != 1:
                raise ValueError(
                    f"Expected exactly one mutation annotation per clade, "
                    f"but found {len(mutation_annotations)} in clade '{name}'. "
                    f"Comment: {comment}"
                )
            muts = []
            for mut_str in m.group("mutations").split(","):
                if not re.fullmatch(r"[A-Z]\d+[A-Z]", mut_str):
                    raise ValueError(
                        f"Invalid mutation format '{mut_str}' in clade '{name}'. "
                        f"Expected format: single capital letter, followed by digits, followed by single capital letter "
                        f"(e.g., 'A123B'). Full comment: {comment}"
                    )
                muts.append(mut_str)
            if muts:
                nodes.setdefault(name, {})["aa_muts"] = {gene: muts}

# Get the root (reference) sequence
root = tree.root
common_anc = tree.common_ancestor(tree.get_terminals())
if tree.root != common_anc:
    raise ValueError(
        f"Tree root is not the common ancestor of all terminal nodes. "
        f"This indicates the tree structure is invalid. "
        f"Root: {tree.root.name if tree.root.name else tree.root}, "
        f"Common ancestor: {common_anc.name if common_anc.name else common_anc}"
    )
root_seqs = [s for s in ancestral_sequences if s.id == root.name]
if len(root_seqs) != 1:
    raise ValueError(
        f"Expected exactly 1 ancestral sequence for root '{root.name}', "
        f"but found {len(root_seqs)}. "
        f"Each node in the tree must have exactly one corresponding ancestral sequence."
    )
out = {"nodes": nodes, "reference": {gene: str(root_seqs[0].seq)}}

with open(aa_muts, "w") as f:
    f.write(json.dumps(out, indent=2))
