"""Process results of ``treetime`` output to formats used by ``augur``."""

import json
import re
import sys
from io import StringIO

import Bio.Phylo
import Bio.SeqIO

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

set_branch_lengths_to_n_mutations = snakemake.params.set_branch_lengths_to_n_mutations

print("Reading the trees.")
trees = {}
root_name = None
for treetype in ["timetree", "divtree"]:
    intree = getattr(snakemake.input, f"{treetype}_w_outgroup")

    print(f"Reading tree from {intree}, pruning outgroup")
    tree_w_outgroup = Bio.Phylo.read(intree, "nexus")

    # Bio.Phylo can incorrectly parse names as confidence; fix that.
    for clade in tree_w_outgroup.find_clades(order="preorder"):
        if (clade.name is None) and (clade.confidence is not None):
            clade.name = clade.confidence
            clade.confidence = None

    # Remove outgroup
    if not snakemake.params.no_outgroup:
        print("Removing outgroup")
        outgroup = [c for c in tree_w_outgroup.find_clades() if c.name == "outgroup"]
        if len(outgroup) != 1:
            raise ValueError(
                f"Expected exactly 1 sequence named 'outgroup' in the tree, "
                f"but found {len(outgroup)}. "
                f"Check that the outgroup was correctly added to the alignment."
            )
        # get common ancestor of all non-outgroup tips
        other_tips = [
            t for t in tree_w_outgroup.get_terminals() if t.name != "outgroup"
        ]
        n_total_tips = len(tree_w_outgroup.get_terminals())
        n_other_tips = len(other_tips)
        if n_other_tips + 1 != n_total_tips:
            raise ValueError(
                f"Expected {n_total_tips} total tips to equal {n_other_tips} non-outgroup tips + 1 outgroup tip. "
                f"This indicates the tree structure is invalid or there are multiple sequences named 'outgroup'."
            )
        root = tree_w_outgroup.common_ancestor(other_tips)
        root.branch_length = 0
        root.comment = None
        tree = Bio.Phylo.BaseTree.Tree(root=root, rooted=True)
        if any(t.name == "outgroup" for t in tree.get_terminals()):
            raise ValueError("the outgroup is not actually outgroup to all other tips")
    else:
        tree = tree_w_outgroup
        root = tree.root

    if not root.name:
        raise ValueError(
            f"Tree root has no name. All clades (nodes and tips) in the tree must have names. "
            f"Root object: {root}"
        )
    if root_name:
        if root.name != root_name:
            raise ValueError("inconsistent root names between timetree and divtree")
    else:
        root_name = root.name

    # save tree w comments for later use
    trees[treetype] = tree

# Get the root (reference) sequence
root_seqs = [
    s
    for s in Bio.SeqIO.parse(snakemake.input.ancestral_sequences, "fasta")
    if s.id == root_name
]
if len(root_seqs) != 1:
    raise ValueError(
        f"Expected exactly 1 ancestral sequence for root '{root_name}', "
        f"but found {len(root_seqs)}. "
        f"Each node in the tree must have exactly one corresponding ancestral sequence."
    )
refseq = str(root_seqs[0].seq)
print(f"Using as the reference sequence {root_name}; {len(refseq)=}")

# Get the site numbering map and do some error checking on it
print(f"Reading site numbering map from {snakemake.input.site_numbering_map}")
site_numbering_map = pd.read_csv(snakemake.input.site_numbering_map, sep="\t")
if len(refseq) != len(site_numbering_map):
    raise ValueError(f"{len(refseq)=}, {len(site_numbering_map)=}")
req_cols = {"sequential_site", "protein", "protein_site"}
if "nuc" in set(site_numbering_map["protein"]):
    raise ValueError("Cannot have a protein called 'nuc' in site_numbering_map")
if not req_cols.issubset(site_numbering_map.columns):
    raise ValueError(f"{site_numbering_map.columns=} lacks {req_cols=}")
for col in ["sequential_site", "protein_site"]:
    if not pd.api.types.is_integer_dtype(site_numbering_map[col]):
        raise ValueError(f"{col=} in site_numbering_map not integer")
if not (
    len(site_numbering_map)
    == site_numbering_map["sequential_site"].nunique()
    == site_numbering_map["sequential_site"].max()
    == len(
        range(
            site_numbering_map["sequential_site"].min(),
            site_numbering_map["sequential_site"].max() + 1,
        )
    )
):
    raise ValueError("sequential_site in site_numbering map not 1, 2, ... integers")
if not (
    len(site_numbering_map[["protein", "protein_site"]].drop_duplicates())
    == len(site_numbering_map)
):
    raise ValueError("each sequential_site not mapped to unique protein, protein_site")
refseq_proteins = {}
protein_coords = {}
for protein, protein_df in site_numbering_map.groupby("protein"):
    protein_start = int(protein_df["sequential_site"].min())
    protein_end = int(protein_df["sequential_site"].max())
    protein_coords[protein] = {"start": 3 * protein_start - 2, "end": protein_end * 3}
    if protein_end - protein_start + 1 != len(protein_df):
        raise ValueError(f"{protein=} in site_numbering_map not sequential consecutive")
    if not (protein_start >= 1 and protein_end <= len(refseq)):
        raise ValueError(
            f"Protein '{protein}' has invalid coordinates in site_numbering_map. "
            f"Start site: {protein_start}, End site: {protein_end}. "
            f"Must be in range [1, {len(refseq)}] (reference sequence length)."
        )
    refseq_proteins[protein] = refseq[protein_start - 1 : protein_end]
    if not (
        (protein_df["protein_site"].min() == 1)
        and (protein_df["protein_site"].max() == len(protein_df))
        and (protein_df["protein_site"].nunique() == len(protein_df))
    ):
        raise ValueError(
            f"{protein=} does not have sequential consecutive protein_site"
        )
site_numbering_map = (
    site_numbering_map.assign(protein_site=lambda x: x["protein_site"].astype(int))
    .set_index("sequential_site")
    .to_dict(orient="index")
)

# Get amino-acid mutations in JSON node format and count mutations per node
# This must be done BEFORE writing trees so we can optionally reset branch lengths
print("Parsing mutations from divtree")
aa_muts_nodes = {}
node_mutation_counts = {}

mut_pat = re.compile(r'mutations="(?P<mutations>[^"]*)"')

for clade in trees["divtree"].find_clades(order="preorder"):
    if not clade.name:
        raise ValueError(
            f"Tree contains a clade without a name. All clades (nodes and tips) in the tree "
            f"must have names. Clade object: {clade}"
        )
    # Biopython puts comment text on .comment for Nexus
    comment = getattr(clade, "comment", None)
    n_mutations = 0
    if comment:
        m = mut_pat.search(comment)
        if m:
            mutation_annotations = mut_pat.findall(comment)
            if len(mutation_annotations) != 1:
                raise ValueError(
                    f"Expected exactly one mutation annotation per clade, "
                    f"but found {len(mutation_annotations)} in clade '{clade.name}'. "
                    f"Comment: {comment}"
                )
            muts = {}
            for mut_str in (
                m.group("mutations").split(",") if m.group("mutations") else []
            ):
                m_match = re.fullmatch(
                    r"(?P<parent>[A-Z\-])(?P<site>\d+)(?P<mut>[A-Z\-])", mut_str
                )
                if not m_match:
                    raise ValueError(
                        f"Invalid mutation format '{mut_str}' in clade '{clade.name}'. "
                        f"Expected format: single capital letter or '-', followed by digits, "
                        f"followed by single capital letter or '-' (e.g., 'A123B' or 'A123-'). "
                        f"Full comment: {comment}"
                    )
                site = int(m_match.group("site"))
                if site not in site_numbering_map:
                    raise ValueError(
                        f"invalid {site=} in {mut_str=} not in {site_numbering_map=}"
                    )
                prot = site_numbering_map[site]["protein"]
                prot_site = site_numbering_map[site]["protein_site"]
                if prot not in muts:
                    muts[prot] = []
                muts[prot].append(
                    f"{m_match.group('parent')}{prot_site}{m_match.group('mut')}"
                )
                n_mutations += 1
            if muts:
                aa_muts_nodes.setdefault(clade.name, {})["aa_muts"] = muts
    node_mutation_counts[clade.name] = n_mutations

# Now write the trees with optionally reset branch lengths
print("Writing trees")
for treetype in ["timetree", "divtree"]:
    outtree = getattr(snakemake.output, treetype)
    # Copy the tree by writing to a string buffer and re-reading, rather than using
    # copy.deepcopy, which hits RecursionError on large trees due to the deeply
    # nested Biopython Clade structure exceeding Python's recursion limit.
    buf = StringIO()
    Bio.Phylo.write(trees[treetype], buf, format="newick")
    buf.seek(0)
    tree_nocomments = Bio.Phylo.read(buf, format="newick")

    # Remove comments
    for clade in tree_nocomments.find_clades(order="preorder"):
        clade.comment = None

    # Optionally reset branch lengths to n_mutations / len(refseq) for divtree
    if treetype == "divtree" and set_branch_lengths_to_n_mutations:
        print(f"Resetting branch lengths to mutation counts for {treetype}")
        for clade in tree_nocomments.find_clades(order="preorder"):
            n_muts = node_mutation_counts.get(clade.name, 0)
            clade.branch_length = n_muts / len(refseq) if len(refseq) > 0 else 0

    # Write tree to Newick
    print(f"Writing {treetype} to {outtree}")
    Bio.Phylo.write(tree_nocomments, outtree, format="newick")


# Write the amino-acid mutations and reference sequence
print(f"Writing aa mutations to {snakemake.output.aa_muts}")
aa_muts_d = {
    "annotations": {
        "nuc": {"start": 1, "end": len(refseq) * 3, "strand": "+", "type": "source"},
    }
    | {
        prot: coords | {"strand": "+", "type": "gene"}
        for prot, coords in protein_coords.items()
    },
    "nodes": aa_muts_nodes,
    "reference": {"nuc": "N" * len(refseq) * 3} | refseq_proteins,
}
with open(snakemake.output.aa_muts, "w") as f:
    f.write(json.dumps(aa_muts_d, indent=2))

# Get the branch lengths in JSON node format, in units of mutations
brlens_nodes = {}

dates = (
    # node names are identifiers, and are looked up below by the tree's string names
    pd.read_csv(snakemake.input.dates, sep="\t", dtype={"#node": str})
    .rename(columns={"numeric date": "num_date"})
    .set_index("#node")
    .to_dict(orient="index")
)
for clade in trees["divtree"].find_clades(order="preorder"):
    if clade.branch_length is None:
        raise ValueError(
            f"Clade '{clade.name}' has no branch length. "
            f"All clades in the tree must have branch lengths assigned by treetime. "
            f"Clade object: {clade}"
        )
    # Handle '--' values from treetime (dates that couldn't be inferred)
    num_date = dates[clade.name]["num_date"]
    if num_date == "--":
        num_date = None
    else:
        num_date = float(num_date)

    # Use integer mutation count if set_branch_lengths_to_n_mutations is True,
    # otherwise use the ML-optimized branch length from TreeTime
    if set_branch_lengths_to_n_mutations:
        mutation_length = node_mutation_counts.get(clade.name, 0)
    else:
        mutation_length = clade.branch_length * len(refseq)

    brlens_nodes[clade.name] = {
        "mutation_length": mutation_length,
        "date": dates[clade.name]["date"],
        "num_date": num_date,
    }

print(f"Writing branch lengths to {snakemake.output.brlens}")
with open(snakemake.output.brlens, "w") as f:
    f.write(json.dumps({"nodes": brlens_nodes}, indent=2))
