"""Strain-name handling shared by `prep_alignment.py` and `prep_and_sanitize_data.py`.

Both rename strains with `sanitize`, the alignment in one and the metadata and titers in
the other, so the names they write must come from this one function to match.

"""

import re

n_show = 10  # how many offending names an error message lists


def sanitize(name):
    """Name with quotes dropped and other special characters replaced by `_`."""
    return re.sub(r"[^A-Za-z0-9._\-/]", "_", name.replace("'", ""))


def renames(names):
    """Map each of `names` to its sanitized name, which must be unique."""
    renamed = {name: sanitize(name) for name in names}
    if len(renamed) != len(set(renamed.values())):
        raise ValueError(f"re-named strains not unique:\n{renamed=}")
    return renamed


def listed(names):
    """One-per-line listing of `names` for an error message, showing `n_show` at most."""
    names = sorted(names)
    truncated = (
        f" (showing first {n_show} of {len(names)})" if len(names) > n_show else ""
    )
    return truncated + ":\n  " + "\n  ".join(names[:n_show])
