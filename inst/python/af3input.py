"""Optional Python helpers for AlphaFold 3 JSON generation.

The R package is the primary interface. This module is intentionally small so it
can be called from reticulate or command-line glue later without duplicating the
R package's table-handling logic.
"""

from __future__ import annotations

import json
from pathlib import Path


_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                            "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def reverse_complement_dna(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return sequence.translate(_COMPLEMENT)[::-1].upper()


def write_json(job: dict, path: str | Path, pretty: bool = True) -> Path:
    """Write one AlphaFold 3 job dictionary to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(job, handle, indent=2 if pretty else None)
        handle.write("\n")
    return path
