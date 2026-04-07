"""MolProp Toolkit - Similarity Search Module.

This subpackage provides fingerprint generation + similarity metrics + search utilities.

Design notes
------------
- Importing this package should not require RDKit; RDKit is only required when
  you actually *compute* fingerprints or similarities.
- The public API re-exports the functions used by the CLI and documented in
  docs/tools/similarity.md.
"""

from __future__ import annotations

from .fingerprints import FINGERPRINT_TYPES, get_fingerprint, get_fingerprint_bulk
from .metrics import (
    SIMILARITY_METRICS,
    bulk_similarity,
    cosine_similarity,
    dice_similarity,
    get_similarity_function,
    tanimoto_similarity,
)
from .search import (
    cluster_by_similarity,
    diversity_pick,
    find_nearest_neighbors,
    pairwise_similarity,
    similarity_search,
)

__all__ = [
    # fingerprints
    "get_fingerprint",
    "get_fingerprint_bulk",
    "FINGERPRINT_TYPES",
    # metrics
    "tanimoto_similarity",
    "dice_similarity",
    "cosine_similarity",
    "get_similarity_function",
    "bulk_similarity",
    "SIMILARITY_METRICS",
    # search utilities
    "similarity_search",
    "pairwise_similarity",
    "diversity_pick",
    "cluster_by_similarity",
    "find_nearest_neighbors",
]
