"""
MolProp Toolkit - Similarity Search Module

Provides fingerprint-based similarity searching for compound libraries.
Supports multiple fingerprint types and similarity metrics.
"""

from .fingerprints import (
    get_fingerprint,
    get_fingerprint_bulk,
    FINGERPRINT_TYPES,
)
from .search import (
    similarity_search,
    pairwise_similarity,
    diversity_pick,
)
from .metrics import (
    tanimoto_similarity,
    dice_similarity,
    cosine_similarity,
    SIMILARITY_METRICS,
)

__all__ = [
    "get_fingerprint",
    "get_fingerprint_bulk",
    "FINGERPRINT_TYPES",
    "similarity_search",
    "pairwise_similarity",
    "diversity_pick",
    "tanimoto_similarity",
    "dice_similarity",
    "cosine_similarity",
    "SIMILARITY_METRICS",
]
