"""Sparse matrix helpers.

We avoid adding SciPy as a hard dependency. Instead we produce CSR arrays directly and
serialize them to an NPZ bundle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np


@dataclass
class CSRBuilder:
    n_rows: int
    n_cols: int
    data: List[float] = field(default_factory=list)
    indices: List[int] = field(default_factory=list)
    indptr: List[int] = field(default_factory=lambda: [0])

    def add_row(self, feats: Dict[int, float]) -> None:
        # Sort columns for deterministic output.
        items = sorted(((int(k), float(v)) for k, v in feats.items() if float(v) != 0.0), key=lambda x: x[0])
        for k, v in items:
            if k < 0 or k >= int(self.n_cols):
                continue
            self.indices.append(int(k))
            self.data.append(float(v))
        self.indptr.append(len(self.data))

    def to_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        data = np.asarray(self.data, dtype=np.float32)
        indices = np.asarray(self.indices, dtype=np.int32)
        indptr = np.asarray(self.indptr, dtype=np.int32)
        shape = np.asarray([int(self.n_rows), int(self.n_cols)], dtype=np.int64)
        return data, indices, indptr, shape


def save_csr_npz(path: str, *, data: np.ndarray, indices: np.ndarray, indptr: np.ndarray, shape: np.ndarray) -> None:
    np.savez_compressed(path, data=data, indices=indices, indptr=indptr, shape=shape)

