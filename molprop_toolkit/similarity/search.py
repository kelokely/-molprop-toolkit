"""
Similarity search functions for molecular libraries.

Provides efficient similarity searching, pairwise comparison matrices,
and diversity-based compound selection.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .fingerprints import get_fingerprint, get_fingerprint_bulk
from .metrics import bulk_similarity, get_similarity_function

try:
    from rdkit import Chem

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@dataclass
class SimilarityResult:
    """Container for a single similarity search result."""

    rank: int
    compound_id: str
    smiles: str
    similarity: float
    fingerprint: Any = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Rank": self.rank,
            "Compound_ID": self.compound_id,
            "SMILES": self.smiles,
            "Similarity": self.similarity,
        }


@dataclass
class SearchResults:
    """Container for similarity search results."""

    query_smiles: str
    query_id: Optional[str]
    fp_type: str
    metric: str
    threshold: float
    results: List[SimilarityResult]
    n_searched: int
    n_valid: int

    def to_dataframe(self) -> pd.DataFrame:
        """Convert results to pandas DataFrame."""
        if not self.results:
            return pd.DataFrame(columns=["Rank", "Compound_ID", "SMILES", "Similarity"])
        return pd.DataFrame([r.to_dict() for r in self.results])

    def to_csv(self, path: str, **kwargs):
        """Save results to CSV file."""
        df = self.to_dataframe()
        df.to_csv(path, index=False, **kwargs)

    def __len__(self):
        return len(self.results)

    def __iter__(self):
        return iter(self.results)


def similarity_search(
    query: Union[str, "Chem.Mol"],
    library: Union[str, pd.DataFrame, List[str]],
    fp_type: str = "morgan",
    metric: str = "tanimoto",
    threshold: float = 0.0,
    top_n: Optional[int] = None,
    smiles_col: str = "SMILES",
    id_col: Optional[str] = "Compound_ID",
    query_id: Optional[str] = None,
    include_query: bool = True,
    fp_params: Optional[Dict] = None,
    n_jobs: int = 1,
    show_progress: bool = False,
) -> SearchResults:
    """
    Search for similar compounds in a library.

    Parameters
    ----------
    query : str or Mol
        Query SMILES or RDKit Mol object
    library : str, DataFrame, or list
        Library to search. Can be:
        - Path to SMILES file or CSV
        - pandas DataFrame with SMILES column
        - List of SMILES strings
    fp_type : str
        Fingerprint type (default: 'morgan')
    metric : str
        Similarity metric (default: 'tanimoto')
    threshold : float
        Minimum similarity threshold (default: 0.0)
    top_n : int, optional
        Return only top N results
    smiles_col : str
        Column name for SMILES in DataFrame (default: 'SMILES')
    id_col : str, optional
        Column name for compound IDs (default: 'Compound_ID')
    query_id : str, optional
        Identifier for the query compound
    include_query : bool
        Include query in results if found in library (default: False)
    fp_params : dict, optional
        Override fingerprint parameters
    n_jobs : int
        Number of parallel jobs for fingerprint generation
    show_progress : bool
        Show progress bar

    Returns
    -------
    SearchResults
        Search results object with ranked compounds

    Examples
    --------
    >>> results = similarity_search(
    ...     "CCO",
    ...     "library.smi",
    ...     fp_type="morgan",
    ...     threshold=0.7,
    ...     top_n=10
    ... )
    >>> print(results.to_dataframe())
    >>> results.to_csv("similar_compounds.csv")
    """
    fp_params = fp_params or {}

    # Generate query fingerprint
    query_fp = get_fingerprint(query, fp_type=fp_type, **fp_params)
    if query_fp is None:
        raise ValueError(f"Could not generate fingerprint for query: {query}")

    # Normalize query SMILES
    query_smiles = query if isinstance(query, str) else Chem.MolToSmiles(query)

    # Load library
    smiles_list, id_list = _load_library(library, smiles_col, id_col)
    n_searched = len(smiles_list)

    # Generate library fingerprints
    library_fps = get_fingerprint_bulk(
        smiles_list,
        fp_type=fp_type,
        n_jobs=n_jobs,
        show_progress=show_progress,
        **fp_params,
    )

    # Calculate similarities
    similarities = bulk_similarity(query_fp, library_fps, metric=metric)

    # Count valid molecules
    n_valid = np.sum(~np.isnan(similarities))

    # Create results list
    results = []
    for idx in np.argsort(similarities)[::-1]:  # Sort descending
        sim = similarities[idx]

        # Skip invalid molecules
        if np.isnan(sim):
            continue

        # Apply threshold
        if sim < threshold:
            break

        # Optionally skip query compound
        if not include_query and smiles_list[idx] == query_smiles:
            continue

        compound_id = id_list[idx] if id_list else f"Compound_{idx+1}"

        results.append(
            SimilarityResult(
                rank=len(results) + 1,
                compound_id=compound_id,
                smiles=smiles_list[idx],
                similarity=float(sim),
                fingerprint=library_fps[idx],
            )
        )

        # Apply top_n limit
        if top_n is not None and len(results) >= top_n:
            break

    return SearchResults(
        query_smiles=query_smiles,
        query_id=query_id,
        fp_type=fp_type,
        metric=metric,
        threshold=threshold,
        results=results,
        n_searched=n_searched,
        n_valid=n_valid,
    )


def pairwise_similarity(
    molecules: Union[str, pd.DataFrame, List[str]],
    fp_type: str = "morgan",
    metric: str = "tanimoto",
    smiles_col: str = "SMILES",
    fp_params: Optional[Dict] = None,
    n_jobs: int = 1,
    show_progress: bool = False,
) -> Tuple[np.ndarray, List[str]]:
    """
    Calculate pairwise similarity matrix for a set of molecules.

    Parameters
    ----------
    molecules : str, DataFrame, or list
        Set of molecules to compare
    fp_type : str
        Fingerprint type
    metric : str
        Similarity metric
    smiles_col : str
        Column name for SMILES
    fp_params : dict, optional
        Override fingerprint parameters
    n_jobs : int
        Number of parallel jobs
    show_progress : bool
        Show progress bar

    Returns
    -------
    tuple
        (similarity_matrix, smiles_list) where similarity_matrix[i,j] is
        the similarity between molecule i and molecule j

    Examples
    --------
    >>> matrix, smiles = pairwise_similarity(["CCO", "CCCO", "CCCCO"])
    >>> print(matrix)
    """
    fp_params = fp_params or {}

    # Load molecules
    smiles_list, _ = _load_library(molecules, smiles_col, None)
    n = len(smiles_list)

    # Generate fingerprints
    fps = get_fingerprint_bulk(
        smiles_list,
        fp_type=fp_type,
        n_jobs=n_jobs,
        show_progress=show_progress,
        **fp_params,
    )

    # Initialize similarity matrix
    sim_matrix = np.zeros((n, n))
    np.fill_diagonal(sim_matrix, 1.0)

    # Calculate pairwise similarities (upper triangle)
    sim_func = get_similarity_function(metric)

    for i in range(n):
        if fps[i] is None:
            sim_matrix[i, :] = np.nan
            sim_matrix[:, i] = np.nan
            continue

        for j in range(i + 1, n):
            if fps[j] is None:
                sim_matrix[i, j] = np.nan
                sim_matrix[j, i] = np.nan
            else:
                sim = sim_func(fps[i], fps[j])
                sim_matrix[i, j] = sim
                sim_matrix[j, i] = sim

    return sim_matrix, smiles_list


def diversity_pick(
    molecules: Union[str, pd.DataFrame, List[str]],
    n_pick: int,
    seed: Optional[Union[str, int]] = None,
    fp_type: str = "morgan",
    metric: str = "tanimoto",
    smiles_col: str = "SMILES",
    id_col: Optional[str] = "Compound_ID",
    fp_params: Optional[Dict] = None,
    n_jobs: int = 1,
    show_progress: bool = False,
) -> pd.DataFrame:
    """
    Select a diverse subset of compounds using MaxMin algorithm.

    Uses a greedy algorithm to iteratively select the compound that
    maximizes the minimum distance to already-selected compounds.

    Parameters
    ----------
    molecules : str, DataFrame, or list
        Set of molecules to pick from
    n_pick : int
        Number of compounds to select
    seed : str, int, or None
        Seed compound (SMILES or index). If None, starts with random.
    fp_type : str
        Fingerprint type
    metric : str
        Similarity metric (will use 1-similarity as distance)
    smiles_col : str
        Column name for SMILES
    id_col : str, optional
        Column name for compound IDs
    fp_params : dict, optional
        Override fingerprint parameters
    n_jobs : int
        Number of parallel jobs
    show_progress : bool
        Show progress bar

    Returns
    -------
    DataFrame
        Selected compounds with columns: Compound_ID, SMILES, Pick_Order, Min_Distance

    Examples
    --------
    >>> diverse_set = diversity_pick("library.smi", n_pick=50)
    >>> print(diverse_set)
    """
    fp_params = fp_params or {}

    # Load molecules
    smiles_list, id_list = _load_library(molecules, smiles_col, id_col)
    n = len(smiles_list)

    if n_pick > n:
        raise ValueError(f"Cannot pick {n_pick} compounds from library of {n}")

    # Generate fingerprints
    fps = get_fingerprint_bulk(
        smiles_list,
        fp_type=fp_type,
        n_jobs=n_jobs,
        show_progress=show_progress,
        **fp_params,
    )

    # Filter out invalid molecules
    valid_indices = [i for i, fp in enumerate(fps) if fp is not None]
    if len(valid_indices) < n_pick:
        raise ValueError(
            f"Only {len(valid_indices)} valid molecules, cannot pick {n_pick}"
        )

    sim_func = get_similarity_function(metric)

    # Initialize selection
    if seed is None:
        # Random seed
        np.random.seed(42)
        picked_idx = [valid_indices[np.random.randint(len(valid_indices))]]
    elif isinstance(seed, int):
        picked_idx = [seed]
    else:
        # Find seed SMILES
        try:
            picked_idx = [smiles_list.index(seed)]
        except ValueError as err:
            raise ValueError(f"Seed SMILES not found in library: {seed}") from err

    picked_set = set(picked_idx)
    min_distances = np.full(n, np.inf)

    # Update min distances from first pick
    for i in valid_indices:
        if fps[i] is not None and i not in picked_set:
            sim = sim_func(fps[picked_idx[0]], fps[i])
            min_distances[i] = 1.0 - sim

    # Iteratively pick compounds
    while len(picked_idx) < n_pick:
        # Find compound with maximum minimum distance
        best_idx = None
        best_dist = -1

        for i in valid_indices:
            if i not in picked_set and min_distances[i] > best_dist:
                best_dist = min_distances[i]
                best_idx = i

        if best_idx is None:
            break

        picked_idx.append(best_idx)
        picked_set.add(best_idx)

        # Update min distances
        for i in valid_indices:
            if fps[i] is not None and i not in picked_set:
                sim = sim_func(fps[best_idx], fps[i])
                dist = 1.0 - sim
                if dist < min_distances[i]:
                    min_distances[i] = dist

    # Build results DataFrame
    results = []
    for order, idx in enumerate(picked_idx, 1):
        compound_id = id_list[idx] if id_list else f"Compound_{idx+1}"
        results.append(
            {
                "Compound_ID": compound_id,
                "SMILES": smiles_list[idx],
                "Pick_Order": order,
                "Min_Distance": min_distances[idx] if order > 1 else 0.0,
            }
        )

    return pd.DataFrame(results)


def cluster_by_similarity(
    molecules: Union[str, pd.DataFrame, List[str]],
    threshold: float = 0.7,
    fp_type: str = "morgan",
    metric: str = "tanimoto",
    smiles_col: str = "SMILES",
    id_col: Optional[str] = "Compound_ID",
    fp_params: Optional[Dict] = None,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """
    Cluster molecules by similarity using leader-follower algorithm.

    Parameters
    ----------
    molecules : str, DataFrame, or list
        Set of molecules to cluster
    threshold : float
        Similarity threshold for cluster membership
    fp_type : str
        Fingerprint type
    metric : str
        Similarity metric
    smiles_col : str
        Column name for SMILES
    id_col : str, optional
        Column name for compound IDs
    fp_params : dict, optional
        Override fingerprint parameters
    n_jobs : int
        Number of parallel jobs

    Returns
    -------
    DataFrame
        Molecules with cluster assignments
    """
    fp_params = fp_params or {}

    # Load molecules
    smiles_list, id_list = _load_library(molecules, smiles_col, id_col)
    n = len(smiles_list)

    # Generate fingerprints
    fps = get_fingerprint_bulk(smiles_list, fp_type=fp_type, n_jobs=n_jobs, **fp_params)

    sim_func = get_similarity_function(metric)

    # Leader-follower clustering
    cluster_ids = [-1] * n
    cluster_centroids = []  # (index, fingerprint)

    for i in range(n):
        if fps[i] is None:
            continue

        # Check similarity to existing centroids
        assigned = False
        for cluster_id, (_centroid_idx, centroid_fp) in enumerate(cluster_centroids):
            sim = sim_func(fps[i], centroid_fp)
            if sim >= threshold:
                cluster_ids[i] = cluster_id
                assigned = True
                break

        # Create new cluster if not assigned
        if not assigned:
            cluster_ids[i] = len(cluster_centroids)
            cluster_centroids.append((i, fps[i]))

    # Build results
    results = []
    for i in range(n):
        compound_id = id_list[i] if id_list else f"Compound_{i+1}"
        results.append(
            {
                "Compound_ID": compound_id,
                "SMILES": smiles_list[i],
                "Cluster_ID": cluster_ids[i],
                "Is_Centroid": any(idx == i for idx, _ in cluster_centroids),
            }
        )

    return pd.DataFrame(results)


def find_nearest_neighbors(
    query: Union[str, "Chem.Mol"],
    library: Union[str, pd.DataFrame, List[str]],
    k: int = 5,
    fp_type: str = "morgan",
    metric: str = "tanimoto",
    smiles_col: str = "SMILES",
    id_col: Optional[str] = "Compound_ID",
    fp_params: Optional[Dict] = None,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """
    Find k nearest neighbors to a query compound.

    Convenience wrapper around similarity_search for KNN queries.

    Parameters
    ----------
    query : str or Mol
        Query SMILES or Mol
    library : str, DataFrame, or list
        Library to search
    k : int
        Number of neighbors to return
    fp_type : str
        Fingerprint type
    metric : str
        Similarity metric
    smiles_col : str
        Column name for SMILES
    id_col : str, optional
        Column name for compound IDs
    fp_params : dict, optional
        Override fingerprint parameters
    n_jobs : int
        Number of parallel jobs

    Returns
    -------
    DataFrame
        K nearest neighbors
    """
    results = similarity_search(
        query=query,
        library=library,
        fp_type=fp_type,
        metric=metric,
        threshold=0.0,
        top_n=k,
        smiles_col=smiles_col,
        id_col=id_col,
        fp_params=fp_params,
        n_jobs=n_jobs,
    )

    return results.to_dataframe()


def _load_library(
    library: Union[str, pd.DataFrame, List[str]],
    smiles_col: str,
    id_col: Optional[str],
) -> Tuple[List[str], Optional[List[str]]]:
    """
    Load library from various input formats.

    Returns
    -------
    tuple
        (smiles_list, id_list) where id_list may be None
    """
    if isinstance(library, list):
        return library, None

    if isinstance(library, pd.DataFrame):
        df = library
    else:
        # Load from file
        path = str(library)
        if path.endswith(".csv"):
            df = pd.read_csv(path)
        elif path.endswith(".parquet"):
            df = pd.read_parquet(path)
        else:
            # Assume SMILES file
            df = _load_smiles_file(path)

    # Extract SMILES column.
    # If the caller specifies --smiles-col, prefer that, but still fall back
    # to the toolkit's "structure-of-record" priority chain.
    if smiles_col and smiles_col in df.columns:
        smiles_list = df[smiles_col].tolist()
    else:
        # MolProp structure-of-record priority (matches Quickstart docs)
        priority = [
            "Calc_Canonical_SMILES",
            "Calc_Base_SMILES",
            "Canonical_SMILES",
            "Input_Canonical_SMILES",
            "SMILES",
            "smiles",
        ]
        for col in priority:
            if col in df.columns:
                smiles_list = df[col].tolist()
                break
        else:
            # Try first column that looks like SMILES
            for col in df.columns:
                if df[col].dtype == object:
                    sample = str(df[col].iloc[0])
                    if any(c in sample for c in "CNOcnoPSFClBrI()=#"):
                        smiles_list = df[col].tolist()
                        break
            else:
                raise ValueError("Could not find SMILES column in data")

    # Extract ID column
    id_list = None
    if id_col is not None and id_col in df.columns:
        id_list = df[id_col].astype(str).tolist()
    elif "Compound_ID" in df.columns:
        id_list = df["Compound_ID"].astype(str).tolist()
    elif "Name" in df.columns:
        id_list = df["Name"].astype(str).tolist()
    elif "name" in df.columns:
        id_list = df["name"].astype(str).tolist()

    return smiles_list, id_list


def _load_smiles_file(path: str) -> pd.DataFrame:
    """Load a SMILES file into DataFrame."""
    records = []

    with open(path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Split on whitespace or tab
            parts = line.split()

            if len(parts) == 1:
                # SMILES only
                records.append({"SMILES": parts[0], "Name": f"Compound_{line_num}"})
            elif len(parts) >= 2:
                # Detect which is SMILES (contains typical SMILES characters)
                if any(c in parts[0] for c in "CNOcnoPSFClBrI()=#"):
                    records.append({"SMILES": parts[0], "Name": parts[1]})
                else:
                    records.append({"Name": parts[0], "SMILES": parts[1]})

    df = pd.DataFrame(records)
    df = df.rename(columns={"Name": "Compound_ID"})
    return df
