"""
Similarity metrics for molecular fingerprint comparison.

Provides implementations of common similarity coefficients used in
cheminformatics for comparing molecular fingerprints.
"""

from typing import Union, Optional
import numpy as np

try:
    from rdkit import DataStructs
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# Available similarity metrics
SIMILARITY_METRICS = {
    "tanimoto": "Tanimoto coefficient (Jaccard index) - most common for fingerprints",
    "dice": "Dice coefficient - weights matches more heavily",
    "cosine": "Cosine similarity - treats fingerprints as vectors",
    "sokal": "Sokal-Michener coefficient - includes both matches",
    "russel": "Russel-Rao coefficient - considers only shared bits",
    "kulczynski": "Kulczynski coefficient - asymmetric measure",
    "mcconnaughey": "McConnaughey coefficient - balanced measure",
    "asymmetric": "Asymmetric similarity (Tversky a=1, b=0)",
    "braun_blanquet": "Braun-Blanquet coefficient",
}


def _to_numpy(fp) -> np.ndarray:
    """Convert fingerprint to numpy array if needed."""
    if isinstance(fp, np.ndarray):
        return fp.astype(np.float64)
    
    if HAS_RDKIT:
        arr = np.zeros((fp.GetNumBits(),), dtype=np.float64)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    
    raise TypeError(f"Cannot convert fingerprint of type {type(fp)} to numpy array")


def tanimoto_similarity(fp1, fp2) -> float:
    """
    Calculate Tanimoto coefficient (Jaccard index) between two fingerprints.
    
    Tc = c / (a + b - c)
    
    Where:
    - a = bits on in fp1
    - b = bits on in fp2  
    - c = bits on in both
    
    Parameters
    ----------
    fp1, fp2 : BitVect or ndarray
        Molecular fingerprints to compare
        
    Returns
    -------
    float
        Tanimoto coefficient in range [0, 1]
    """
    # Use RDKit's optimized implementation if available
    if HAS_RDKIT and not isinstance(fp1, np.ndarray) and not isinstance(fp2, np.ndarray):
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    
    # Numpy implementation
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    intersection = np.sum(np.minimum(fp1, fp2))
    union = np.sum(np.maximum(fp1, fp2))
    
    if union == 0:
        return 0.0
    
    return float(intersection / union)


def dice_similarity(fp1, fp2) -> float:
    """
    Calculate Dice coefficient between two fingerprints.
    
    Dc = 2c / (a + b)
    
    Parameters
    ----------
    fp1, fp2 : BitVect or ndarray
        Molecular fingerprints to compare
        
    Returns
    -------
    float
        Dice coefficient in range [0, 1]
    """
    if HAS_RDKIT and not isinstance(fp1, np.ndarray) and not isinstance(fp2, np.ndarray):
        return DataStructs.DiceSimilarity(fp1, fp2)
    
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    intersection = np.sum(np.minimum(fp1, fp2))
    total = np.sum(fp1) + np.sum(fp2)
    
    if total == 0:
        return 0.0
    
    return float(2 * intersection / total)


def cosine_similarity(fp1, fp2) -> float:
    """
    Calculate cosine similarity between two fingerprints.
    
    cos = (fp1 · fp2) / (||fp1|| × ||fp2||)
    
    Parameters
    ----------
    fp1, fp2 : BitVect or ndarray
        Molecular fingerprints to compare
        
    Returns
    -------
    float
        Cosine similarity in range [0, 1] for binary fingerprints
    """
    if HAS_RDKIT and not isinstance(fp1, np.ndarray) and not isinstance(fp2, np.ndarray):
        return DataStructs.CosineSimilarity(fp1, fp2)
    
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    dot_product = np.dot(fp1, fp2)
    norm1 = np.linalg.norm(fp1)
    norm2 = np.linalg.norm(fp2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    return float(dot_product / (norm1 * norm2))


def sokal_similarity(fp1, fp2) -> float:
    """
    Calculate Sokal-Michener coefficient.
    
    Considers both shared presence and shared absence of bits.
    
    SM = (c + d) / (a + b + c + d)
    
    Where d = bits off in both
    """
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    both_on = np.sum((fp1 > 0) & (fp2 > 0))
    both_off = np.sum((fp1 == 0) & (fp2 == 0))
    total = len(fp1)
    
    return float((both_on + both_off) / total)


def russel_similarity(fp1, fp2) -> float:
    """
    Calculate Russel-Rao coefficient.
    
    RR = c / n
    
    Only considers bits on in both, relative to total bits.
    """
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    both_on = np.sum((fp1 > 0) & (fp2 > 0))
    total = len(fp1)
    
    return float(both_on / total)


def kulczynski_similarity(fp1, fp2) -> float:
    """
    Calculate Kulczynski coefficient.
    
    K = 0.5 × (c/a + c/b)
    
    Average of conditional probabilities.
    """
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    c = np.sum((fp1 > 0) & (fp2 > 0))
    a = np.sum(fp1 > 0)
    b = np.sum(fp2 > 0)
    
    if a == 0 or b == 0:
        return 0.0
    
    return float(0.5 * (c / a + c / b))


def mcconnaughey_similarity(fp1, fp2) -> float:
    """
    Calculate McConnaughey coefficient.
    
    Mc = (c² - (a-c)(b-c)) / (a × b)
    
    A balanced measure that can be negative.
    """
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    c = np.sum((fp1 > 0) & (fp2 > 0))
    a = np.sum(fp1 > 0)
    b = np.sum(fp2 > 0)
    
    if a == 0 or b == 0:
        return 0.0
    
    return float((c * c - (a - c) * (b - c)) / (a * b))


def asymmetric_similarity(fp1, fp2, alpha: float = 1.0, beta: float = 0.0) -> float:
    """
    Calculate Tversky index (asymmetric similarity).
    
    Tv = c / (c + α(a-c) + β(b-c))
    
    Parameters
    ----------
    fp1 : query fingerprint
    fp2 : target fingerprint
    alpha : weight for query-only bits (default=1.0)
    beta : weight for target-only bits (default=0.0)
    
    With alpha=beta=1, equals Tanimoto.
    With alpha=beta=0.5, equals Dice.
    With alpha=1, beta=0, measures substructure similarity.
    """
    if HAS_RDKIT and not isinstance(fp1, np.ndarray) and not isinstance(fp2, np.ndarray):
        return DataStructs.TverskySimilarity(fp1, fp2, alpha, beta)
    
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    c = np.sum((fp1 > 0) & (fp2 > 0))
    a_only = np.sum((fp1 > 0) & (fp2 == 0))
    b_only = np.sum((fp1 == 0) & (fp2 > 0))
    
    denominator = c + alpha * a_only + beta * b_only
    
    if denominator == 0:
        return 0.0
    
    return float(c / denominator)


def braun_blanquet_similarity(fp1, fp2) -> float:
    """
    Calculate Braun-Blanquet coefficient.
    
    BB = c / max(a, b)
    """
    fp1 = _to_numpy(fp1)
    fp2 = _to_numpy(fp2)
    
    c = np.sum((fp1 > 0) & (fp2 > 0))
    a = np.sum(fp1 > 0)
    b = np.sum(fp2 > 0)
    
    if max(a, b) == 0:
        return 0.0
    
    return float(c / max(a, b))


def get_similarity_function(metric: str):
    """
    Get similarity function by name.
    
    Parameters
    ----------
    metric : str
        Name of similarity metric
        
    Returns
    -------
    callable
        Similarity function
    """
    metric = metric.lower()
    
    functions = {
        "tanimoto": tanimoto_similarity,
        "dice": dice_similarity,
        "cosine": cosine_similarity,
        "sokal": sokal_similarity,
        "russel": russel_similarity,
        "kulczynski": kulczynski_similarity,
        "mcconnaughey": mcconnaughey_similarity,
        "asymmetric": asymmetric_similarity,
        "braun_blanquet": braun_blanquet_similarity,
    }
    
    if metric not in functions:
        raise ValueError(
            f"Unknown similarity metric: {metric}. "
            f"Available: {list(functions.keys())}"
        )
    
    return functions[metric]


def bulk_similarity(query_fp, target_fps, metric: str = "tanimoto") -> np.ndarray:
    """
    Calculate similarity between a query and multiple targets.
    
    Parameters
    ----------
    query_fp : BitVect or ndarray
        Query fingerprint
    target_fps : list
        List of target fingerprints
    metric : str
        Similarity metric to use
        
    Returns
    -------
    ndarray
        Array of similarity values
    """
    sim_func = get_similarity_function(metric)
    
    # Use RDKit bulk calculation if available and applicable
    if (HAS_RDKIT and metric == "tanimoto" and 
        not isinstance(query_fp, np.ndarray) and
        all(not isinstance(fp, np.ndarray) for fp in target_fps if fp is not None)):
        
        similarities = []
        for fp in target_fps:
            if fp is None:
                similarities.append(np.nan)
            else:
                similarities.append(DataStructs.TanimotoSimilarity(query_fp, fp))
        return np.array(similarities)
    
    # General implementation
    similarities = []
    for fp in target_fps:
        if fp is None:
            similarities.append(np.nan)
        else:
            similarities.append(sim_func(query_fp, fp))
    
    return np.array(similarities)
