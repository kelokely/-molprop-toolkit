"""
Fingerprint generation utilities for MolProp Toolkit.

Supports multiple fingerprint types commonly used in medicinal chemistry:
- Morgan (ECFP-like circular fingerprints)
- MACCS (166-bit structural keys)
- RDKit (RDKit topological fingerprint)
- Atom Pair
- Topological Torsion
- Pattern (SMARTS-based substructure)
"""

from typing import List, Optional, Union

import numpy as np

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# Fingerprint type configurations
FINGERPRINT_TYPES = {
    "morgan": {
        "description": "Morgan circular fingerprint (ECFP-like)",
        "default_params": {
            "radius": 2,
            "nBits": 2048,
            "useFeatures": False,
            "useChirality": False,
        },
    },
    "morgan_feat": {
        "description": "Morgan fingerprint with pharmacophoric features (FCFP-like)",
        "default_params": {
            "radius": 2,
            "nBits": 2048,
            "useFeatures": True,
            "useChirality": False,
        },
    },
    "maccs": {
        "description": "MACCS 166-bit structural keys",
        "default_params": {},
    },
    "rdkit": {
        "description": "RDKit topological fingerprint",
        "default_params": {"minPath": 1, "maxPath": 7, "fpSize": 2048},
    },
    "atompair": {
        "description": "Atom pair fingerprint",
        "default_params": {"nBits": 2048},
    },
    "torsion": {
        "description": "Topological torsion fingerprint",
        "default_params": {"nBits": 2048},
    },
    "pattern": {
        "description": "Pattern fingerprint (substructure keys)",
        "default_params": {"fpSize": 2048},
    },
}


def _check_rdkit():
    """Raise ImportError if RDKit is not available."""
    if not HAS_RDKIT:
        raise ImportError(
            "RDKit is required for fingerprint generation. "
            "Install with: conda install -c conda-forge rdkit"
        )


def _mol_from_input(mol_input: Union[str, "Chem.Mol"]) -> Optional["Chem.Mol"]:
    """
    Convert input to RDKit Mol object.

    Parameters
    ----------
    mol_input : str or Mol
        SMILES string or RDKit Mol object

    Returns
    -------
    Mol or None
        RDKit Mol object, or None if conversion fails
    """
    _check_rdkit()

    if mol_input is None:
        return None
    if isinstance(mol_input, Chem.Mol):
        return mol_input
    if isinstance(mol_input, str):
        s = mol_input.strip()
        if not s:
            return None
        mol = Chem.MolFromSmiles(s)
        return mol
    return None


def get_fingerprint(
    mol_input: Union[str, "Chem.Mol"],
    fp_type: str = "morgan",
    as_numpy: bool = False,
    **kwargs,
) -> Optional[Union["DataStructs.ExplicitBitVect", np.ndarray]]:
    """
    Generate a molecular fingerprint.

    Parameters
    ----------
    mol_input : str or Mol
        SMILES string or RDKit Mol object
    fp_type : str
        Fingerprint type. One of: morgan, morgan_feat, maccs, rdkit,
        atompair, torsion, pattern
    as_numpy : bool
        If True, return as numpy array instead of RDKit BitVect
    **kwargs
        Override default fingerprint parameters

    Returns
    -------
    BitVect or ndarray or None
        Fingerprint object, or None if molecule is invalid

    Examples
    --------
    >>> fp = get_fingerprint("CCO", fp_type="morgan", radius=2)
    >>> fp_array = get_fingerprint("CCO", fp_type="maccs", as_numpy=True)
    """
    _check_rdkit()

    mol = _mol_from_input(mol_input)
    if mol is None:
        return None

    fp_type = fp_type.lower()
    if fp_type not in FINGERPRINT_TYPES:
        raise ValueError(
            f"Unknown fingerprint type: {fp_type}. "
            f"Available: {list(FINGERPRINT_TYPES.keys())}"
        )

    # Merge default params with user overrides
    params = FINGERPRINT_TYPES[fp_type]["default_params"].copy()
    params.update(kwargs)

    # Generate fingerprint based on type
    fp = None

    if fp_type in ("morgan", "morgan_feat"):
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            radius=int(params.get("radius", 2)),
            nBits=int(params.get("nBits", 2048)),
            useFeatures=bool(params.get("useFeatures", fp_type == "morgan_feat")),
            useChirality=bool(params.get("useChirality", False)),
        )

    elif fp_type == "maccs":
        fp = MACCSkeys.GenMACCSKeys(mol)

    elif fp_type == "rdkit":
        fp = Chem.RDKFingerprint(
            mol,
            minPath=params.get("minPath", 1),
            maxPath=params.get("maxPath", 7),
            fpSize=params.get("fpSize", 2048),
        )

    elif fp_type == "atompair":
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol,
            nBits=params.get("nBits", 2048),
        )

    elif fp_type == "torsion":
        fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol,
            nBits=params.get("nBits", 2048),
        )

    elif fp_type == "pattern":
        fp = Chem.PatternFingerprint(
            mol,
            fpSize=params.get("fpSize", 2048),
        )

    if fp is None:
        return None

    if as_numpy:
        arr = np.zeros((fp.GetNumBits(),), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    return fp


def get_fingerprint_bulk(
    mol_inputs: List[Union[str, "Chem.Mol"]],
    fp_type: str = "morgan",
    as_numpy: bool = False,
    n_jobs: int = 1,
    show_progress: bool = False,
    **kwargs,
) -> List[Optional[Union["DataStructs.ExplicitBitVect", np.ndarray]]]:
    """
    Generate fingerprints for multiple molecules.

    Parameters
    ----------
    mol_inputs : list
        List of SMILES strings or RDKit Mol objects
    fp_type : str
        Fingerprint type
    as_numpy : bool
        If True, return as numpy arrays
    n_jobs : int
        Number of parallel jobs (requires joblib). Use -1 for all CPUs.
    show_progress : bool
        Show progress bar (requires tqdm)
    **kwargs
        Override default fingerprint parameters

    Returns
    -------
    list
        List of fingerprints (None for invalid molecules)
    """
    _check_rdkit()

    if n_jobs != 1:
        try:
            from joblib import Parallel, delayed

            use_parallel = True
        except ImportError:
            use_parallel = False
            if n_jobs != 1:
                import warnings

                warnings.warn(
                    "joblib not installed, falling back to serial processing. "
                    "Install with: pip install joblib",
                    stacklevel=2,
                )
    else:
        use_parallel = False

    # Setup progress bar if requested
    iterator = mol_inputs
    if show_progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(mol_inputs, desc="Generating fingerprints")
        except ImportError:
            pass

    if use_parallel and n_jobs != 1:
        fps = Parallel(n_jobs=n_jobs)(
            delayed(get_fingerprint)(mol, fp_type, as_numpy, **kwargs)
            for mol in iterator
        )
    else:
        fps = [get_fingerprint(mol, fp_type, as_numpy, **kwargs) for mol in iterator]

    return fps


def fingerprint_to_bitvect(fp: np.ndarray) -> "DataStructs.ExplicitBitVect":
    """Convert numpy array back to RDKit BitVect."""
    _check_rdkit()
    bitvect = DataStructs.ExplicitBitVect(len(fp))
    on_bits = np.where(fp > 0)[0].tolist()
    bitvect.SetBitsFromList(on_bits)
    return bitvect


def get_on_bits(fp) -> List[int]:
    """Get list of 'on' bit indices from a fingerprint."""
    _check_rdkit()
    if isinstance(fp, np.ndarray):
        return np.where(fp > 0)[0].tolist()
    return list(fp.GetOnBits())


def fingerprint_density(fp) -> float:
    """Calculate the bit density (fraction of bits set) of a fingerprint."""
    if isinstance(fp, np.ndarray):
        return np.sum(fp > 0) / len(fp)
    return fp.GetNumOnBits() / fp.GetNumBits()
