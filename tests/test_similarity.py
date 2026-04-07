"""
Tests for molprop_toolkit.similarity module.

Run with: pytest tests/test_similarity.py -v
"""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

# Skip all tests if RDKit not available
rdkit = pytest.importorskip("rdkit")

from molprop_toolkit.similarity.fingerprints import (
    FINGERPRINT_TYPES,
    fingerprint_density,
    get_fingerprint,
    get_fingerprint_bulk,
    get_on_bits,
)
from molprop_toolkit.similarity.metrics import (
    SIMILARITY_METRICS,
    bulk_similarity,
    cosine_similarity,
    dice_similarity,
    get_similarity_function,
    tanimoto_similarity,
)
from molprop_toolkit.similarity.search import (
    cluster_by_similarity,
    diversity_pick,
    find_nearest_neighbors,
    pairwise_similarity,
    similarity_search,
)

# Test data
TEST_SMILES = [
    "CCO",  # ethanol
    "CCCO",  # propanol
    "CCCCO",  # butanol
    "c1ccccc1",  # benzene
    "c1ccccc1O",  # phenol
    "c1ccccc1N",  # aniline
    "CC(=O)O",  # acetic acid
    "CC(=O)N",  # acetamide
    "CCCCCCCC",  # octane
    "C1CCCCC1",  # cyclohexane
]

INVALID_SMILES = ["not_a_smiles", "XXXXX", ""]


class TestFingerprints:
    """Tests for fingerprint generation."""

    def test_get_fingerprint_morgan(self):
        """Test Morgan fingerprint generation."""
        fp = get_fingerprint("CCO", fp_type="morgan")
        assert fp is not None
        assert fp.GetNumBits() == 2048
        assert fp.GetNumOnBits() > 0

    def test_get_fingerprint_maccs(self):
        """Test MACCS fingerprint generation."""
        fp = get_fingerprint("CCO", fp_type="maccs")
        assert fp is not None
        assert fp.GetNumBits() == 167  # MACCS has 166 keys + 1

    @pytest.mark.parametrize("fp_type", FINGERPRINT_TYPES.keys())
    def test_all_fingerprint_types(self, fp_type):
        """Test all fingerprint types generate valid output."""
        fp = get_fingerprint("c1ccccc1", fp_type=fp_type)
        assert fp is not None
        assert fp.GetNumOnBits() > 0

    def test_get_fingerprint_as_numpy(self):
        """Test conversion to numpy array."""
        fp = get_fingerprint("CCO", fp_type="morgan", as_numpy=True)
        assert isinstance(fp, np.ndarray)
        assert fp.shape == (2048,)
        assert fp.sum() > 0

    def test_get_fingerprint_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        fp = get_fingerprint("not_a_smiles")
        assert fp is None

    def test_get_fingerprint_custom_params(self):
        """Test custom fingerprint parameters.

        Use a structure where radius affects the neighborhood environment.
        Very small molecules (e.g., ethanol) may yield identical fingerprints
        for different radii.
        """
        smi = "CCOc1ccccc1"  # phenetole
        fp1 = get_fingerprint(smi, fp_type="morgan", radius=2, nBits=1024)
        fp2 = get_fingerprint(smi, fp_type="morgan", radius=3, nBits=1024)

        assert fp1.GetNumBits() == 1024
        assert fp2.GetNumBits() == 1024
        # Different radius should give different fingerprints
        assert list(fp1.GetOnBits()) != list(fp2.GetOnBits())

    def test_get_fingerprint_bulk(self):
        """Test bulk fingerprint generation."""
        fps = get_fingerprint_bulk(TEST_SMILES, fp_type="morgan")

        assert len(fps) == len(TEST_SMILES)
        assert all(fp is not None for fp in fps)

    def test_get_fingerprint_bulk_with_invalid(self):
        """Test bulk generation with invalid SMILES."""
        smiles = TEST_SMILES + INVALID_SMILES
        fps = get_fingerprint_bulk(smiles, fp_type="morgan")

        assert len(fps) == len(smiles)
        # Valid SMILES should have fingerprints
        assert all(fps[i] is not None for i in range(len(TEST_SMILES)))
        # Invalid SMILES should be None
        assert all(fps[i] is None for i in range(len(TEST_SMILES), len(smiles)))

    def test_fingerprint_density(self):
        """Test fingerprint density calculation."""
        fp = get_fingerprint("CCO", fp_type="morgan")
        density = fingerprint_density(fp)

        assert 0 <= density <= 1
        assert density == fp.GetNumOnBits() / fp.GetNumBits()

    def test_get_on_bits(self):
        """Test getting on bit indices."""
        fp = get_fingerprint("CCO", fp_type="morgan")
        on_bits = get_on_bits(fp)

        assert isinstance(on_bits, list)
        assert len(on_bits) == fp.GetNumOnBits()
        assert all(isinstance(b, int) for b in on_bits)


class TestMetrics:
    """Tests for similarity metrics."""

    @pytest.fixture
    def fp_pair(self):
        """Generate a pair of fingerprints for testing."""
        fp1 = get_fingerprint("CCO")
        fp2 = get_fingerprint("CCCO")
        return fp1, fp2

    def test_tanimoto_similarity(self, fp_pair):
        """Test Tanimoto coefficient."""
        fp1, fp2 = fp_pair
        sim = tanimoto_similarity(fp1, fp2)

        assert 0 <= sim <= 1
        # Self-similarity should be 1
        assert tanimoto_similarity(fp1, fp1) == pytest.approx(1.0)

    def test_dice_similarity(self, fp_pair):
        """Test Dice coefficient."""
        fp1, fp2 = fp_pair
        sim = dice_similarity(fp1, fp2)

        assert 0 <= sim <= 1
        assert dice_similarity(fp1, fp1) == pytest.approx(1.0)
        # Dice >= Tanimoto for same pair
        assert sim >= tanimoto_similarity(fp1, fp2)

    def test_cosine_similarity(self, fp_pair):
        """Test cosine similarity."""
        fp1, fp2 = fp_pair
        sim = cosine_similarity(fp1, fp2)

        assert 0 <= sim <= 1
        assert cosine_similarity(fp1, fp1) == pytest.approx(1.0)

    @pytest.mark.parametrize("metric", SIMILARITY_METRICS.keys())
    def test_all_metrics(self, fp_pair, metric):
        """Test all similarity metrics."""
        fp1, fp2 = fp_pair
        sim_func = get_similarity_function(metric)
        sim = sim_func(fp1, fp2)

        # Most metrics should be bounded, but McConnaughey can be negative
        if metric != "mcconnaughey":
            assert 0 <= sim <= 1
        else:
            assert -1 <= sim <= 1

    def test_numpy_compatibility(self, fp_pair):
        """Test that metrics work with numpy arrays."""
        fp1_np = get_fingerprint("CCO", as_numpy=True)
        fp2_np = get_fingerprint("CCCO", as_numpy=True)
        fp1, fp2 = fp_pair

        sim_rdkit = tanimoto_similarity(fp1, fp2)
        sim_numpy = tanimoto_similarity(fp1_np, fp2_np)

        assert sim_rdkit == pytest.approx(sim_numpy, rel=1e-5)

    def test_bulk_similarity(self):
        """Test bulk similarity calculation."""
        query = get_fingerprint("CCO")
        targets = get_fingerprint_bulk(TEST_SMILES)

        sims = bulk_similarity(query, targets, metric="tanimoto")

        assert len(sims) == len(TEST_SMILES)
        assert all(0 <= s <= 1 for s in sims)
        # First compound is same as query (CCO)
        assert sims[0] == pytest.approx(1.0)

    def test_get_similarity_function_invalid(self):
        """Test error handling for invalid metric."""
        with pytest.raises(ValueError, match="Unknown similarity metric"):
            get_similarity_function("invalid_metric")


class TestSearch:
    """Tests for similarity search functions."""

    @pytest.fixture
    def library_file(self):
        """Create a temporary library file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            for i, smi in enumerate(TEST_SMILES):
                f.write(f"{smi}\tCompound_{i+1}\n")
            f.flush()
            yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def library_csv(self):
        """Create a temporary CSV library."""
        df = pd.DataFrame(
            {
                "Compound_ID": [f"Cpd_{i+1}" for i in range(len(TEST_SMILES))],
                "SMILES": TEST_SMILES,
            }
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            df.to_csv(f.name, index=False)
            yield f.name
        os.unlink(f.name)

    def test_similarity_search_basic(self, library_file):
        """Test basic similarity search."""
        results = similarity_search(
            query="CCO",
            library=library_file,
            fp_type="morgan",
            top_n=5,
        )

        assert len(results) <= 5
        assert results.n_searched == len(TEST_SMILES)
        assert results.n_valid == len(TEST_SMILES)

        # Results should be sorted by similarity (descending)
        sims = [r.similarity for r in results]
        assert sims == sorted(sims, reverse=True)

    def test_similarity_search_threshold(self, library_file):
        """Test similarity search with threshold."""
        results = similarity_search(
            query="CCO",
            library=library_file,
            threshold=0.5,
        )

        # All results should meet threshold
        assert all(r.similarity >= 0.5 for r in results)

    def test_similarity_search_from_list(self):
        """Test search with list input."""
        results = similarity_search(
            query="CCO",
            library=TEST_SMILES,
            fp_type="morgan",
        )

        assert len(results) > 0
        assert results.n_searched == len(TEST_SMILES)

    def test_similarity_search_from_dataframe(self):
        """Test search with DataFrame input."""
        df = pd.DataFrame(
            {
                "SMILES": TEST_SMILES,
                "Compound_ID": [f"ID_{i}" for i in range(len(TEST_SMILES))],
            }
        )

        results = similarity_search(
            query="c1ccccc1",
            library=df,
            fp_type="morgan",
        )

        assert len(results) > 0
        # Check that IDs are preserved
        assert all("ID_" in r.compound_id for r in results)

    def test_similarity_search_to_dataframe(self, library_file):
        """Test converting results to DataFrame."""
        results = similarity_search(
            query="CCO",
            library=library_file,
            top_n=5,
        )

        df = results.to_dataframe()

        assert isinstance(df, pd.DataFrame)
        assert "Rank" in df.columns
        assert "Compound_ID" in df.columns
        assert "SMILES" in df.columns
        assert "Similarity" in df.columns
        assert len(df) == len(results)

    def test_similarity_search_to_csv(self, library_file):
        """Test saving results to CSV."""
        results = similarity_search(
            query="CCO",
            library=library_file,
            top_n=5,
        )

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            results.to_csv(f.name)
            df = pd.read_csv(f.name)

            assert len(df) == len(results)
            os.unlink(f.name)

    def test_similarity_search_csv_library(self, library_csv):
        """Test search with CSV library."""
        results = similarity_search(
            query="CCO",
            library=library_csv,
            top_n=5,
        )

        assert len(results) > 0
        assert "Cpd_" in results.results[0].compound_id

    def test_pairwise_similarity(self):
        """Test pairwise similarity matrix."""
        smiles = TEST_SMILES[:5]
        matrix, smiles_out = pairwise_similarity(smiles, fp_type="morgan")

        assert matrix.shape == (5, 5)
        assert smiles_out == smiles

        # Diagonal should be 1
        assert np.allclose(np.diag(matrix), 1.0)

        # Matrix should be symmetric
        assert np.allclose(matrix, matrix.T)

    def test_diversity_pick(self, library_file):
        """Test diversity picking."""
        df = diversity_pick(
            molecules=library_file,
            n_pick=5,
            fp_type="morgan",
        )

        assert len(df) == 5
        assert "Compound_ID" in df.columns
        assert "SMILES" in df.columns
        assert "Pick_Order" in df.columns
        assert "Min_Distance" in df.columns

        # Pick order should be 1, 2, 3, ...
        assert list(df["Pick_Order"]) == [1, 2, 3, 4, 5]

    def test_diversity_pick_with_seed(self):
        """Test diversity picking with seed compound."""
        df = diversity_pick(
            molecules=TEST_SMILES,
            n_pick=3,
            seed="c1ccccc1",
            fp_type="morgan",
        )

        # First pick should be the seed
        assert df.iloc[0]["SMILES"] == "c1ccccc1"

    def test_cluster_by_similarity(self, library_file):
        """Test similarity-based clustering."""
        df = cluster_by_similarity(
            molecules=library_file,
            threshold=0.5,
            fp_type="morgan",
        )

        assert len(df) == len(TEST_SMILES)
        assert "Cluster_ID" in df.columns
        assert "Is_Centroid" in df.columns

        # Should have at least one cluster
        assert df["Cluster_ID"].max() >= 0

        # Should have at least one centroid
        assert df["Is_Centroid"].sum() >= 1

    def test_find_nearest_neighbors(self, library_file):
        """Test KNN search."""
        df = find_nearest_neighbors(
            query="CCO",
            library=library_file,
            k=3,
            fp_type="morgan",
        )

        assert len(df) == 3
        assert "Similarity" in df.columns


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_library(self):
        """Test search with empty library."""
        results = similarity_search(
            query="CCO",
            library=[],
            fp_type="morgan",
        )

        assert len(results) == 0

    def test_single_compound_library(self):
        """Test search with single compound."""
        results = similarity_search(
            query="CCO",
            library=["CCO"],
            fp_type="morgan",
        )

        assert len(results) == 1
        assert results.results[0].similarity == pytest.approx(1.0)

    def test_invalid_query(self):
        """Test with invalid query SMILES."""
        with pytest.raises(ValueError, match="Could not generate fingerprint"):
            similarity_search(
                query="not_valid_smiles",
                library=TEST_SMILES,
            )

    def test_diversity_pick_more_than_available(self):
        """Test diversity pick requesting more than available."""
        with pytest.raises(ValueError, match="Cannot pick"):
            diversity_pick(
                molecules=TEST_SMILES[:3],
                n_pick=10,
            )

    def test_unknown_fingerprint_type(self):
        """Test with unknown fingerprint type."""
        with pytest.raises(ValueError, match="Unknown fingerprint type"):
            get_fingerprint("CCO", fp_type="unknown_fp")


class TestIntegration:
    """Integration tests for complete workflows."""

    def test_full_workflow(self):
        """Test complete search workflow."""
        # Generate library
        library = TEST_SMILES

        # Search
        results = similarity_search(
            query="c1ccccc1O",  # phenol
            library=library,
            fp_type="morgan",
            threshold=0.3,
            top_n=5,
        )

        # Convert to DataFrame
        df = results.to_dataframe()

        # Verify results make sense
        assert len(df) > 0
        # Phenol should be similar to benzene and aniline
        smiles_found = set(df["SMILES"])
        assert "c1ccccc1O" in smiles_found  # phenol itself

    def test_diversity_followed_by_search(self):
        """Test diversity picking followed by search."""
        # Pick diverse subset
        diverse = diversity_pick(
            molecules=TEST_SMILES,
            n_pick=5,
            fp_type="morgan",
        )

        # Use diverse set as library for search
        diverse_smiles = diverse["SMILES"].tolist()

        results = similarity_search(
            query="CCCCO",
            library=diverse_smiles,
            fp_type="morgan",
        )

        assert len(results) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
