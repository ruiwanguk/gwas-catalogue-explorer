"""Tests for simplified Mendelian Randomization module."""

import numpy as np
import pandas as pd
import responses

from gwas_explorer.config import OPENGWAS_API_URL
from gwas_explorer.mr_analysis import run_mr_analysis


@responses.activate
def test_wald_ratio_calculation(sample_candidate_genes: pd.DataFrame) -> None:
    """Single-SNP MR should compute correct Wald ratio."""
    # Mock: get instruments for gene (SNP-exposure)
    responses.add(
        responses.GET,
        f"{OPENGWAS_API_URL}/associations",
        json=[
            {
                "rsid": "rs5219",
                "beta": 0.3,
                "se": 0.05,
                "p": 1e-10,
                "ea": "T",
                "nea": "C",
                "eaf": 0.35,
            },
        ],
        status=200,
    )
    # Mock: get outcome associations (SNP-T2D)
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=[
            {
                "rsid": "rs5219",
                "beta": 0.15,
                "se": 0.03,
                "p": 1e-6,
                "ea": "T",
                "nea": "C",
                "eaf": 0.35,
            },
        ],
        status=200,
    )

    # Test with just KCNJ11
    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    # Wald ratio = beta_outcome / beta_exposure = 0.15 / 0.3 = 0.5
    assert np.isclose(row["MR_ESTIMATE"], 0.5, atol=0.01)
    assert row["METHOD"] == "wald_ratio"
    assert row["MR_STATUS"] == "ok"


@responses.activate
def test_insufficient_instruments(sample_candidate_genes: pd.DataFrame) -> None:
    """Should flag genes with no valid instruments."""
    # Mock: no instruments found
    responses.add(
        responses.GET,
        f"{OPENGWAS_API_URL}/associations",
        json=[],
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "TCF7L2"].copy()
    result = run_mr_analysis(genes)

    assert result.iloc[0]["MR_STATUS"] == "insufficient_data"
    assert np.isnan(result.iloc[0]["MR_ESTIMATE"])


@responses.activate
def test_handles_api_error(sample_candidate_genes: pd.DataFrame) -> None:
    """Should not crash on API failures."""
    responses.add(
        responses.GET,
        f"{OPENGWAS_API_URL}/associations",
        status=500,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "PPARG"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    assert result.iloc[0]["MR_STATUS"] == "error"
