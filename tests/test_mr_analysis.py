"""Tests for simplified Mendelian Randomization module."""

import numpy as np
import pandas as pd
import responses

from gwas_explorer.config import OPENGWAS_API_URL
from gwas_explorer.mr_analysis import _harmonize_alleles, _ivw, run_mr_analysis


@responses.activate
def test_eqtl_instruments_from_tophits(sample_candidate_genes: pd.DataFrame) -> None:
    """Should fetch eQTL instruments via /tophits and outcome via /associations."""
    eqtl_instruments = [
        {
            "rsid": "rs10885396",
            "beta": -0.285,
            "se": 0.01,
            "p": 8.3e-132,
            "ea": "A",
            "nea": "G",
            "eaf": 0.4,
            "id": "eqtl-a-ENSG00000187486",
            "trait": "ENSG00000187486",
        },
    ]
    outcome_associations = [
        {
            "rsid": "rs10885396",
            "beta": 0.02,
            "se": 0.008,
            "p": 0.01,
            "ea": "A",
            "nea": "G",
            "eaf": 0.4,
            "id": "ebi-a-GCST006867",
        },
    ]

    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=eqtl_instruments,
        status=200,
    )
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=outcome_associations,
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["MR_STATUS"] == "ok"
    assert row["METHOD"] == "wald_ratio"
    assert row["N_INSTRUMENTS"] == 1
    # Wald ratio = 0.02 / -0.285 ≈ -0.0702
    assert abs(row["MR_ESTIMATE"] - (0.02 / -0.285)) < 0.01


@responses.activate
def test_wald_ratio_calculation(sample_candidate_genes: pd.DataFrame) -> None:
    """Single-SNP MR should compute correct Wald ratio."""
    instruments = [
        {"rsid": "rs5219", "beta": 0.3, "se": 0.05, "p": 1e-10, "ea": "T", "nea": "C", "eaf": 0.35},
    ]
    outcomes = [
        {"rsid": "rs5219", "beta": 0.15, "se": 0.03, "p": 1e-6, "ea": "T", "nea": "C", "eaf": 0.35},
    ]
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=instruments,
        status=200,
    )
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=outcomes,
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    # Wald ratio = beta_outcome / beta_exposure = 0.15 / 0.3 = 0.5
    assert np.isclose(row["MR_ESTIMATE"], 0.5, atol=0.01)
    assert row["METHOD"] == "wald_ratio"
    assert row["MR_STATUS"] == "ok"


@responses.activate
def test_ivw_calculation(sample_candidate_genes: pd.DataFrame) -> None:
    """Multi-SNP MR should compute correct IVW estimate."""
    instruments = [
        {"rsid": "rs5219", "beta": 0.3, "se": 0.05, "p": 1e-10, "ea": "T", "nea": "C", "eaf": 0.35},
        {"rsid": "rs5220", "beta": 0.2, "se": 0.04, "p": 1e-8, "ea": "A", "nea": "G", "eaf": 0.4},
    ]
    outcomes = [
        {"rsid": "rs5219", "beta": 0.15, "se": 0.03, "p": 1e-6, "ea": "T", "nea": "C", "eaf": 0.35},
        {
            "rsid": "rs5220",
            "beta": 0.08,
            "se": 0.025,
            "p": 0.001,
            "ea": "A",
            "nea": "G",
            "eaf": 0.4,
        },
    ]
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=instruments,
        status=200,
    )
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=outcomes,
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["METHOD"] == "ivw"
    assert row["MR_STATUS"] == "ok"
    assert row["N_INSTRUMENTS"] == 2
    # Hand-computed IVW ≈ 0.4610
    assert np.isclose(row["MR_ESTIMATE"], 0.4610, atol=0.01)


@responses.activate
def test_insufficient_instruments(sample_candidate_genes: pd.DataFrame) -> None:
    """Should flag genes with no valid instruments."""
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=[],
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "TCF7L2"].copy()
    result = run_mr_analysis(genes)

    assert result.iloc[0]["MR_STATUS"] == "no_eqtl_instruments"
    assert np.isnan(result.iloc[0]["MR_ESTIMATE"])


@responses.activate
def test_handles_api_error(sample_candidate_genes: pd.DataFrame) -> None:
    """Should not crash on API failures."""
    responses.add(responses.POST, f"{OPENGWAS_API_URL}/tophits", status=500)

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "PPARG"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    assert result.iloc[0]["MR_STATUS"] == "error"


@responses.activate
def test_wald_ratio_zero_beta(sample_candidate_genes: pd.DataFrame) -> None:
    """Should handle zero beta_exp without division by zero."""
    instruments = [
        {"rsid": "rs5219", "beta": 0.0, "se": 0.05, "p": 1e-10, "ea": "T", "nea": "C", "eaf": 0.35},
    ]
    outcomes = [
        {"rsid": "rs5219", "beta": 0.15, "se": 0.03, "p": 1e-6, "ea": "T", "nea": "C", "eaf": 0.35},
    ]
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=instruments,
        status=200,
    )
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=outcomes,
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    assert result.iloc[0]["MR_STATUS"] == "invalid_instrument"


@responses.activate
def test_allele_harmonization_flipped(sample_candidate_genes: pd.DataFrame) -> None:
    """Should flip outcome beta when alleles are swapped."""
    instruments = [
        {"rsid": "rs5219", "beta": 0.3, "se": 0.05, "p": 1e-10, "ea": "T", "nea": "C", "eaf": 0.35},
    ]
    # Outcome has swapped alleles: ea=C, nea=T (opposite of instrument)
    outcomes = [
        {
            "rsid": "rs5219",
            "beta": -0.15,
            "se": 0.03,
            "p": 1e-6,
            "ea": "C",
            "nea": "T",
            "eaf": 0.65,
        },
    ]
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=instruments,
        status=200,
    )
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=outcomes,
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    # After harmonization: beta_out flipped from -0.15 to 0.15
    # Wald ratio = 0.15 / 0.3 = 0.5
    assert np.isclose(row["MR_ESTIMATE"], 0.5, atol=0.01)
    assert row["MR_STATUS"] == "ok"


@responses.activate
def test_no_eqtl_data_for_gene(sample_candidate_genes: pd.DataFrame) -> None:
    """Genes without eQTLGen data should get no_eqtl_instruments status."""
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=[],
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "HLA-A"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    assert result.iloc[0]["MR_STATUS"] == "no_eqtl_instruments"


@responses.activate
def test_empty_candidates_returns_empty() -> None:
    """Should return empty DataFrame with correct columns for empty input."""
    empty = pd.DataFrame(columns=["GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP"])
    result = run_mr_analysis(empty)
    assert len(result) == 0
    expected = {
        "GENE_SYMBOL",
        "MR_ESTIMATE",
        "MR_SE",
        "MR_PVALUE",
        "N_INSTRUMENTS",
        "METHOD",
        "MR_STATUS",
    }
    assert set(result.columns) == expected


def test_palindromic_at_snp_excluded() -> None:
    """A/T palindromic SNP with mismatched outcome alleles should return None."""
    # Instrument is A/T (palindromic) and outcome reports unresolvable alleles
    inst = {"rsid": "rs123", "beta": 0.3, "se": 0.05, "ea": "A", "nea": "T"}
    out = {"rsid": "rs123", "beta": 0.1, "se": 0.02, "ea": "G", "nea": "C"}
    result = _harmonize_alleles(inst, out)
    assert result is None


def test_palindromic_cg_snp_excluded() -> None:
    """C/G palindromic SNP with mismatched outcome alleles should return None."""
    inst = {"rsid": "rs456", "beta": 0.2, "se": 0.04, "ea": "C", "nea": "G"}
    out = {"rsid": "rs456", "beta": 0.05, "se": 0.01, "ea": "A", "nea": "T"}
    result = _harmonize_alleles(inst, out)
    assert result is None


def test_harmonize_missing_allele_passes_through() -> None:
    """Missing allele info should pass through unchanged."""
    inst = {"rsid": "rs789", "beta": 0.3, "se": 0.05, "ea": "A", "nea": ""}
    out = {"rsid": "rs789", "beta": 0.1, "se": 0.02, "ea": "T", "nea": "C"}
    result = _harmonize_alleles(inst, out)
    assert result is out


def test_harmonize_aligned_alleles() -> None:
    """Already aligned alleles should return outcome unchanged."""
    inst = {"rsid": "rs100", "beta": 0.3, "se": 0.05, "ea": "A", "nea": "G"}
    out = {"rsid": "rs100", "beta": 0.1, "se": 0.02, "ea": "A", "nea": "G"}
    result = _harmonize_alleles(inst, out)
    assert result is out


@responses.activate
def test_ivw_excludes_palindromic_snps(sample_candidate_genes: pd.DataFrame) -> None:
    """IVW should skip palindromic SNPs when computing the estimate."""
    instruments = [
        {"rsid": "rs5219", "beta": 0.3, "se": 0.05, "p": 1e-10, "ea": "T", "nea": "C"},
        # This SNP is palindromic (A/T) with mismatched outcome alleles — should be excluded
        {"rsid": "rs9999", "beta": 0.2, "se": 0.04, "p": 1e-8, "ea": "A", "nea": "T"},
    ]
    outcomes = [
        {"rsid": "rs5219", "beta": 0.15, "se": 0.03, "p": 1e-6, "ea": "T", "nea": "C"},
        # Outcome alleles G/C don't match instrument A/T, and instrument is palindromic
        {"rsid": "rs9999", "beta": 0.08, "se": 0.02, "p": 0.001, "ea": "G", "nea": "C"},
    ]
    responses.add(responses.POST, f"{OPENGWAS_API_URL}/tophits", json=instruments, status=200)
    responses.add(responses.POST, f"{OPENGWAS_API_URL}/associations", json=outcomes, status=200)

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    row = result.iloc[0]
    # Palindromic SNP excluded, so only 1 instrument used — falls back to IVW with 1
    assert row["N_INSTRUMENTS"] == 1
    assert row["MR_STATUS"] == "ok"
