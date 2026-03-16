"""Tests for variant-to-gene mapping module."""

import pandas as pd
import pytest
import responses

from gwas_explorer.gene_mapping import map_variants_to_genes


@pytest.fixture(autouse=True)
def mock_ensembl_api():
    """Mock Ensembl REST API for all tests to avoid real HTTP calls."""
    with responses.RequestsMock() as rsps:
        rsps.add(
            responses.POST,
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            json={
                "TCF7L2": {"id": "ENSG00000148737"},
                "PPARG": {"id": "ENSG00000132170"},
                "KCNJ11": {"id": "ENSG00000187486"},
                "ABCC8": {"id": "ENSG00000006071"},
                "HLA-A": {"id": "ENSG00000206503"},
                "ADRB2": {"id": "ENSG00000169252"},
            },
            status=200,
        )
        yield rsps


def test_prioritizes_mapped_gene(sample_t2d_hits: pd.DataFrame) -> None:
    """Should use MAPPED_GENE when present, not REPORTED_GENES."""
    result = map_variants_to_genes(sample_t2d_hits)
    # MAPPED_GENE for rs5219 is "KCNJ11", REPORTED_GENES is "KCNJ11 - ABCC8"
    # Since MAPPED_GENE is present, ABCC8 should NOT appear
    assert "KCNJ11" in result["GENE_SYMBOL"].values
    assert "ABCC8" not in result["GENE_SYMBOL"].values


def test_falls_back_to_reported_genes(sample_t2d_hits: pd.DataFrame) -> None:
    """Should use REPORTED_GENES when MAPPED_GENE is missing."""
    hits = sample_t2d_hits.copy()
    # Clear MAPPED_GENE for the KCNJ11 row, forcing fallback to REPORTED_GENES
    hits.loc[hits["SNP_ID"] == "rs5219", "MAPPED_GENE"] = ""
    result = map_variants_to_genes(hits)
    # Now ABCC8 should appear from the REPORTED_GENES fallback
    assert "KCNJ11" in result["GENE_SYMBOL"].values
    assert "ABCC8" in result["GENE_SYMBOL"].values


def test_deduplicates_genes(sample_t2d_hits: pd.DataFrame) -> None:
    """Should have one row per unique gene."""
    result = map_variants_to_genes(sample_t2d_hits)
    assert result["GENE_SYMBOL"].is_unique


def test_keeps_strongest_association(sample_t2d_hits: pd.DataFrame) -> None:
    """TCF7L2 appears multiple times; should keep lowest p-value."""
    result = map_variants_to_genes(sample_t2d_hits)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert len(tcf) == 1
    assert tcf.iloc[0]["P_VALUE"] == 1e-30


def test_counts_supporting_snps(sample_t2d_hits: pd.DataFrame) -> None:
    """TCF7L2 has 2 SNPs (rs7903146 and rs12255372)."""
    result = map_variants_to_genes(sample_t2d_hits)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert tcf.iloc[0]["N_SUPPORTING_SNPS"] == 2


def test_output_schema(sample_t2d_hits: pd.DataFrame) -> None:
    """Output should have the expected columns."""
    result = map_variants_to_genes(sample_t2d_hits)
    expected = {"GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP", "P_VALUE", "OR_BETA", "N_SUPPORTING_SNPS"}
    assert set(result.columns) == expected


def test_ensembl_id_lookup(sample_t2d_hits: pd.DataFrame) -> None:
    """Should resolve Ensembl IDs via Ensembl REST API."""
    result = map_variants_to_genes(sample_t2d_hits)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert tcf.iloc[0]["ENSEMBL_ID"] == "ENSG00000148737"
