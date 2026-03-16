"""Tests for target prioritization and scoring module."""

import numpy as np
import pandas as pd
import pytest

from gwas_explorer.prioritize import prioritize_targets


@pytest.fixture
def full_pipeline_data() -> pd.DataFrame:
    """Combined data simulating full pipeline output for prioritization."""
    return pd.DataFrame(
        {
            "GENE_SYMBOL": ["KCNJ11", "TCF7L2", "PPARG", "ADRB2", "HLA-A"],
            "ENSEMBL_ID": [
                "ENSG00000187486",
                "ENSG00000148737",
                "ENSG00000132170",
                "ENSG00000169252",
                "ENSG00000206503",
            ],
            "LEAD_SNP": ["rs5219", "rs7903146", "rs1801282", "rs222222", "rs111111"],
            "P_VALUE": [5e-15, 1e-30, 2e-10, 1e-9, 1e-20],
            "OR_BETA": [1.2, 1.4, 1.15, 1.1, 1.5],
            "N_SUPPORTING_SNPS": [1, 2, 1, 1, 1],
            "IS_DRUGGABLE": [True, False, True, True, False],
            "TRACTABILITY_SM": [True, False, True, True, False],
            "TRACTABILITY_AB": [False, False, False, False, False],
            "HAS_EXISTING_DRUG": [True, False, True, True, False],
            "DRUGS": [
                [
                    {
                        "name": "GLIBENCLAMIDE",
                        "indication": "type 2 diabetes mellitus",
                        "phase": 4,
                        "mechanism": "",
                        "status": "Approved",
                    }
                ],
                [],
                [
                    {
                        "name": "PIOGLITAZONE",
                        "indication": "type 2 diabetes mellitus",
                        "phase": 4,
                        "mechanism": "",
                        "status": "Approved",
                    }
                ],
                [
                    {
                        "name": "SALBUTAMOL",
                        "indication": "asthma",
                        "phase": 4,
                        "mechanism": "",
                        "status": "Approved",
                    }
                ],
                [],
            ],
            "MR_ESTIMATE": [0.5, np.nan, 0.3, 0.1, np.nan],
            "MR_SE": [0.1, np.nan, 0.08, 0.2, np.nan],
            "MR_PVALUE": [1e-6, np.nan, 0.001, 0.6, np.nan],
            "N_INSTRUMENTS": [1, 0, 1, 1, 0],
            "METHOD": ["wald_ratio", "none", "wald_ratio", "wald_ratio", "none"],
            "MR_STATUS": ["ok", "insufficient_data", "ok", "ok", "error"],
        }
    )


def test_validated_target(full_pipeline_data: pd.DataFrame) -> None:
    """KCNJ11 has T2D drug -> should be 'validated'."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "KCNJ11"]
    assert row.iloc[0]["TARGET_CATEGORY"] == "validated"


def test_repurposing_candidate(full_pipeline_data: pd.DataFrame) -> None:
    """ADRB2 has drug but not for T2D -> should be 'repurposing_candidate'."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "ADRB2"]
    assert row.iloc[0]["TARGET_CATEGORY"] == "repurposing_candidate"


def test_novel_challenging(full_pipeline_data: pd.DataFrame) -> None:
    """TCF7L2 has no drugs and is not druggable -> should be 'novel_challenging'."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert row.iloc[0]["TARGET_CATEGORY"] == "novel_challenging"


def test_mr_causality_flag(full_pipeline_data: pd.DataFrame) -> None:
    """KCNJ11 MR p < 0.05 -> MR_SUPPORTS_CAUSALITY should be True."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "KCNJ11"]
    assert row.iloc[0]["MR_SUPPORTS_CAUSALITY"] is True


def test_mr_no_causality_flag(full_pipeline_data: pd.DataFrame) -> None:
    """ADRB2 MR p > 0.05 -> MR_SUPPORTS_CAUSALITY should be False."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "ADRB2"]
    assert row.iloc[0]["MR_SUPPORTS_CAUSALITY"] is False


def test_output_sorted_by_category_and_pvalue(full_pipeline_data: pd.DataFrame) -> None:
    """Output should be sorted: validated first, then by p-value within category."""
    result = prioritize_targets(full_pipeline_data)
    categories = result["TARGET_CATEGORY"].tolist()
    category_order = ["validated", "repurposing_candidate", "novel_druggable", "novel_challenging"]
    cat_indices = [category_order.index(c) for c in categories]
    assert cat_indices == sorted(cat_indices)
