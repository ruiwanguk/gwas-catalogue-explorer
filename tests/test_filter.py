"""Tests for T2D association filtering module."""

import pandas as pd

from gwas_explorer.filter import filter_t2d_associations


def test_filters_by_trait(sample_gwas_associations: pd.DataFrame) -> None:
    """Should keep only T2D-related associations."""
    result = filter_t2d_associations(sample_gwas_associations)
    # "Height" row (rs999999) should be excluded
    assert "rs999999" not in result["SNP_ID"].values


def test_filters_by_pvalue(sample_gwas_associations: pd.DataFrame) -> None:
    """Should keep only genome-wide significant hits."""
    result = filter_t2d_associations(sample_gwas_associations)
    assert (result["P_VALUE"] < 5e-8).all()


def test_deduplicates_by_snp(sample_gwas_associations: pd.DataFrame) -> None:
    """Should keep one row per SNP (strongest p-value)."""
    result = filter_t2d_associations(sample_gwas_associations)
    assert result["SNP_ID"].is_unique


def test_output_schema(sample_gwas_associations: pd.DataFrame) -> None:
    """Output should have the expected columns."""
    result = filter_t2d_associations(sample_gwas_associations)
    expected_cols = {
        "SNP_ID",
        "CHR",
        "POSITION",
        "P_VALUE",
        "OR_BETA",
        "MAPPED_GENE",
        "REPORTED_GENES",
        "STUDY_ACCESSION",
    }
    assert set(result.columns) == expected_cols


def test_rs7903146_kept_with_best_pvalue(sample_gwas_associations: pd.DataFrame) -> None:
    """rs7903146 appears twice for T2D; should keep the one with p=1e-30."""
    result = filter_t2d_associations(sample_gwas_associations)
    row = result[result["SNP_ID"] == "rs7903146"]
    assert len(row) == 1
    assert row.iloc[0]["P_VALUE"] == 1e-30


def test_handles_non_numeric_pvalue(sample_gwas_associations: pd.DataFrame) -> None:
    """Should skip rows with non-numeric p-values instead of crashing."""
    df = sample_gwas_associations.copy()
    # Add a row with a non-numeric p-value
    bad_row = df.iloc[0].copy()
    bad_row["P-VALUE"] = "NS"
    bad_row["SNPS"] = "rs_bad"
    df = pd.concat([df, bad_row.to_frame().T], ignore_index=True)
    result = filter_t2d_associations(df)
    assert "rs_bad" not in result["SNP_ID"].values


def test_handles_nan_disease_trait(sample_gwas_associations: pd.DataFrame) -> None:
    """Should not crash when DISEASE/TRAIT is NaN."""
    df = sample_gwas_associations.copy()
    df.loc[0, "DISEASE/TRAIT"] = float("nan")
    result = filter_t2d_associations(df)
    assert len(result) > 0
