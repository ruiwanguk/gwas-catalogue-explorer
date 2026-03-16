"""Tests for Open Targets druggability lookup module."""

import pandas as pd
import responses

from gwas_explorer.config import OPEN_TARGETS_GRAPHQL_URL
from gwas_explorer.druggability import query_druggability


@responses.activate
def test_identifies_druggable_target(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """KCNJ11 should be flagged as druggable with known drugs."""
    for _, row in sample_candidate_genes.iterrows():
        if row["GENE_SYMBOL"] == "KCNJ11":
            body = mock_open_targets_response
        else:
            body = mock_open_targets_novel_response
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, json=body, status=200)

    result = query_druggability(sample_candidate_genes)
    kcnj11 = result[result["GENE_SYMBOL"] == "KCNJ11"]
    assert kcnj11.iloc[0]["IS_DRUGGABLE"] is True
    assert kcnj11.iloc[0]["HAS_EXISTING_DRUG"] is True


@responses.activate
def test_identifies_novel_target(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """TCF7L2 should show no existing drugs."""
    for _, row in sample_candidate_genes.iterrows():
        if row["GENE_SYMBOL"] == "KCNJ11":
            body = mock_open_targets_response
        else:
            body = mock_open_targets_novel_response
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, json=body, status=200)

    result = query_druggability(sample_candidate_genes)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert tcf.iloc[0]["HAS_EXISTING_DRUG"] is False


@responses.activate
def test_extracts_drug_details(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """Should extract drug names and indications."""
    for _, row in sample_candidate_genes.iterrows():
        if row["GENE_SYMBOL"] == "KCNJ11":
            body = mock_open_targets_response
        else:
            body = mock_open_targets_novel_response
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, json=body, status=200)

    result = query_druggability(sample_candidate_genes)
    kcnj11 = result[result["GENE_SYMBOL"] == "KCNJ11"]
    drugs = kcnj11.iloc[0]["DRUGS"]
    assert any("GLIBENCLAMIDE" in d["name"] for d in drugs)


@responses.activate
def test_handles_api_error_gracefully(sample_candidate_genes: pd.DataFrame) -> None:
    """Should not crash on API errors; mark gene as unknown druggability."""
    for _ in range(len(sample_candidate_genes)):
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, status=500)

    result = query_druggability(sample_candidate_genes)
    assert len(result) == len(sample_candidate_genes)
    assert (result["IS_DRUGGABLE"] == False).all()  # noqa: E712
