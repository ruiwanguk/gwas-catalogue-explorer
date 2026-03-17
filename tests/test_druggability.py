"""Tests for Open Targets druggability lookup module."""

import json

import pandas as pd
import responses

from gwas_explorer.config import OPEN_TARGETS_GRAPHQL_URL
from gwas_explorer.druggability import query_druggability


def _make_ot_callback(druggable_response: dict, novel_response: dict):
    """Return a responses callback that routes by Ensembl ID in the request body."""
    kcnj11_id = "ENSG00000187486"

    def callback(request):
        body = json.loads(request.body)
        ensembl_id = body.get("variables", {}).get("ensemblId", "")
        resp = druggable_response if ensembl_id == kcnj11_id else novel_response
        return (200, {}, json.dumps(resp))

    return callback


@responses.activate
def test_identifies_druggable_target(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """KCNJ11 should be flagged as druggable with known drugs."""
    responses.add_callback(
        responses.POST,
        OPEN_TARGETS_GRAPHQL_URL,
        callback=_make_ot_callback(mock_open_targets_response, mock_open_targets_novel_response),
    )

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
    responses.add_callback(
        responses.POST,
        OPEN_TARGETS_GRAPHQL_URL,
        callback=_make_ot_callback(mock_open_targets_response, mock_open_targets_novel_response),
    )

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
    responses.add_callback(
        responses.POST,
        OPEN_TARGETS_GRAPHQL_URL,
        callback=_make_ot_callback(mock_open_targets_response, mock_open_targets_novel_response),
    )

    result = query_druggability(sample_candidate_genes)
    kcnj11 = result[result["GENE_SYMBOL"] == "KCNJ11"]
    drugs = kcnj11.iloc[0]["DRUGS"]
    assert any("GLIBENCLAMIDE" in d["name"] for d in drugs)


@responses.activate
def test_handles_api_error_gracefully(sample_candidate_genes: pd.DataFrame) -> None:
    """Should not crash on API errors; mark gene as unknown druggability."""
    responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, status=500)

    result = query_druggability(sample_candidate_genes)
    assert len(result) == len(sample_candidate_genes)
    assert (result["IS_DRUGGABLE"] == False).all()  # noqa: E712


@responses.activate
def test_empty_candidates_returns_empty() -> None:
    """Should return empty DataFrame with correct columns for empty input."""
    empty = pd.DataFrame(columns=["GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP",
                                   "P_VALUE", "OR_BETA", "N_SUPPORTING_SNPS"])
    result = query_druggability(empty)
    assert len(result) == 0
    assert "IS_DRUGGABLE" in result.columns
    assert "DRUGS" in result.columns
