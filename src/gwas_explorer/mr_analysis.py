"""Simplified two-sample Mendelian Randomization using OpenGWAS API."""

import logging
import math

import numpy as np
import pandas as pd
from scipy.stats import norm

from gwas_explorer.config import (
    EQTLGEN_STUDY_PREFIX,
    MR_INSTRUMENT_PVALUE,
    OPENGWAS_API_URL,
    OPENGWAS_TOKEN,
    T2D_OUTCOME_ID,
)
from gwas_explorer.http_utils import RateLimitedExecutor, create_session

logger = logging.getLogger(__name__)

_session = create_session()

if OPENGWAS_TOKEN:
    _session.headers["Authorization"] = f"Bearer {OPENGWAS_TOKEN}"
else:
    logger.warning(
        "OPENGWAS_TOKEN not set. OpenGWAS API requires authentication — "
        "MR instrument lookups will likely fail. "
        "See https://api.opengwas.io for token registration."
    )

_INSUFFICIENT = {
    "MR_ESTIMATE": np.nan,
    "MR_SE": np.nan,
    "MR_PVALUE": np.nan,
    "N_INSTRUMENTS": 0,
    "METHOD": "none",
    "MR_STATUS": "insufficient_data",
}


def _get_eqtl_instruments(ensembl_id: str) -> list[dict] | None:
    """Retrieve genome-wide significant cis-eQTL instruments from eQTLGen via OpenGWAS.

    Uses the /tophits endpoint to get the strongest eQTL SNPs for the gene.
    Returns a list of instruments, or None if the API request failed.
    """
    eqtl_study_id = f"{EQTLGEN_STUDY_PREFIX}{ensembl_id}"
    try:
        response = _session.post(
            f"{OPENGWAS_API_URL}/tophits",
            json={"id": [eqtl_study_id]},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        if not isinstance(data, list):
            return None
        return [d for d in data if float(d.get("p", 1)) < MR_INSTRUMENT_PVALUE]
    except Exception:
        logger.warning(
            "Failed to get eQTL instruments for %s (%s)", ensembl_id, eqtl_study_id, exc_info=True
        )
        return None


def _get_outcome_associations(
    snp_ids: list[str], outcome_id: str = T2D_OUTCOME_ID
) -> dict[str, dict]:
    """Retrieve SNP-outcome associations for T2D from OpenGWAS."""
    try:
        response = _session.post(
            f"{OPENGWAS_API_URL}/associations",
            json={"variant": snp_ids, "id": [outcome_id], "proxies": 0},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        return {d["rsid"]: d for d in data}
    except Exception:
        logger.warning("Failed to get outcome associations for %s", snp_ids, exc_info=True)
        return {}


def _harmonize_alleles(inst: dict, out: dict) -> dict | None:
    """Harmonize effect alleles between exposure and outcome.

    Returns the outcome dict with beta sign-flipped if alleles are swapped,
    or None if alleles cannot be aligned (e.g., palindromic ambiguity).
    """
    exp_ea = inst.get("ea", "").upper()
    exp_nea = inst.get("nea", "").upper()
    out_ea = out.get("ea", "").upper()
    out_nea = out.get("nea", "").upper()

    if not all([exp_ea, exp_nea, out_ea, out_nea]):
        return out  # Missing allele info — pass through unchanged

    # Already aligned
    if exp_ea == out_ea and exp_nea == out_nea:
        return out

    # Alleles swapped — flip outcome beta
    if exp_ea == out_nea and exp_nea == out_ea:
        flipped = dict(out)
        flipped["beta"] = -float(out["beta"])
        return flipped

    # Palindromic SNPs (A/T or C/G) — cannot resolve without frequency info
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if exp_ea == complement.get(exp_nea, ""):
        logger.debug("Excluding palindromic SNP %s", inst.get("rsid"))
        return None

    return out  # Non-matching, non-palindromic — pass through


def _wald_ratio(beta_exp: float, se_exp: float, beta_out: float, se_out: float) -> dict:
    """Compute Wald ratio MR estimate from a single SNP."""
    if abs(beta_exp) < 1e-10:
        return {
            "MR_ESTIMATE": np.nan,
            "MR_SE": np.nan,
            "MR_PVALUE": np.nan,
            "N_INSTRUMENTS": 0,
            "METHOD": "wald_ratio",
            "MR_STATUS": "invalid_instrument",
        }

    mr_estimate = beta_out / beta_exp
    mr_se = abs(se_out / beta_exp) * math.sqrt(1 + (se_exp / beta_exp) ** 2)
    z = abs(mr_estimate / mr_se) if mr_se > 0 else 0
    mr_pvalue = 2 * norm.sf(z)
    return {
        "MR_ESTIMATE": mr_estimate,
        "MR_SE": mr_se,
        "MR_PVALUE": mr_pvalue,
        "N_INSTRUMENTS": 1,
        "METHOD": "wald_ratio",
        "MR_STATUS": "ok",
    }


def _ivw(instruments: list[dict], outcomes: dict[str, dict]) -> dict:
    """Inverse-variance weighted MR estimate from multiple SNPs."""
    estimates = []
    weights = []

    for inst in instruments:
        rsid = inst["rsid"]
        if rsid not in outcomes:
            continue
        out = _harmonize_alleles(inst, outcomes[rsid])
        if out is None:
            continue
        beta_exp = float(inst["beta"])
        beta_out = float(out["beta"])
        se_out = float(out["se"])

        if abs(beta_exp) < 1e-10:
            continue

        ratio = beta_out / beta_exp
        ratio_se = abs(se_out / beta_exp)
        if ratio_se > 0:
            w = 1 / (ratio_se**2)
            estimates.append(ratio)
            weights.append(w)

    if not estimates:
        return {
            "MR_ESTIMATE": np.nan,
            "MR_SE": np.nan,
            "MR_PVALUE": np.nan,
            "N_INSTRUMENTS": 0,
            "METHOD": "ivw",
            "MR_STATUS": "insufficient_data",
        }

    weights_arr = np.array(weights)
    estimates_arr = np.array(estimates)
    ivw_estimate = np.sum(weights_arr * estimates_arr) / np.sum(weights_arr)
    ivw_se = math.sqrt(1 / np.sum(weights_arr))

    z = abs(ivw_estimate / ivw_se) if ivw_se > 0 else 0
    mr_pvalue = 2 * norm.sf(z)

    return {
        "MR_ESTIMATE": float(ivw_estimate),
        "MR_SE": float(ivw_se),
        "MR_PVALUE": float(mr_pvalue),
        "N_INSTRUMENTS": len(estimates),
        "METHOD": "ivw",
        "MR_STATUS": "ok",
    }


def _run_mr_for_gene(gene: str, ensembl_id: str) -> dict:
    """Run two-sample MR for a single gene using eQTL exposure and T2D outcome."""
    instruments = _get_eqtl_instruments(ensembl_id)

    if instruments is None:
        return {"GENE_SYMBOL": gene, **_INSUFFICIENT, "MR_STATUS": "error"}

    if not instruments:
        return {"GENE_SYMBOL": gene, **_INSUFFICIENT, "MR_STATUS": "no_eqtl_instruments"}

    snp_ids = [inst["rsid"] for inst in instruments]
    outcomes = _get_outcome_associations(snp_ids)

    if not outcomes:
        return {"GENE_SYMBOL": gene, **_INSUFFICIENT}

    if len(instruments) == 1:
        inst = instruments[0]
        rsid = inst["rsid"]
        if rsid in outcomes:
            out = _harmonize_alleles(inst, outcomes[rsid])
            if out is None:
                mr = dict(_INSUFFICIENT)
            else:
                mr = _wald_ratio(
                    float(inst["beta"]),
                    float(inst["se"]),
                    float(out["beta"]),
                    float(out["se"]),
                )
        else:
            mr = dict(_INSUFFICIENT)
    else:
        mr = _ivw(instruments, outcomes)

    mr["GENE_SYMBOL"] = gene
    return mr


def run_mr_analysis(candidate_genes: pd.DataFrame, max_workers: int = 4) -> pd.DataFrame:
    """Run two-sample MR (eQTL exposure → T2D outcome) for each candidate gene.

    Args:
        candidate_genes: DataFrame with GENE_SYMBOL and ENSEMBL_ID columns.
        max_workers: Number of concurrent API requests (default 4).

    Returns:
        DataFrame with MR results per gene.
    """
    if candidate_genes.empty:
        return pd.DataFrame(
            columns=[
                "GENE_SYMBOL",
                "MR_ESTIMATE",
                "MR_SE",
                "MR_PVALUE",
                "N_INSTRUMENTS",
                "METHOD",
                "MR_STATUS",
            ]
        )

    genes = list(zip(candidate_genes["GENE_SYMBOL"], candidate_genes["ENSEMBL_ID"]))
    executor = RateLimitedExecutor(max_workers=max_workers, min_delay=0.1)
    results = executor.map(_run_mr_for_gene, genes)

    return pd.DataFrame(results)[
        ["GENE_SYMBOL", "MR_ESTIMATE", "MR_SE", "MR_PVALUE", "N_INSTRUMENTS", "METHOD", "MR_STATUS"]
    ].reset_index(drop=True)
