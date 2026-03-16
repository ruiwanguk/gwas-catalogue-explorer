"""Simplified two-sample Mendelian Randomization using OpenGWAS API."""

import math

import numpy as np
import pandas as pd
import requests

from gwas_explorer.config import (
    MR_INSTRUMENT_PVALUE,
    OPENGWAS_API_URL,
    T2D_OUTCOME_ID,
)


def _get_instruments(snp_id: str, exposure_id: str | None = None) -> list[dict] | None:
    """Retrieve SNP-exposure associations from OpenGWAS.

    Returns a list of instruments, or None if the API request failed.
    """
    try:
        params: dict = {"rsid": snp_id, "proxies": 0}
        if exposure_id:
            params["id"] = exposure_id
        response = requests.get(
            f"{OPENGWAS_API_URL}/associations",
            params=params,
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        return [d for d in data if float(d.get("p", 1)) < MR_INSTRUMENT_PVALUE]
    except requests.RequestException:
        return None
    except ValueError:
        return []


def _get_outcome_associations(
    snp_ids: list[str], outcome_id: str = T2D_OUTCOME_ID
) -> dict[str, dict]:
    """Retrieve SNP-outcome associations for T2D from OpenGWAS."""
    try:
        response = requests.post(
            f"{OPENGWAS_API_URL}/associations",
            json={"rsid": snp_ids, "id": [outcome_id], "proxies": 0},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        return {d["rsid"]: d for d in data}
    except requests.RequestException, ValueError:
        return {}


def _wald_ratio(beta_exp: float, se_exp: float, beta_out: float, se_out: float) -> dict:
    """Compute Wald ratio MR estimate from a single SNP."""
    mr_estimate = beta_out / beta_exp
    mr_se = abs(se_out / beta_exp) * math.sqrt(1 + (se_exp / beta_exp) ** 2)
    z = abs(mr_estimate / mr_se) if mr_se > 0 else 0
    from scipy.stats import norm

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
        out = outcomes[rsid]
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

    from scipy.stats import norm

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


def run_mr_analysis(candidate_genes: pd.DataFrame) -> pd.DataFrame:
    """Run simplified MR for each candidate gene.

    Args:
        candidate_genes: DataFrame with GENE_SYMBOL and LEAD_SNP columns.

    Returns:
        DataFrame with MR results per gene.
    """
    results = []
    for _, row in candidate_genes.iterrows():
        gene = row["GENE_SYMBOL"]
        lead_snp = row["LEAD_SNP"]

        instruments = _get_instruments(lead_snp)

        # None means the API call itself failed
        if instruments is None:
            results.append(
                {
                    "GENE_SYMBOL": gene,
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "error",
                }
            )
            continue

        if not instruments:
            results.append(
                {
                    "GENE_SYMBOL": gene,
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "insufficient_data",
                }
            )
            continue

        snp_ids = [inst["rsid"] for inst in instruments]
        outcomes = _get_outcome_associations(snp_ids)

        if not outcomes:
            results.append(
                {
                    "GENE_SYMBOL": gene,
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "insufficient_data",
                }
            )
            continue

        if len(instruments) == 1:
            inst = instruments[0]
            rsid = inst["rsid"]
            if rsid in outcomes:
                out = outcomes[rsid]
                mr = _wald_ratio(
                    float(inst["beta"]),
                    float(inst["se"]),
                    float(out["beta"]),
                    float(out["se"]),
                )
            else:
                mr = {
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "insufficient_data",
                }
        else:
            mr = _ivw(instruments, outcomes)

        mr["GENE_SYMBOL"] = gene
        results.append(mr)

    return pd.DataFrame(results)[
        ["GENE_SYMBOL", "MR_ESTIMATE", "MR_SE", "MR_PVALUE", "N_INSTRUMENTS", "METHOD", "MR_STATUS"]
    ].reset_index(drop=True)
