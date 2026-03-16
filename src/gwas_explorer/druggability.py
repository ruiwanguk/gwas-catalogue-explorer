"""Query Open Targets Platform for target druggability information."""

import pandas as pd
import requests

from gwas_explorer.config import OPEN_TARGETS_GRAPHQL_URL

DRUGGABILITY_QUERY = """
query TargetDruggability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability {
      modality
      id
      value
    }
    knownDrugs {
      uniqueDrugs
      rows {
        drug {
          name
          mechanismOfAction
        }
        phase
        status
        disease {
          name
        }
      }
    }
  }
}
"""


def _query_single_target(ensembl_id: str) -> dict:
    """Query Open Targets for a single gene's druggability data."""
    try:
        response = requests.post(
            OPEN_TARGETS_GRAPHQL_URL,
            json={"query": DRUGGABILITY_QUERY, "variables": {"ensemblId": ensembl_id}},
            timeout=30,
        )
        response.raise_for_status()
        return response.json()
    except requests.RequestException, ValueError:
        return {}


def _parse_druggability(response: dict) -> dict:
    """Parse Open Targets response into druggability summary."""
    target = response.get("data", {}).get("target")
    if not target:
        return {
            "IS_DRUGGABLE": False,
            "TRACTABILITY_SM": False,
            "TRACTABILITY_AB": False,
            "HAS_EXISTING_DRUG": False,
            "DRUGS": [],
        }

    tractability = target.get("tractability") or []
    sm_tractable = any(t["value"] for t in tractability if t["modality"] == "SM")
    ab_tractable = any(t["value"] for t in tractability if t["modality"] == "AB")

    known_drugs = target.get("knownDrugs")
    drugs = []
    if known_drugs and known_drugs.get("rows"):
        for row in known_drugs["rows"]:
            drugs.append(
                {
                    "name": row["drug"]["name"],
                    "mechanism": row["drug"].get("mechanismOfAction", ""),
                    "phase": row.get("phase", 0),
                    "status": row.get("status", ""),
                    "indication": row.get("disease", {}).get("name", ""),
                }
            )

    return {
        "IS_DRUGGABLE": sm_tractable or ab_tractable or len(drugs) > 0,
        "TRACTABILITY_SM": sm_tractable,
        "TRACTABILITY_AB": ab_tractable,
        "HAS_EXISTING_DRUG": len(drugs) > 0,
        "DRUGS": drugs,
    }


def query_druggability(candidate_genes: pd.DataFrame) -> pd.DataFrame:
    """Query Open Targets for druggability of all candidate genes.

    Args:
        candidate_genes: DataFrame with GENE_SYMBOL and ENSEMBL_ID columns.

    Returns:
        Input DataFrame merged with druggability columns.
    """
    records = []
    for _, row in candidate_genes.iterrows():
        ensembl_id = row["ENSEMBL_ID"]
        response = _query_single_target(ensembl_id)
        drug_info = _parse_druggability(response)
        drug_info["GENE_SYMBOL"] = row["GENE_SYMBOL"]
        records.append(drug_info)

    drug_df = pd.DataFrame(records)
    result = candidate_genes.merge(drug_df, on="GENE_SYMBOL", how="left")
    # Store as Python bools (object dtype) so identity checks work.
    for col in ("IS_DRUGGABLE", "TRACTABILITY_SM", "TRACTABILITY_AB", "HAS_EXISTING_DRUG"):
        result[col] = pd.array([bool(x) for x in result[col]], dtype=object)
    return result
