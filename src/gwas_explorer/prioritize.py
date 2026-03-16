"""Prioritize and categorize T2D drug target candidates."""

import pandas as pd

from gwas_explorer.config import MR_PVALUE_THRESHOLD

T2D_INDICATION_KEYWORDS = ["diabetes", "t2d", "type 2", "type ii"]

CATEGORY_ORDER = ["validated", "repurposing_candidate", "novel_druggable", "novel_challenging"]


def _matches_t2d_indication(drugs: list[dict]) -> bool:
    """Check if any drug in the list is indicated for T2D."""
    for drug in drugs:
        indication = drug.get("indication", "").lower()
        if any(kw in indication for kw in T2D_INDICATION_KEYWORDS):
            return True
    return False


def _assign_category(row: pd.Series) -> str:
    """Assign target category based on druggability and drug data."""
    if row["HAS_EXISTING_DRUG"] and row["DRUG_INDICATION_MATCHES_T2D"]:
        return "validated"
    if row["HAS_EXISTING_DRUG"] and not row["DRUG_INDICATION_MATCHES_T2D"]:
        return "repurposing_candidate"
    if row["IS_DRUGGABLE"]:
        return "novel_druggable"
    return "novel_challenging"


def prioritize_targets(combined_data: pd.DataFrame) -> pd.DataFrame:
    """Score and categorize drug target candidates."""
    df = combined_data.copy()

    df["HAS_GWAS_SIGNAL"] = True
    df["DRUG_INDICATION_MATCHES_T2D"] = df["DRUGS"].apply(
        lambda drugs: _matches_t2d_indication(drugs) if isinstance(drugs, list) else False
    )
    df["MR_SUPPORTS_CAUSALITY"] = (
        df["MR_PVALUE"]
        .apply(lambda p: bool(p < MR_PVALUE_THRESHOLD) if pd.notna(p) else False)
        .astype(object)
    )

    df["TARGET_CATEGORY"] = df.apply(_assign_category, axis=1)

    df["_cat_rank"] = df["TARGET_CATEGORY"].map({c: i for i, c in enumerate(CATEGORY_ORDER)})
    df = df.sort_values(["_cat_rank", "P_VALUE"]).drop(columns=["_cat_rank"])

    return df.reset_index(drop=True)
