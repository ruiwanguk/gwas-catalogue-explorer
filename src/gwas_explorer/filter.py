"""Filter GWAS Catalog associations for T2D genome-wide significant hits."""

import logging

import pandas as pd

from gwas_explorer.config import PVALUE_THRESHOLD, T2D_EFO_TERM, T2D_TRAIT_KEYWORDS

logger = logging.getLogger(__name__)


def filter_t2d_associations(associations: pd.DataFrame) -> pd.DataFrame:
    """Filter associations for T2D trait and genome-wide significance.

    Args:
        associations: Raw GWAS Catalog associations DataFrame.

    Returns:
        Filtered DataFrame with standardized column names.
    """
    df = associations.copy()

    # Filter by T2D trait: match EFO URI or keyword in disease/trait
    trait_col = df["DISEASE/TRAIT"].fillna("").str.lower()
    uri_col = df["MAPPED_TRAIT_URI"].fillna("")

    is_t2d_keyword = trait_col.apply(lambda x: any(kw in x for kw in T2D_TRAIT_KEYWORDS))
    is_t2d_efo = uri_col.str.contains(T2D_EFO_TERM, na=False)
    df = df[is_t2d_keyword | is_t2d_efo]

    # Filter by p-value (coerce non-numeric to NaN, which gets excluded by <)
    df["P-VALUE"] = pd.to_numeric(df["P-VALUE"], errors="coerce")
    df = df[df["P-VALUE"] < PVALUE_THRESHOLD]

    # Rename columns
    df = df.rename(
        columns={
            "SNPS": "SNP_ID",
            "CHR_ID": "CHR",
            "CHR_POS": "POSITION",
            "P-VALUE": "P_VALUE",
            "OR or BETA": "OR_BETA",
            "MAPPED_GENE": "MAPPED_GENE",
            "REPORTED GENE(S)": "REPORTED_GENES",
            "STUDY ACCESSION": "STUDY_ACCESSION",
        }
    )

    # Convert types
    df["POSITION"] = pd.to_numeric(df["POSITION"], errors="coerce")
    df["OR_BETA"] = pd.to_numeric(df["OR_BETA"], errors="coerce")

    # Deduplicate: keep strongest p-value per SNP
    df = df.sort_values("P_VALUE").drop_duplicates(subset=["SNP_ID"], keep="first")

    # Select output columns
    output_cols = [
        "SNP_ID",
        "CHR",
        "POSITION",
        "P_VALUE",
        "OR_BETA",
        "MAPPED_GENE",
        "REPORTED_GENES",
        "STUDY_ACCESSION",
    ]
    return df[output_cols].reset_index(drop=True)
