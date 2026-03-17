"""Map GWAS variants to genes using catalog annotations."""

import logging

import pandas as pd

from gwas_explorer.http_utils import create_session

logger = logging.getLogger(__name__)

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"

_session = create_session()


def _lookup_ensembl_ids(gene_symbols: list[str]) -> dict[str, str]:
    """Batch lookup Ensembl gene IDs for a list of gene symbols."""
    result: dict[str, str] = {}
    batch_size = 1000
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i : i + batch_size]
        try:
            response = _session.post(
                ENSEMBL_LOOKUP_URL,
                json={"symbols": batch},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=60,
            )
            response.raise_for_status()
            data = response.json()
            for symbol in batch:
                if symbol in data and "id" in data[symbol]:
                    result[symbol] = data[symbol]["id"]
                else:
                    result[symbol] = ""
        except Exception:
            logger.warning("Ensembl lookup failed for batch starting at index %d", i)
            for symbol in batch:
                result[symbol] = ""
    return result


def map_variants_to_genes(t2d_hits: pd.DataFrame) -> pd.DataFrame:
    """Map filtered T2D hits to unique genes with Ensembl IDs.

    Prioritizes MAPPED_GENE over REPORTED_GENES. Splits multi-gene entries.
    Deduplicates to one row per gene, keeping the strongest association.
    """
    if t2d_hits.empty:
        return pd.DataFrame(
            columns=["GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP", "P_VALUE",
                      "OR_BETA", "N_SUPPORTING_SNPS"]
        )

    df = t2d_hits.copy()

    # Choose gene source: MAPPED_GENE preferred, fall back to REPORTED_GENES
    mapped = df["MAPPED_GENE"].fillna("").astype(str).str.strip()
    reported = df["REPORTED_GENES"].fillna("").astype(str).str.strip()
    df["_GENE_SOURCE"] = mapped.where(mapped != "", reported)
    df = df[df["_GENE_SOURCE"] != ""]

    if df.empty:
        return pd.DataFrame(
            columns=["GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP", "P_VALUE",
                      "OR_BETA", "N_SUPPORTING_SNPS"]
        )

    # Normalize separators and explode multi-gene entries
    df["_GENE_SOURCE"] = df["_GENE_SOURCE"].str.replace(";", " - ", regex=False)
    df["_GENE_SOURCE"] = df["_GENE_SOURCE"].str.replace(",", " - ", regex=False)
    df["GENE_SYMBOL"] = df["_GENE_SOURCE"].str.split(" - ")
    df = df.explode("GENE_SYMBOL")
    df["GENE_SYMBOL"] = df["GENE_SYMBOL"].str.strip()
    df = df[df["GENE_SYMBOL"] != ""]

    # Count supporting SNPs per gene
    snp_counts = df.groupby("GENE_SYMBOL")["SNP_ID"].nunique().reset_index()
    snp_counts.columns = ["GENE_SYMBOL", "N_SUPPORTING_SNPS"]

    # Deduplicate: keep strongest association per gene
    deduped = df.sort_values("P_VALUE").drop_duplicates(subset=["GENE_SYMBOL"], keep="first")
    deduped = deduped.rename(columns={"SNP_ID": "LEAD_SNP"})
    deduped = deduped.merge(snp_counts, on="GENE_SYMBOL")

    # Resolve Ensembl IDs
    gene_symbols = deduped["GENE_SYMBOL"].tolist()
    ensembl_map = _lookup_ensembl_ids(gene_symbols)
    deduped["ENSEMBL_ID"] = deduped["GENE_SYMBOL"].map(ensembl_map)

    output_cols = [
        "GENE_SYMBOL",
        "ENSEMBL_ID",
        "LEAD_SNP",
        "P_VALUE",
        "OR_BETA",
        "N_SUPPORTING_SNPS",
    ]
    return deduped[output_cols].reset_index(drop=True)
