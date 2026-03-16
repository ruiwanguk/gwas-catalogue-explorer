"""Map GWAS variants to genes using catalog annotations."""

import pandas as pd
import requests

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"


def _lookup_ensembl_ids(gene_symbols: list[str]) -> dict[str, str]:
    """Batch lookup Ensembl gene IDs for a list of gene symbols."""
    result: dict[str, str] = {}
    batch_size = 1000
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i : i + batch_size]
        response = requests.post(
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
    return result


def map_variants_to_genes(t2d_hits: pd.DataFrame) -> pd.DataFrame:
    """Map filtered T2D hits to unique genes with Ensembl IDs.

    Prioritizes MAPPED_GENE over REPORTED_GENES. Splits multi-gene entries.
    Deduplicates to one row per gene, keeping the strongest association.
    """
    rows = []
    for _, hit in t2d_hits.iterrows():
        # Prioritize MAPPED_GENE; fall back to REPORTED_GENES only when missing
        gene_str = hit.get("MAPPED_GENE", "")
        if pd.isna(gene_str) or not str(gene_str).strip():
            gene_str = hit.get("REPORTED_GENES", "")
        if pd.isna(gene_str) or not str(gene_str).strip():
            continue

        gene_str = str(gene_str)
        genes = [
            g.strip()
            for g in gene_str.replace(";", " - ").replace(",", " - ").split(" - ")
            if g.strip()
        ]

        for gene in genes:
            if gene:
                rows.append(
                    {
                        "GENE_SYMBOL": gene,
                        "SNP_ID": hit["SNP_ID"],
                        "P_VALUE": hit["P_VALUE"],
                        "OR_BETA": hit["OR_BETA"],
                    }
                )

    expanded = pd.DataFrame(rows)

    snp_counts = expanded.groupby("GENE_SYMBOL")["SNP_ID"].nunique().reset_index()
    snp_counts.columns = ["GENE_SYMBOL", "N_SUPPORTING_SNPS"]

    deduped = expanded.sort_values("P_VALUE").drop_duplicates(subset=["GENE_SYMBOL"], keep="first")
    deduped = deduped.rename(columns={"SNP_ID": "LEAD_SNP"})
    deduped = deduped.merge(snp_counts, on="GENE_SYMBOL")

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
