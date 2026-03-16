"""Run the full T2D drug target identification pipeline."""

import pandas as pd

from gwas_explorer.config import PROCESSED_DIR, RESULTS_DIR
from gwas_explorer.download import download_gwas_catalog
from gwas_explorer.druggability import query_druggability
from gwas_explorer.filter import filter_t2d_associations
from gwas_explorer.gene_mapping import map_variants_to_genes
from gwas_explorer.mr_analysis import run_mr_analysis
from gwas_explorer.prioritize import prioritize_targets


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Download
    print("Step 1: Downloading GWAS Catalog...")
    catalog_path = download_gwas_catalog()

    # Step 2: Filter
    print("Step 2: Filtering T2D associations...")
    raw = pd.read_csv(catalog_path, sep="\t", low_memory=False)
    t2d_hits = filter_t2d_associations(raw)
    t2d_hits.to_parquet(PROCESSED_DIR / "t2d_significant_hits.parquet", index=False)
    print(f"  Found {len(t2d_hits)} significant hits")

    # Step 3: Gene mapping
    print("Step 3: Mapping variants to genes...")
    candidate_genes = map_variants_to_genes(t2d_hits)
    candidate_genes.to_parquet(PROCESSED_DIR / "t2d_candidate_genes.parquet", index=False)
    print(f"  Identified {len(candidate_genes)} candidate genes")

    # Step 4: Druggability
    print("Step 4: Querying druggability...")
    druggability = query_druggability(candidate_genes)
    druggability.to_parquet(PROCESSED_DIR / "t2d_druggability.parquet", index=False)
    print(f"  Druggable: {druggability['IS_DRUGGABLE'].sum()} / {len(druggability)}")

    # Step 5: Mendelian Randomization
    print("Step 5: Running Mendelian Randomization...")
    mr_results = run_mr_analysis(candidate_genes)
    mr_results.to_parquet(PROCESSED_DIR / "t2d_mr_results.parquet", index=False)
    ok_count = (mr_results["MR_STATUS"] == "ok").sum()
    print(f"  MR completed for {ok_count} / {len(mr_results)} genes")

    # Step 6: Prioritize
    print("Step 6: Prioritizing targets...")
    combined = druggability.merge(mr_results, on="GENE_SYMBOL", how="left")
    prioritized = prioritize_targets(combined)
    prioritized.to_parquet(RESULTS_DIR / "t2d_target_prioritization.parquet", index=False)
    csv_cols = [c for c in prioritized.columns if c != "DRUGS"]
    prioritized[csv_cols].to_csv(RESULTS_DIR / "t2d_target_prioritization.csv", index=False)

    # Summary
    print("\n=== Results ===")
    print(prioritized["TARGET_CATEGORY"].value_counts().to_string())
    print(f"\nTotal targets: {len(prioritized)}")
    print(f"Results saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
