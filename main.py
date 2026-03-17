"""Run the full T2D drug target identification pipeline."""

import argparse
import logging

import pandas as pd

from gwas_explorer.config import CACHE_STALENESS_DAYS, PROCESSED_DIR, RESULTS_DIR
from gwas_explorer.download import download_gwas_catalog
from gwas_explorer.druggability import query_druggability
from gwas_explorer.filter import filter_t2d_associations
from gwas_explorer.gene_mapping import map_variants_to_genes
from gwas_explorer.mr_analysis import run_mr_analysis
from gwas_explorer.prioritize import prioritize_targets

logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description="T2D drug target identification pipeline")
    parser.add_argument("--max-age-days", type=int, default=CACHE_STALENESS_DAYS,
                        help="Skip download if cached file is newer than this (days)")
    parser.add_argument("--max-workers", type=int, default=4,
                        help="Number of concurrent API requests")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Download
    logger.info("Step 1: Downloading GWAS Catalog...")
    catalog_path = download_gwas_catalog(max_age_days=args.max_age_days)

    # Step 2: Filter
    logger.info("Step 2: Filtering T2D associations...")
    raw = pd.read_csv(catalog_path, sep="\t", low_memory=False)
    t2d_hits = filter_t2d_associations(raw)
    t2d_hits.to_parquet(PROCESSED_DIR / "t2d_significant_hits.parquet", index=False)
    logger.info("  Found %d significant hits", len(t2d_hits))

    # Step 3: Gene mapping
    logger.info("Step 3: Mapping variants to genes...")
    candidate_genes = map_variants_to_genes(t2d_hits)
    candidate_genes.to_parquet(PROCESSED_DIR / "t2d_candidate_genes.parquet", index=False)
    logger.info("  Identified %d candidate genes", len(candidate_genes))

    # Step 4: Druggability
    logger.info("Step 4: Querying druggability...")
    druggability = query_druggability(candidate_genes, max_workers=args.max_workers)
    druggability.to_parquet(PROCESSED_DIR / "t2d_druggability.parquet", index=False)
    logger.info("  Druggable: %d / %d", druggability["IS_DRUGGABLE"].sum(), len(druggability))

    # Step 5: Mendelian Randomization
    logger.info("Step 5: Running Mendelian Randomization...")
    mr_results = run_mr_analysis(candidate_genes, max_workers=args.max_workers)
    mr_results.to_parquet(PROCESSED_DIR / "t2d_mr_results.parquet", index=False)
    ok_count = (mr_results["MR_STATUS"] == "ok").sum()
    logger.info("  MR completed for %d / %d genes", ok_count, len(mr_results))

    # Step 6: Prioritize
    logger.info("Step 6: Prioritizing targets...")
    combined = druggability.merge(mr_results, on="GENE_SYMBOL", how="left")
    prioritized = prioritize_targets(combined)
    prioritized.to_parquet(RESULTS_DIR / "t2d_target_prioritization.parquet", index=False)
    csv_cols = [c for c in prioritized.columns if c != "DRUGS"]
    prioritized[csv_cols].to_csv(RESULTS_DIR / "t2d_target_prioritization.csv", index=False)

    # Summary
    logger.info("=== Results ===")
    logger.info("\n%s", prioritized["TARGET_CATEGORY"].value_counts().to_string())
    logger.info("Total targets: %d", len(prioritized))
    logger.info("Results saved to %s", RESULTS_DIR)


if __name__ == "__main__":
    main()
