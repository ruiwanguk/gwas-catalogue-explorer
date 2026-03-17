# GWAS Catalogue Explorer

T2D drug target identification and validation pipeline. Uses GWAS Catalog data to identify genome-wide significant variants for Type 2 Diabetes, maps them to genes, assesses druggability, runs simplified Mendelian Randomization for causal evidence, and prioritizes targets into actionable categories.

## Architecture

The project has two layers:

- **`src/gwas_explorer/`** — Core logic modules with DataFrame-in/DataFrame-out interfaces
- **`notebooks/`** — Step-by-step orchestration and visualization, importing from `src/`

Intermediate results are cached as parquet files in `data/` to avoid redundant downloads and reprocessing.

### Pipeline Flow

```
GWAS Catalog (EBI FTP)
  → download.py    Download & extract associations ZIP
  → filter.py      Filter for T2D trait (EFO_0001360) + genome-wide significance (p < 5e-8)
  → gene_mapping.py    Split multi-gene entries, deduplicate, resolve Ensembl IDs
  → druggability.py    Query Open Targets GraphQL for tractability & known drugs
  → mr_analysis.py     Simplified two-sample MR (Wald ratio / IVW) via OpenGWAS API
  → prioritize.py      Categorize targets into four tiers
```

### Target Categories

| Category | Criteria |
|---|---|
| **Validated** | Existing drug indicated for T2D (e.g., KCNJ11/sulfonylureas) |
| **Repurposing candidate** | Existing drug, but for a different indication |
| **Novel druggable** | No existing drug, but tractable protein (SM or antibody) |
| **Novel challenging** | No existing drug, not easily tractable |

### External APIs

| Service | Purpose |
|---|---|
| [GWAS Catalog FTP](https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/) | Bulk association data download |
| [Ensembl REST API](https://rest.ensembl.org) | Gene symbol to Ensembl ID mapping |
| [Open Targets GraphQL](https://api.platform.opentargets.org/api/v4/graphql) | Tractability assessment and known drugs |
| [OpenGWAS API](https://api.opengwas.io/api) | SNP-exposure and SNP-outcome associations for MR |

## Project Structure

```
src/gwas_explorer/
  config.py           API endpoints, thresholds, file paths
  download.py         Bulk download + ZIP extraction with staleness-based caching
  filter.py           T2D trait + significance filtering
  gene_mapping.py     Variant-to-gene mapping with Ensembl ID resolution
  druggability.py     Open Targets druggability queries (concurrent)
  mr_analysis.py      Two-sample MR analysis (concurrent)
  prioritize.py       Target scoring and categorization

notebooks/
  01_extract_gwas_hits.ipynb
  02_map_variants_to_genes.ipynb
  03_druggability_lookup.ipynb
  04_mendelian_randomization.ipynb
  05_prioritize_targets.ipynb

tests/                pytest tests with mocked HTTP (no network access)
data/
  raw/                Downloaded bulk files
  processed/          Intermediate parquet files
  results/            Final prioritization output (parquet + CSV)
```

## Getting Started

Requires Python 3.14+ and [uv](https://docs.astral.sh/uv/).

```bash
# Install dependencies
uv sync

# Run the full pipeline
uv run python main.py

# Or run step-by-step via notebooks
uv run jupyter notebook
```

## Development

```bash
# Run tests
uv run pytest

# Lint and format
uv run ruff check src/ tests/
uv run ruff format src/ tests/

# Type check
uv run pyright src/
```

## Existing Features

- **Bulk data download** — Downloads GWAS Catalog associations ZIP from EBI FTP with configurable staleness-based caching (default 30 days)
- **T2D filtering** — Filters by EFO ontology term (`EFO_0001360`) and trait keywords, applies genome-wide significance threshold (p < 5e-8), deduplicates by lead SNP
- **Gene mapping** — Splits multi-gene annotations (`;`, `,`, ` - ` delimiters), prioritizes `MAPPED_GENE` over reported genes, batch-resolves Ensembl IDs
- **Druggability assessment** — Queries Open Targets for small molecule and antibody tractability, extracts known drugs with phase, mechanism, and indication; concurrent API requests
- **Mendelian Randomization** — Simplified two-sample MR using OpenGWAS; Wald ratio for single instruments, inverse-variance weighted (IVW) for multiple; concurrent API requests
- **Target prioritization** — Four-tier categorization combining druggability, drug indication matching, and MR causal evidence
- **Notebooks** — Five Jupyter notebooks covering each pipeline step with visualizations (Manhattan-style plots, bar charts, pie charts)
- **Test suite** — 30 pytest tests with full HTTP mocking via `responses` library; no network access during tests

## Missing Features / Known Limitations

- **LD clumping** — No linkage disequilibrium clumping; nearby variants in the same locus may map to different genes, inflating candidate counts
- **Colocalization** — No colocalization analysis (e.g., coloc, eCAVIAR) to distinguish shared vs. distinct causal variants between GWAS and eQTL signals
- **eQTL integration** — No expression QTL data to link variants to gene expression changes and strengthen causal gene assignment
- **MR diagnostics** — No MR-Egger, weighted median, or MR-PRESSO to assess pleiotropy and instrument validity; current MR is simplified
- **Pathway enrichment** — No gene set enrichment or pathway analysis on prioritized targets
- **Tissue-specific expression** — No filtering by pancreas/adipose/liver expression to contextualize target relevance
- **PheWAS** — No phenome-wide association scan to assess pleiotropy risk of candidate targets
- **Progress reporting** — No progress bars for long-running API queries (druggability and MR steps)
- **Rate limiting** — No adaptive rate limiting or retry-with-backoff for API failures
- **Multi-trait support** — Hardcoded for T2D; not configurable for other diseases without code changes
