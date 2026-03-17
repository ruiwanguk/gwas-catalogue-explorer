# T2D Drug Target Identification & Validation Pipeline

## Goal

Identify and validate drug targets for Type 2 Diabetes using GWAS Catalog data, druggability databases, and Mendelian Randomization. The pipeline produces a prioritized list of gene targets categorized as validated, repurposing candidates, novel druggable, or novel challenging.

## Architecture

Two layers:

- **`src/gwas_explorer/`** — Core logic modules, each with well-defined DataFrame inputs/outputs
- **`notebooks/`** — Orchestration and visualization, importing from `src/`

Intermediate results cached as parquet files in `data/` to avoid re-downloading or re-processing.

### Data Flow

```
GWAS Catalog bulk TSV
  → filter.py (T2D, p < 5e-8)
  → gene_mapping.py (variant-to-gene, Ensembl IDs)
  → druggability.py (Open Targets GraphQL)
  → mr_analysis.py (OpenGWAS / IEU API)
  → prioritize.py (scoring & categorization)
```

## Module Specifications

### `config.py`

Constants and configuration:

- EBI FTP URL for GWAS Catalog full associations file
- EFO ontology term for T2D: `EFO_0001360`
- P-value threshold: `5e-8`
- Open Targets Platform GraphQL endpoint
- OpenGWAS (IEU) API endpoint
- Cache directory paths (`data/raw/`, `data/processed/`, `data/results/`)

### `download.py`

Downloads the GWAS Catalog full associations TSV from the EBI FTP site.

- Caches in `data/raw/` with date stamp
- Skips download if a recent file exists (configurable staleness threshold)
- Returns path to local file

### `filter.py`

Loads raw associations and filters for T2D genome-wide significant hits.

**Filters:**
- Trait matching: `DISEASE/TRAIT` text match + `MAPPED_TRAIT_URI` containing `EFO_0001360`
- Significance: `P-VALUE` < 5e-8
- Deduplication by lead SNP per locus

**Output schema:** `SNP_ID`, `CHR`, `POSITION`, `P_VALUE`, `OR_BETA`, `MAPPED_GENE`, `STUDY_ACCESSION`

**Output file:** `data/processed/t2d_significant_hits.parquet`

### `gene_mapping.py`

Extracts gene assignments from GWAS Catalog annotations.

- Splits multi-gene entries (e.g., `"KCNJ11 - ABCC8"`) into individual rows
- Prioritizes `MAPPED_GENE` over `REPORTED GENE(S)`
- Deduplicates to unique gene list with strongest association per gene
- Maps gene symbols to Ensembl gene IDs via lookup

**Output schema:** `GENE_SYMBOL`, `ENSEMBL_ID`, `LEAD_SNP`, `P_VALUE`, `OR_BETA`, `N_SUPPORTING_SNPS`

**Output file:** `data/processed/t2d_candidate_genes.parquet`

### `druggability.py`

Queries Open Targets Platform GraphQL API for each candidate gene.

**Retrieves:**
- Tractability assessments (small molecule, antibody, other modalities)
- Known drugs: name, phase, mechanism of action, indication
- Target safety information

Batches queries to respect rate limits.

**Output file:** `data/processed/t2d_druggability.parquet`

### `mr_analysis.py`

Simplified two-sample Mendelian Randomization using OpenGWAS (IEU) API.

**Per candidate gene:**
- Selects lead SNP(s) near gene as instruments
- Retrieves SNP-exposure and SNP-outcome associations from OpenGWAS
- Computes Wald ratio (single SNP) or IVW estimate (multiple SNPs)
- Returns causal effect estimate, CI, and p-value

Genes with insufficient instruments flagged as `MR_STATUS: insufficient_data`.

**Output schema:** `GENE_SYMBOL`, `MR_ESTIMATE`, `MR_SE`, `MR_PVALUE`, `N_INSTRUMENTS`, `METHOD`, `MR_STATUS`

**Output file:** `data/processed/t2d_mr_results.parquet`

### `prioritize.py`

Combines all upstream results into scored ranking.

**Per-gene flags:**
- `HAS_GWAS_SIGNAL`: always true
- `IS_DRUGGABLE`: from tractability assessment
- `HAS_EXISTING_DRUG`: boolean + drug details
- `DRUG_INDICATION_MATCHES_T2D`: existing drug indicated for diabetes
- `MR_SUPPORTS_CAUSALITY`: MR p-value < 0.05

**Target categories:**
- `validated` — existing T2D drug targets (e.g., KCNJ11/sulfonylureas)
- `repurposing_candidate` — drug exists but for other indication
- `novel_druggable` — no drug, but tractable protein
- `novel_challenging` — no drug, not easily tractable

**Output files:** `data/results/t2d_target_prioritization.parquet`, `data/results/t2d_target_prioritization.csv`

## Project Structure

```
src/
  gwas_explorer/
    __init__.py
    config.py
    download.py
    filter.py
    gene_mapping.py
    druggability.py
    mr_analysis.py
    prioritize.py

tests/
  conftest.py
  test_download.py
  test_filter.py
  test_gene_mapping.py
  test_druggability.py
  test_mr_analysis.py
  test_prioritize.py

notebooks/
  01_extract_gwas_hits.ipynb
  02_map_variants_to_genes.ipynb
  03_druggability_lookup.ipynb
  04_mendelian_randomization.ipynb
  05_prioritize_targets.ipynb

data/
  raw/          # Downloaded bulk files
  processed/    # Intermediate parquet files
  results/      # Final prioritization output
```

## Testing Strategy

- **Fixtures in `conftest.py`:** sample GWAS DataFrame (~20 rows, mix of T2D and non-T2D), mock Open Targets GraphQL responses, mock OpenGWAS responses
- **Unit tests per module:** API calls mocked with `responses` library or `monkeypatch`
- **Key assertions:** correct filtering counts, gene deduplication, handling of missing data and API errors, MR calculation correctness vs. hand-computed values
- **Integration test in `test_prioritize.py`:** full scoring logic on fixture data, verifying categorization labels

## Tooling

**Dependencies to add:**
- `pytest` — test runner
- `responses` — HTTP mocking for tests
- `ruff` — linting and formatting
- `pyright` — type checking

**Ruff config:** target Python 3.14, line length 100, standard rules + import sorting

**Dev commands:**
```bash
uv run pytest                        # all tests
uv run pytest tests/test_filter.py   # single file
uv run ruff check src/ tests/        # lint
uv run ruff format src/ tests/       # format
uv run pyright src/                  # type check
```
