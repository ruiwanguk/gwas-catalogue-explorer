# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

T2D drug target identification and validation pipeline. Uses GWAS Catalog data to identify genome-wide significant variants, maps them to genes, assesses druggability via Open Targets, runs simplified Mendelian Randomization for causal evidence, and prioritizes targets into categories: validated, repurposing candidate, novel druggable, and novel challenging.

## Development Commands

This project uses `uv` as the package manager with Python 3.14.

```bash
# Install all dependencies (including dev)
uv sync

# Run the full pipeline
uv run python main.py

# Run Jupyter notebooks
uv run jupyter notebook

# Run tests
uv run pytest                        # all tests
uv run pytest tests/test_filter.py   # single file
uv run pytest -k "test_name"         # single test

# Lint and format
uv run ruff check src/ tests/
uv run ruff format src/ tests/

# Type check
uv run pyright src/
```

## Architecture

Core logic lives in `src/gwas_explorer/` modules with DataFrame-in/DataFrame-out interfaces. Notebooks in `notebooks/` orchestrate and visualize each step. Intermediate results are cached as parquet files in `data/`.

**Pipeline flow:**
```
download.py → filter.py → gene_mapping.py → druggability.py → mr_analysis.py → prioritize.py
```

- `config.py` — API endpoints, thresholds (p-value, EFO terms), file paths
- `http_utils.py` — Retry-capable requests session and `RateLimitedExecutor` for concurrent API calls
- `download.py` — Bulk download of GWAS Catalog associations TSV with staleness-based caching
- `filter.py` — Filters for T2D trait (EFO_0001360 + keywords) and genome-wide significance (p < 5e-8)
- `gene_mapping.py` — Splits multi-gene entries, deduplicates, resolves Ensembl IDs via REST API
- `druggability.py` — Queries Open Targets GraphQL for tractability and known drugs per gene
- `mr_analysis.py` — Simplified two-sample MR (Wald ratio / IVW) using OpenGWAS API
- `prioritize.py` — Categorizes targets: validated, repurposing_candidate, novel_druggable, novel_challenging

## Testing

Tests use `pytest` with `responses` for HTTP mocking. Shared fixtures in `tests/conftest.py` provide sample DataFrames and mock API responses. All API calls are mocked — no network access during tests.

## Key Resources

- GWAS Catalogue REST API: https://www.ebi.ac.uk/gwas/rest/api
- Open Targets GraphQL: https://api.platform.opentargets.org/api/v4/graphql
- OpenGWAS API: https://api.opengwas.io/api
- Ensembl REST API: https://rest.ensembl.org
