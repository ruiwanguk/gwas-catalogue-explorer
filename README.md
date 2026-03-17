# GWAS Catalogue Explorer

**A demonstration pipeline for identifying potential drug targets for Type 2 Diabetes using publicly available genetic data.**

> **Note:** This project is for **educational and demonstration purposes only**. It illustrates how GWAS data can be used in a drug target identification workflow, but uses simplified methods and lacks the statistical rigor required for real-world research or clinical decision-making. Do not use these results to draw biological or therapeutic conclusions.

This pipeline takes genome-wide association study (GWAS) data and works through a systematic process: starting from thousands of genetic variants linked to Type 2 Diabetes (T2D), it narrows them down to a prioritized list of potential drug targets — ranked by whether they are already druggable, already targeted by existing drugs, and whether there is causal evidence linking them to the disease.

## Background

### What problem does this solve?

Identifying new drug targets is one of the biggest bottlenecks in drug development. Genetics can help: if a genetic variant is strongly associated with a disease, the gene it affects might be a good drug target. But going from a list of thousands of GWAS "hits" to a shortlist of actionable targets requires several steps — mapping variants to genes, checking if those genes encode proteins we can actually drug, and gathering causal evidence. This pipeline automates that process for T2D.

### Key concepts

- **GWAS (Genome-Wide Association Study)** — Large-scale studies that scan the genome for genetic variants statistically associated with a disease or trait. A "hit" means a variant occurs more often in people with the disease than without.
- **Genome-wide significance** — A strict statistical threshold (p < 5 × 10⁻⁸) used to filter out false positives, given the millions of variants tested.
- **Druggability / Tractability** — Whether a protein can be targeted by a drug. Some proteins have binding pockets suitable for small molecules; others can be reached by antibodies. Many cannot be easily drugged with current technology.
- **Mendelian Randomization (MR)** — A statistical method that uses genetic variants as "natural experiments" to test whether a gene's activity causally affects disease risk, rather than just being correlated with it.

## What the pipeline does

The pipeline runs six steps, each building on the last:

```
Step 1: Download       Fetch the full GWAS Catalog (~300k associations) from EBI
Step 2: Filter         Keep only T2D-associated variants at genome-wide significance
Step 3: Gene mapping   Map variants to their nearest genes and resolve gene IDs
Step 4: Druggability   Check each gene against Open Targets for drug tractability
Step 5: MR analysis    Run two-sample Mendelian Randomization for causal evidence
Step 6: Prioritize     Assign each target to one of four categories
```

### Target categories

Each gene ends up in one of these categories:

| Category | What it means | Example |
|---|---|---|
| **Validated** | A drug already exists for this gene *and* is used to treat T2D | KCNJ11 (targeted by sulfonylureas) |
| **Repurposing candidate** | A drug exists but is approved for a different disease — it could potentially be repurposed for T2D | A gene targeted by an oncology drug |
| **Novel druggable** | No drug exists yet, but the protein is tractable (small molecule or antibody) | A newly identified T2D gene with a good binding pocket |
| **Novel challenging** | No drug exists and the protein is hard to drug with current approaches | An intracellular protein with no known binding site |

## Getting started

### Prerequisites

- Python 3.14+
- [uv](https://docs.astral.sh/uv/) (Python package manager)

### Installation

```bash
git clone https://github.com/ruiwanguk/gwas-catalogue-explorer.git
cd gwas-catalogue-explorer
uv sync
```

### Run the full pipeline

```bash
uv run python main.py
```

This takes a few minutes (mostly waiting on API calls). Results are saved to `data/results/`.

Options:
```bash
uv run python main.py --max-workers 8      # More concurrent API requests (default: 4)
uv run python main.py --max-age-days 7      # Re-download if data is older than 7 days
uv run python main.py --log-level DEBUG     # Verbose logging
```

### Run step-by-step with notebooks

If you prefer to explore interactively, five Jupyter notebooks walk through each pipeline step with visualizations:

```bash
uv run jupyter notebook
```

| Notebook | Step |
|---|---|
| `01_extract_gwas_hits.ipynb` | Download and filter GWAS data |
| `02_map_variants_to_genes.ipynb` | Map variants to genes |
| `03_druggability_lookup.ipynb` | Assess druggability via Open Targets |
| `04_mendelian_randomization.ipynb` | Run Mendelian Randomization |
| `05_prioritize_targets.ipynb` | Score and categorize targets |

## Project structure

```
main.py                          CLI entry point — runs the full pipeline

src/gwas_explorer/
  config.py                      API endpoints, thresholds, file paths
  http_utils.py                  Retry-capable HTTP session and rate-limited concurrency
  download.py                    GWAS Catalog bulk download with caching
  filter.py                      T2D trait and significance filtering
  gene_mapping.py                Variant-to-gene mapping with Ensembl ID resolution
  druggability.py                Open Targets druggability and known drug queries
  mr_analysis.py                 Two-sample Mendelian Randomization (Wald ratio / IVW)
  prioritize.py                  Target scoring and four-tier categorization

notebooks/                       Step-by-step Jupyter notebooks with visualizations
tests/                           pytest test suite (fully mocked, no network access)
data/
  raw/                           Downloaded bulk files
  processed/                     Intermediate results (parquet)
  results/                       Final prioritized target list (parquet + CSV)
```

## External data sources

This pipeline pulls data from four public APIs — no API keys required:

| Service | What it provides |
|---|---|
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) (EBI) | Curated catalog of genome-wide association studies |
| [Ensembl REST API](https://rest.ensembl.org) | Gene symbol to Ensembl ID mapping |
| [Open Targets Platform](https://platform.opentargets.org/) | Druggability/tractability assessments and known drug data |
| [OpenGWAS](https://gwas.mrcieu.ac.uk/) (MRC IEU) | GWAS summary statistics for Mendelian Randomization |

## Development

```bash
uv run pytest                           # Run tests
uv run ruff check src/ tests/           # Lint
uv run ruff format src/ tests/          # Format
uv run pyright src/                     # Type check
```

## Limitations

**This is an educational project, not a production-grade genetics platform.** The results should not be used for research conclusions, clinical decisions, or drug development without substantially more rigorous analysis. Key limitations:

- **No LD clumping** — Nearby correlated variants may inflate candidate counts
- **Simplified MR** — Uses Wald ratio and IVW only; no pleiotropy-robust methods (MR-Egger, MR-PRESSO)
- **No colocalization** — Cannot distinguish shared vs. distinct causal variants between GWAS and eQTL signals
- **T2D only** — Hardcoded for Type 2 Diabetes; adapting to other diseases requires code changes
- **No tissue-specific filtering** — Does not filter by expression in disease-relevant tissues (pancreas, liver, adipose)

## Disclaimer

This project is provided as-is for **educational and demonstration purposes only**. It is not intended for production use, clinical decision-making, or as a basis for therapeutic development. The methods used are simplified and do not reflect the full rigor of professional genetic epidemiology or drug discovery workflows.
