# T2D Drug Target Pipeline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a pipeline that identifies and prioritizes drug targets for Type 2 Diabetes using GWAS Catalog data, Open Targets druggability info, and simplified Mendelian Randomization.

**Architecture:** Core logic in `src/gwas_explorer/` modules with DataFrame-in/DataFrame-out interfaces. Notebooks in `notebooks/` orchestrate and visualize. Intermediate results cached as parquet in `data/`. Tests mock all HTTP calls.

**Tech Stack:** Python 3.14, pandas, requests, pytest, responses, ruff, pyright, statsmodels (for MR stats)

---

### Task 1: Project Scaffolding & Tooling

**Files:**
- Modify: `pyproject.toml`
- Modify: `.gitignore`
- Create: `src/gwas_explorer/__init__.py`
- Create: `src/gwas_explorer/config.py`
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`
- Create: `data/raw/.gitkeep`
- Create: `data/processed/.gitkeep`
- Create: `data/results/.gitkeep`

**Step 1: Update pyproject.toml**

Add dev dependencies and ruff config. Add the `src/gwas_explorer` package.

```toml
[project]
name = "gwas-catalogue-explorer"
version = "0.1.0"
description = "Exploring GWAS Catalogue data: associations, traits, and summary statistics"
readme = "README.md"
requires-python = ">=3.14"
dependencies = [
    "requests>=2.32.5",
    "pandas>=3.0.1",
    "numpy>=2.4.1",
    "matplotlib>=3.10.8",
    "seaborn>=0.13.2",
    "scipy>=1.17.1",
    "statsmodels>=0.14.6",
    "jupyter>=1.1.1",
    "pyarrow>=20.0.0",
]

[dependency-groups]
dev = [
    "pytest>=8.0",
    "responses>=0.25",
    "ruff>=0.11",
    "pyright>=1.1",
]

[tool.ruff]
target-version = "py314"
line-length = 100

[tool.ruff.lint]
select = ["E", "F", "I", "W", "UP"]

[tool.pytest.ini_options]
testpaths = ["tests"]

[tool.pyright]
pythonVersion = "3.14"
typeCheckingMode = "basic"
```

**Step 2: Update .gitignore**

Add data subdirectory patterns:

```gitignore
# Python-generated files
__pycache__/
*.py[oc]
build/
dist/
wheels/
*.egg-info

# Virtual environments
.venv

# Jupyter
.ipynb_checkpoints/

# Data (keep structure, ignore large files)
data/**/*.csv
data/**/*.tsv
data/**/*.gz
data/**/*.zip
data/**/*.parquet
!data/**/.gitkeep

# IDE
.vscode/
.idea/

# OS
.DS_Store
Thumbs.db
```

**Step 3: Create directory structure and config**

`src/gwas_explorer/__init__.py`:
```python
"""GWAS Catalogue Explorer — T2D drug target identification pipeline."""
```

`src/gwas_explorer/config.py`:
```python
"""Pipeline configuration: endpoints, thresholds, and paths."""

from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"

# GWAS Catalog
GWAS_CATALOG_FTP_URL = (
    "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
)
GWAS_ASSOCIATIONS_FILENAME = "gwas_catalog_associations.tsv"

# T2D filtering
T2D_EFO_TERM = "EFO_0001360"
T2D_TRAIT_KEYWORDS = ["type 2 diabetes", "type ii diabetes", "t2d"]
PVALUE_THRESHOLD = 5e-8

# Open Targets Platform
OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# OpenGWAS (IEU) API
OPENGWAS_API_URL = "https://api.opengwas.io/api"
T2D_OUTCOME_ID = "ebi-a-GCST006867"  # Large T2D GWAS outcome dataset

# MR thresholds
MR_PVALUE_THRESHOLD = 0.05
MR_INSTRUMENT_PVALUE = 5e-8

# Download settings
CACHE_STALENESS_DAYS = 30
```

`tests/__init__.py`: empty file

`tests/conftest.py`:
```python
"""Shared test fixtures for GWAS Explorer tests."""

import pandas as pd
import pytest


@pytest.fixture
def sample_gwas_associations() -> pd.DataFrame:
    """Small GWAS associations DataFrame mimicking catalog TSV structure."""
    return pd.DataFrame(
        {
            "SNPS": [
                "rs7903146",
                "rs1801282",
                "rs5219",
                "rs12255372",
                "rs999999",
                "rs111111",
                "rs7903146",
                "rs222222",
            ],
            "CHR_ID": ["10", "3", "11", "10", "1", "6", "10", "5"],
            "CHR_POS": [
                "112998590",
                "12393125",
                "17409572",
                "112998590",
                "50000000",
                "30000000",
                "112998590",
                "70000000",
            ],
            "P-VALUE": [1e-30, 2e-10, 5e-15, 3e-12, 1e-3, 1e-20, 1e-25, 1e-9],
            "OR or BETA": [1.4, 1.15, 1.2, 1.3, 1.01, 1.5, 1.35, 1.1],
            "MAPPED_GENE": [
                "TCF7L2",
                "PPARG",
                "KCNJ11",
                "TCF7L2",
                "GENEX",
                "HLA-A",
                "TCF7L2",
                "ADRB2",
            ],
            "REPORTED GENE(S)": [
                "TCF7L2",
                "PPARG",
                "KCNJ11 - ABCC8",
                "TCF7L2",
                "GENEX",
                "HLA-A",
                "TCF7L2",
                "ADRB2",
            ],
            "DISEASE/TRAIT": [
                "Type 2 diabetes",
                "Type 2 diabetes",
                "Type 2 diabetes",
                "Type 2 diabetes",
                "Height",
                "Type 2 diabetes",
                "Type 2 diabetes",
                "Type 2 diabetes",
            ],
            "MAPPED_TRAIT_URI": [
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://www.ebi.ac.uk/efo/EFO_0004339",
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://www.ebi.ac.uk/efo/EFO_0001360",
            ],
            "STUDY ACCESSION": [
                "GCST001",
                "GCST001",
                "GCST002",
                "GCST001",
                "GCST999",
                "GCST003",
                "GCST002",
                "GCST003",
            ],
        }
    )


@pytest.fixture
def sample_t2d_hits() -> pd.DataFrame:
    """Filtered T2D significant hits (output of filter module)."""
    return pd.DataFrame(
        {
            "SNP_ID": ["rs7903146", "rs1801282", "rs5219", "rs12255372", "rs111111", "rs222222"],
            "CHR": ["10", "3", "11", "10", "6", "5"],
            "POSITION": [112998590, 12393125, 17409572, 112998590, 30000000, 70000000],
            "P_VALUE": [1e-30, 2e-10, 5e-15, 3e-12, 1e-20, 1e-9],
            "OR_BETA": [1.4, 1.15, 1.2, 1.3, 1.5, 1.1],
            "MAPPED_GENE": ["TCF7L2", "PPARG", "KCNJ11", "TCF7L2", "HLA-A", "ADRB2"],
            "REPORTED_GENES": [
                "TCF7L2",
                "PPARG",
                "KCNJ11 - ABCC8",
                "TCF7L2",
                "HLA-A",
                "ADRB2",
            ],
            "STUDY_ACCESSION": [
                "GCST001",
                "GCST001",
                "GCST002",
                "GCST001",
                "GCST003",
                "GCST003",
            ],
        }
    )


@pytest.fixture
def sample_candidate_genes() -> pd.DataFrame:
    """Candidate gene list (output of gene_mapping module)."""
    return pd.DataFrame(
        {
            "GENE_SYMBOL": ["TCF7L2", "PPARG", "KCNJ11", "ABCC8", "HLA-A", "ADRB2"],
            "ENSEMBL_ID": [
                "ENSG00000148737",
                "ENSG00000132170",
                "ENSG00000187486",
                "ENSG00000006071",
                "ENSG00000206503",
                "ENSG00000169252",
            ],
            "LEAD_SNP": [
                "rs7903146",
                "rs1801282",
                "rs5219",
                "rs5219",
                "rs111111",
                "rs222222",
            ],
            "P_VALUE": [1e-30, 2e-10, 5e-15, 5e-15, 1e-20, 1e-9],
            "OR_BETA": [1.4, 1.15, 1.2, 1.2, 1.5, 1.1],
            "N_SUPPORTING_SNPS": [2, 1, 1, 1, 1, 1],
        }
    )


@pytest.fixture
def mock_open_targets_response() -> dict:
    """Mock Open Targets GraphQL response for a druggable gene (KCNJ11)."""
    return {
        "data": {
            "target": {
                "id": "ENSG00000187486",
                "approvedSymbol": "KCNJ11",
                "tractability": [
                    {
                        "modality": "SM",
                        "id": "High_Quality_ChEMBL_compounds",
                        "value": True,
                    },
                    {
                        "modality": "AB",
                        "id": "UniProt_loc_med_conf",
                        "value": False,
                    },
                ],
                "knownDrugs": {
                    "uniqueDrugs": 2,
                    "rows": [
                        {
                            "drug": {
                                "name": "GLIBENCLAMIDE",
                                "mechanismOfAction": "Potassium channel blocker",
                            },
                            "phase": 4,
                            "status": "Approved",
                            "disease": {"name": "type 2 diabetes mellitus"},
                        },
                        {
                            "drug": {
                                "name": "GLIMEPIRIDE",
                                "mechanismOfAction": "Potassium channel blocker",
                            },
                            "phase": 4,
                            "status": "Approved",
                            "disease": {"name": "type 2 diabetes mellitus"},
                        },
                    ],
                },
            }
        }
    }


@pytest.fixture
def mock_open_targets_novel_response() -> dict:
    """Mock Open Targets response for a novel gene (TCF7L2) — no drugs."""
    return {
        "data": {
            "target": {
                "id": "ENSG00000148737",
                "approvedSymbol": "TCF7L2",
                "tractability": [
                    {
                        "modality": "SM",
                        "id": "High_Quality_ChEMBL_compounds",
                        "value": False,
                    },
                ],
                "knownDrugs": None,
            }
        }
    }


@pytest.fixture
def mock_opengwas_associations() -> list[dict]:
    """Mock OpenGWAS association response for MR instruments."""
    return [
        {
            "rsid": "rs5219",
            "beta": 0.05,
            "se": 0.01,
            "p": 1e-10,
            "ea": "T",
            "nea": "C",
            "eaf": 0.35,
        },
    ]
```

**Step 4: Create data subdirectory .gitkeep files**

Create empty `.gitkeep` files in `data/raw/`, `data/processed/`, `data/results/`.

**Step 5: Install dependencies and verify**

Run:
```bash
uv sync
uv run pytest --collect-only
uv run ruff check src/ tests/
```
Expected: dependencies install, pytest discovers `tests/` (no tests yet), ruff passes.

**Step 6: Commit**

```bash
git add pyproject.toml .gitignore src/ tests/ data/raw/.gitkeep data/processed/.gitkeep data/results/.gitkeep
git commit -m "feat: project scaffolding with config, test fixtures, and dev tooling"
```

---

### Task 2: Download Module

**Files:**
- Create: `src/gwas_explorer/download.py`
- Create: `tests/test_download.py`

**Step 1: Write the failing tests**

`tests/test_download.py`:
```python
"""Tests for GWAS Catalog bulk download module."""

from pathlib import Path
from unittest.mock import patch

import responses

from gwas_explorer.config import GWAS_CATALOG_FTP_URL, RAW_DIR
from gwas_explorer.download import download_gwas_catalog


@responses.activate
def test_download_creates_file(tmp_path: Path) -> None:
    """Download should save TSV to the raw data directory."""
    fake_tsv = "col1\tcol2\nval1\tval2\n"
    responses.add(responses.GET, GWAS_CATALOG_FTP_URL, body=fake_tsv, status=200)

    with patch("gwas_explorer.download.RAW_DIR", tmp_path):
        path = download_gwas_catalog()

    assert path.exists()
    assert path.suffix == ".tsv"
    assert "val1" in path.read_text()


@responses.activate
def test_download_skips_if_recent(tmp_path: Path) -> None:
    """Download should skip if a recent file already exists."""
    existing = tmp_path / "gwas_catalog_associations.tsv"
    existing.write_text("existing data")

    with patch("gwas_explorer.download.RAW_DIR", tmp_path):
        path = download_gwas_catalog(max_age_days=999)

    assert path == existing
    assert path.read_text() == "existing data"
    assert len(responses.calls) == 0  # no HTTP request made


@responses.activate
def test_download_replaces_stale_file(tmp_path: Path) -> None:
    """Download should re-fetch if file is older than max_age_days."""
    import time

    existing = tmp_path / "gwas_catalog_associations.tsv"
    existing.write_text("old data")
    # Set mtime to 60 days ago
    old_time = time.time() - (60 * 86400)
    import os

    os.utime(existing, (old_time, old_time))

    fresh_tsv = "col1\tcol2\nfresh\tdata\n"
    responses.add(responses.GET, GWAS_CATALOG_FTP_URL, body=fresh_tsv, status=200)

    with patch("gwas_explorer.download.RAW_DIR", tmp_path):
        path = download_gwas_catalog(max_age_days=30)

    assert "fresh" in path.read_text()
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_download.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'gwas_explorer.download'`

**Step 3: Write implementation**

`src/gwas_explorer/download.py`:
```python
"""Download GWAS Catalog bulk association data."""

import time
from pathlib import Path

import requests

from gwas_explorer.config import (
    CACHE_STALENESS_DAYS,
    GWAS_ASSOCIATIONS_FILENAME,
    GWAS_CATALOG_FTP_URL,
    RAW_DIR,
)


def download_gwas_catalog(max_age_days: int = CACHE_STALENESS_DAYS) -> Path:
    """Download GWAS Catalog associations TSV, using cache if recent.

    Args:
        max_age_days: Skip download if existing file is newer than this.

    Returns:
        Path to the local TSV file.
    """
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    filepath = RAW_DIR / GWAS_ASSOCIATIONS_FILENAME

    if filepath.exists():
        age_days = (time.time() - filepath.stat().st_mtime) / 86400
        if age_days < max_age_days:
            return filepath

    response = requests.get(GWAS_CATALOG_FTP_URL, stream=True, timeout=300)
    response.raise_for_status()

    with open(filepath, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    return filepath
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_download.py -v`
Expected: 3 passed

**Step 5: Lint and commit**

```bash
uv run ruff check src/gwas_explorer/download.py tests/test_download.py
uv run ruff format src/gwas_explorer/download.py tests/test_download.py
git add src/gwas_explorer/download.py tests/test_download.py
git commit -m "feat: add GWAS Catalog bulk download with caching"
```

---

### Task 3: Filter Module

**Files:**
- Create: `src/gwas_explorer/filter.py`
- Create: `tests/test_filter.py`

**Step 1: Write the failing tests**

`tests/test_filter.py`:
```python
"""Tests for T2D association filtering module."""

import pandas as pd

from gwas_explorer.filter import filter_t2d_associations


def test_filters_by_trait(sample_gwas_associations: pd.DataFrame) -> None:
    """Should keep only T2D-related associations."""
    result = filter_t2d_associations(sample_gwas_associations)
    # "Height" row (rs999999) should be excluded
    assert "rs999999" not in result["SNP_ID"].values


def test_filters_by_pvalue(sample_gwas_associations: pd.DataFrame) -> None:
    """Should keep only genome-wide significant hits."""
    result = filter_t2d_associations(sample_gwas_associations)
    assert (result["P_VALUE"] < 5e-8).all()


def test_deduplicates_by_snp(sample_gwas_associations: pd.DataFrame) -> None:
    """Should keep one row per SNP (strongest p-value)."""
    result = filter_t2d_associations(sample_gwas_associations)
    assert result["SNP_ID"].is_unique


def test_output_schema(sample_gwas_associations: pd.DataFrame) -> None:
    """Output should have the expected columns."""
    result = filter_t2d_associations(sample_gwas_associations)
    expected_cols = {
        "SNP_ID",
        "CHR",
        "POSITION",
        "P_VALUE",
        "OR_BETA",
        "MAPPED_GENE",
        "REPORTED_GENES",
        "STUDY_ACCESSION",
    }
    assert set(result.columns) == expected_cols


def test_rs7903146_kept_with_best_pvalue(sample_gwas_associations: pd.DataFrame) -> None:
    """rs7903146 appears twice for T2D; should keep the one with p=1e-30."""
    result = filter_t2d_associations(sample_gwas_associations)
    row = result[result["SNP_ID"] == "rs7903146"]
    assert len(row) == 1
    assert row.iloc[0]["P_VALUE"] == 1e-30
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_filter.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write implementation**

`src/gwas_explorer/filter.py`:
```python
"""Filter GWAS Catalog associations for T2D genome-wide significant hits."""

import pandas as pd

from gwas_explorer.config import PVALUE_THRESHOLD, T2D_EFO_TERM, T2D_TRAIT_KEYWORDS


def filter_t2d_associations(associations: pd.DataFrame) -> pd.DataFrame:
    """Filter associations for T2D trait and genome-wide significance.

    Args:
        associations: Raw GWAS Catalog associations DataFrame.

    Returns:
        Filtered DataFrame with standardized column names.
    """
    df = associations.copy()

    # Filter by T2D trait: match EFO URI or keyword in disease/trait
    trait_col = df["DISEASE/TRAIT"].str.lower()
    uri_col = df["MAPPED_TRAIT_URI"].fillna("")

    is_t2d_keyword = trait_col.apply(
        lambda x: any(kw in x for kw in T2D_TRAIT_KEYWORDS)
    )
    is_t2d_efo = uri_col.str.contains(T2D_EFO_TERM, na=False)
    df = df[is_t2d_keyword | is_t2d_efo]

    # Filter by p-value
    df = df[df["P-VALUE"].astype(float) < PVALUE_THRESHOLD]

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
    df["P_VALUE"] = df["P_VALUE"].astype(float)
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
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_filter.py -v`
Expected: 5 passed

**Step 5: Lint and commit**

```bash
uv run ruff check src/gwas_explorer/filter.py tests/test_filter.py
uv run ruff format src/gwas_explorer/filter.py tests/test_filter.py
git add src/gwas_explorer/filter.py tests/test_filter.py
git commit -m "feat: add T2D association filtering with p-value and trait matching"
```

---

### Task 4: Gene Mapping Module

**Files:**
- Create: `src/gwas_explorer/gene_mapping.py`
- Create: `tests/test_gene_mapping.py`

**Step 1: Write the failing tests**

`tests/test_gene_mapping.py`:
```python
"""Tests for variant-to-gene mapping module."""

import pandas as pd
import responses

from gwas_explorer.gene_mapping import map_variants_to_genes


def test_splits_multi_gene_entries(sample_t2d_hits: pd.DataFrame) -> None:
    """Should split 'KCNJ11 - ABCC8' into two separate gene rows."""
    result = map_variants_to_genes(sample_t2d_hits)
    assert "KCNJ11" in result["GENE_SYMBOL"].values
    assert "ABCC8" in result["GENE_SYMBOL"].values


def test_deduplicates_genes(sample_t2d_hits: pd.DataFrame) -> None:
    """Should have one row per unique gene."""
    result = map_variants_to_genes(sample_t2d_hits)
    assert result["GENE_SYMBOL"].is_unique


def test_keeps_strongest_association(sample_t2d_hits: pd.DataFrame) -> None:
    """TCF7L2 appears multiple times; should keep lowest p-value."""
    result = map_variants_to_genes(sample_t2d_hits)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert len(tcf) == 1
    assert tcf.iloc[0]["P_VALUE"] == 1e-30


def test_counts_supporting_snps(sample_t2d_hits: pd.DataFrame) -> None:
    """TCF7L2 has 2 SNPs (rs7903146 and rs12255372)."""
    result = map_variants_to_genes(sample_t2d_hits)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert tcf.iloc[0]["N_SUPPORTING_SNPS"] == 2


def test_output_schema(sample_t2d_hits: pd.DataFrame) -> None:
    """Output should have the expected columns."""
    result = map_variants_to_genes(sample_t2d_hits)
    expected = {"GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP", "P_VALUE", "OR_BETA", "N_SUPPORTING_SNPS"}
    assert set(result.columns) == expected


@responses.activate
def test_ensembl_id_lookup(sample_t2d_hits: pd.DataFrame) -> None:
    """Should resolve Ensembl IDs via Ensembl REST API."""
    # Mock the Ensembl POST lookup endpoint
    responses.add(
        responses.POST,
        "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
        json={
            "TCF7L2": {"id": "ENSG00000148737"},
            "PPARG": {"id": "ENSG00000132170"},
            "KCNJ11": {"id": "ENSG00000187486"},
            "ABCC8": {"id": "ENSG00000006071"},
            "HLA-A": {"id": "ENSG00000206503"},
            "ADRB2": {"id": "ENSG00000169252"},
        },
        status=200,
    )

    result = map_variants_to_genes(sample_t2d_hits)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert tcf.iloc[0]["ENSEMBL_ID"] == "ENSG00000148737"
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_gene_mapping.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write implementation**

`src/gwas_explorer/gene_mapping.py`:
```python
"""Map GWAS variants to genes using catalog annotations."""

import requests
import pandas as pd


ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"


def _lookup_ensembl_ids(gene_symbols: list[str]) -> dict[str, str]:
    """Batch lookup Ensembl gene IDs for a list of gene symbols.

    Args:
        gene_symbols: List of HGNC gene symbols.

    Returns:
        Dict mapping gene symbol to Ensembl ID.
    """
    result: dict[str, str] = {}
    # Ensembl POST endpoint accepts up to 1000 symbols per request
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

    Args:
        t2d_hits: Filtered T2D associations DataFrame.

    Returns:
        DataFrame with one row per candidate gene.
    """
    rows = []
    for _, hit in t2d_hits.iterrows():
        # Prioritize MAPPED_GENE, fall back to REPORTED_GENES
        gene_str = hit["MAPPED_GENE"]
        if pd.isna(gene_str) or gene_str.strip() == "":
            gene_str = hit.get("REPORTED_GENES", "")
        if pd.isna(gene_str) or gene_str.strip() == "":
            continue

        # Split multi-gene entries (separated by " - ", ", ", or ";")
        genes = [g.strip() for g in gene_str.replace(";", " - ").replace(",", " - ").split(" - ")]
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

    # Count supporting SNPs per gene
    snp_counts = expanded.groupby("GENE_SYMBOL")["SNP_ID"].nunique().reset_index()
    snp_counts.columns = ["GENE_SYMBOL", "N_SUPPORTING_SNPS"]

    # Keep strongest association per gene
    deduped = expanded.sort_values("P_VALUE").drop_duplicates(subset=["GENE_SYMBOL"], keep="first")
    deduped = deduped.rename(columns={"SNP_ID": "LEAD_SNP"})
    deduped = deduped.merge(snp_counts, on="GENE_SYMBOL")

    # Lookup Ensembl IDs
    gene_symbols = deduped["GENE_SYMBOL"].tolist()
    ensembl_map = _lookup_ensembl_ids(gene_symbols)
    deduped["ENSEMBL_ID"] = deduped["GENE_SYMBOL"].map(ensembl_map)

    output_cols = ["GENE_SYMBOL", "ENSEMBL_ID", "LEAD_SNP", "P_VALUE", "OR_BETA", "N_SUPPORTING_SNPS"]
    return deduped[output_cols].reset_index(drop=True)
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_gene_mapping.py -v`
Expected: 6 passed

**Step 5: Lint and commit**

```bash
uv run ruff check src/gwas_explorer/gene_mapping.py tests/test_gene_mapping.py
uv run ruff format src/gwas_explorer/gene_mapping.py tests/test_gene_mapping.py
git add src/gwas_explorer/gene_mapping.py tests/test_gene_mapping.py
git commit -m "feat: add variant-to-gene mapping with Ensembl ID lookup"
```

---

### Task 5: Druggability Module

**Files:**
- Create: `src/gwas_explorer/druggability.py`
- Create: `tests/test_druggability.py`

**Step 1: Write the failing tests**

`tests/test_druggability.py`:
```python
"""Tests for Open Targets druggability lookup module."""

import pandas as pd
import responses

from gwas_explorer.config import OPEN_TARGETS_GRAPHQL_URL
from gwas_explorer.druggability import query_druggability


@responses.activate
def test_identifies_druggable_target(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """KCNJ11 should be flagged as druggable with known drugs."""
    # Mock responses for each gene
    for _, row in sample_candidate_genes.iterrows():
        if row["GENE_SYMBOL"] == "KCNJ11":
            body = mock_open_targets_response
        else:
            body = mock_open_targets_novel_response
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, json=body, status=200)

    result = query_druggability(sample_candidate_genes)
    kcnj11 = result[result["GENE_SYMBOL"] == "KCNJ11"]
    assert kcnj11.iloc[0]["IS_DRUGGABLE"] is True
    assert kcnj11.iloc[0]["HAS_EXISTING_DRUG"] is True


@responses.activate
def test_identifies_novel_target(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """TCF7L2 should show no existing drugs."""
    for _, row in sample_candidate_genes.iterrows():
        if row["GENE_SYMBOL"] == "KCNJ11":
            body = mock_open_targets_response
        else:
            body = mock_open_targets_novel_response
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, json=body, status=200)

    result = query_druggability(sample_candidate_genes)
    tcf = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert tcf.iloc[0]["HAS_EXISTING_DRUG"] is False


@responses.activate
def test_extracts_drug_details(
    sample_candidate_genes: pd.DataFrame,
    mock_open_targets_response: dict,
    mock_open_targets_novel_response: dict,
) -> None:
    """Should extract drug names and indications."""
    for _, row in sample_candidate_genes.iterrows():
        if row["GENE_SYMBOL"] == "KCNJ11":
            body = mock_open_targets_response
        else:
            body = mock_open_targets_novel_response
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, json=body, status=200)

    result = query_druggability(sample_candidate_genes)
    kcnj11 = result[result["GENE_SYMBOL"] == "KCNJ11"]
    drugs = kcnj11.iloc[0]["DRUGS"]
    assert any("GLIBENCLAMIDE" in d["name"] for d in drugs)


@responses.activate
def test_handles_api_error_gracefully(sample_candidate_genes: pd.DataFrame) -> None:
    """Should not crash on API errors; mark gene as unknown druggability."""
    for _ in range(len(sample_candidate_genes)):
        responses.add(responses.POST, OPEN_TARGETS_GRAPHQL_URL, status=500)

    result = query_druggability(sample_candidate_genes)
    assert len(result) == len(sample_candidate_genes)
    assert (result["IS_DRUGGABLE"] == False).all()  # noqa: E712
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_druggability.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write implementation**

`src/gwas_explorer/druggability.py`:
```python
"""Query Open Targets Platform for target druggability information."""

import requests
import pandas as pd

from gwas_explorer.config import OPEN_TARGETS_GRAPHQL_URL

DRUGGABILITY_QUERY = """
query TargetDruggability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability {
      modality
      id
      value
    }
    knownDrugs {
      uniqueDrugs
      rows {
        drug {
          name
          mechanismOfAction
        }
        phase
        status
        disease {
          name
        }
      }
    }
  }
}
"""


def _query_single_target(ensembl_id: str) -> dict:
    """Query Open Targets for a single gene's druggability data."""
    try:
        response = requests.post(
            OPEN_TARGETS_GRAPHQL_URL,
            json={"query": DRUGGABILITY_QUERY, "variables": {"ensemblId": ensembl_id}},
            timeout=30,
        )
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError):
        return {}


def _parse_druggability(response: dict) -> dict:
    """Parse Open Targets response into druggability summary."""
    target = response.get("data", {}).get("target")
    if not target:
        return {
            "IS_DRUGGABLE": False,
            "TRACTABILITY_SM": False,
            "TRACTABILITY_AB": False,
            "HAS_EXISTING_DRUG": False,
            "DRUGS": [],
        }

    # Parse tractability
    tractability = target.get("tractability") or []
    sm_tractable = any(t["value"] for t in tractability if t["modality"] == "SM")
    ab_tractable = any(t["value"] for t in tractability if t["modality"] == "AB")

    # Parse known drugs
    known_drugs = target.get("knownDrugs")
    drugs = []
    if known_drugs and known_drugs.get("rows"):
        for row in known_drugs["rows"]:
            drugs.append(
                {
                    "name": row["drug"]["name"],
                    "mechanism": row["drug"].get("mechanismOfAction", ""),
                    "phase": row.get("phase", 0),
                    "status": row.get("status", ""),
                    "indication": row.get("disease", {}).get("name", ""),
                }
            )

    return {
        "IS_DRUGGABLE": sm_tractable or ab_tractable or len(drugs) > 0,
        "TRACTABILITY_SM": sm_tractable,
        "TRACTABILITY_AB": ab_tractable,
        "HAS_EXISTING_DRUG": len(drugs) > 0,
        "DRUGS": drugs,
    }


def query_druggability(candidate_genes: pd.DataFrame) -> pd.DataFrame:
    """Query Open Targets for druggability of all candidate genes.

    Args:
        candidate_genes: DataFrame with GENE_SYMBOL and ENSEMBL_ID columns.

    Returns:
        Input DataFrame merged with druggability columns.
    """
    records = []
    for _, row in candidate_genes.iterrows():
        ensembl_id = row["ENSEMBL_ID"]
        response = _query_single_target(ensembl_id)
        drug_info = _parse_druggability(response)
        drug_info["GENE_SYMBOL"] = row["GENE_SYMBOL"]
        records.append(drug_info)

    drug_df = pd.DataFrame(records)
    return candidate_genes.merge(drug_df, on="GENE_SYMBOL", how="left")
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_druggability.py -v`
Expected: 4 passed

**Step 5: Lint and commit**

```bash
uv run ruff check src/gwas_explorer/druggability.py tests/test_druggability.py
uv run ruff format src/gwas_explorer/druggability.py tests/test_druggability.py
git add src/gwas_explorer/druggability.py tests/test_druggability.py
git commit -m "feat: add Open Targets druggability lookup"
```

---

### Task 6: Mendelian Randomization Module

**Files:**
- Create: `src/gwas_explorer/mr_analysis.py`
- Create: `tests/test_mr_analysis.py`

**Step 1: Write the failing tests**

`tests/test_mr_analysis.py`:
```python
"""Tests for simplified Mendelian Randomization module."""

import numpy as np
import pandas as pd
import responses

from gwas_explorer.config import OPENGWAS_API_URL
from gwas_explorer.mr_analysis import run_mr_analysis


@responses.activate
def test_wald_ratio_calculation(sample_candidate_genes: pd.DataFrame) -> None:
    """Single-SNP MR should compute correct Wald ratio."""
    # Mock: get instruments for gene (SNP-exposure)
    responses.add(
        responses.GET,
        f"{OPENGWAS_API_URL}/associations",
        json=[
            {"rsid": "rs5219", "beta": 0.3, "se": 0.05, "p": 1e-10, "ea": "T", "nea": "C", "eaf": 0.35},
        ],
        status=200,
    )
    # Mock: get outcome associations (SNP-T2D)
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=[
            {"rsid": "rs5219", "beta": 0.15, "se": 0.03, "p": 1e-6, "ea": "T", "nea": "C", "eaf": 0.35},
        ],
        status=200,
    )

    # Test with just KCNJ11
    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    # Wald ratio = beta_outcome / beta_exposure = 0.15 / 0.3 = 0.5
    assert np.isclose(row["MR_ESTIMATE"], 0.5, atol=0.01)
    assert row["METHOD"] == "wald_ratio"
    assert row["MR_STATUS"] == "ok"


@responses.activate
def test_insufficient_instruments(sample_candidate_genes: pd.DataFrame) -> None:
    """Should flag genes with no valid instruments."""
    # Mock: no instruments found
    responses.add(
        responses.GET,
        f"{OPENGWAS_API_URL}/associations",
        json=[],
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "TCF7L2"].copy()
    result = run_mr_analysis(genes)

    assert result.iloc[0]["MR_STATUS"] == "insufficient_data"
    assert np.isnan(result.iloc[0]["MR_ESTIMATE"])


@responses.activate
def test_handles_api_error(sample_candidate_genes: pd.DataFrame) -> None:
    """Should not crash on API failures."""
    responses.add(
        responses.GET,
        f"{OPENGWAS_API_URL}/associations",
        status=500,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "PPARG"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    assert result.iloc[0]["MR_STATUS"] == "error"
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_mr_analysis.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write implementation**

`src/gwas_explorer/mr_analysis.py`:
```python
"""Simplified two-sample Mendelian Randomization using OpenGWAS API."""

import math

import numpy as np
import pandas as pd
import requests

from gwas_explorer.config import (
    MR_INSTRUMENT_PVALUE,
    OPENGWAS_API_URL,
    T2D_OUTCOME_ID,
)


def _get_instruments(snp_id: str, exposure_id: str | None = None) -> list[dict]:
    """Retrieve SNP-exposure associations from OpenGWAS."""
    try:
        params = {"rsid": snp_id, "proxies": 0}
        if exposure_id:
            params["id"] = exposure_id
        response = requests.get(
            f"{OPENGWAS_API_URL}/associations",
            params=params,
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        # Filter for strong instruments
        return [d for d in data if float(d.get("p", 1)) < MR_INSTRUMENT_PVALUE]
    except (requests.RequestException, ValueError):
        return []


def _get_outcome_associations(snp_ids: list[str], outcome_id: str = T2D_OUTCOME_ID) -> dict[str, dict]:
    """Retrieve SNP-outcome associations for T2D from OpenGWAS."""
    try:
        response = requests.post(
            f"{OPENGWAS_API_URL}/associations",
            json={"rsid": snp_ids, "id": [outcome_id], "proxies": 0},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        return {d["rsid"]: d for d in data}
    except (requests.RequestException, ValueError):
        return {}


def _wald_ratio(beta_exp: float, se_exp: float, beta_out: float, se_out: float) -> dict:
    """Compute Wald ratio MR estimate from a single SNP."""
    mr_estimate = beta_out / beta_exp
    # Delta method SE
    mr_se = abs(se_out / beta_exp) * math.sqrt(1 + (se_exp / beta_exp) ** 2)
    # Z-test p-value
    z = abs(mr_estimate / mr_se) if mr_se > 0 else 0
    from scipy.stats import norm

    mr_pvalue = 2 * norm.sf(z)
    return {
        "MR_ESTIMATE": mr_estimate,
        "MR_SE": mr_se,
        "MR_PVALUE": mr_pvalue,
        "N_INSTRUMENTS": 1,
        "METHOD": "wald_ratio",
        "MR_STATUS": "ok",
    }


def _ivw(instruments: list[dict], outcomes: dict[str, dict]) -> dict:
    """Inverse-variance weighted MR estimate from multiple SNPs."""
    estimates = []
    weights = []

    for inst in instruments:
        rsid = inst["rsid"]
        if rsid not in outcomes:
            continue
        out = outcomes[rsid]
        beta_exp = float(inst["beta"])
        se_exp = float(inst["se"])
        beta_out = float(out["beta"])
        se_out = float(out["se"])

        if abs(beta_exp) < 1e-10:
            continue

        ratio = beta_out / beta_exp
        ratio_se = abs(se_out / beta_exp)
        if ratio_se > 0:
            w = 1 / (ratio_se**2)
            estimates.append(ratio)
            weights.append(w)

    if not estimates:
        return {
            "MR_ESTIMATE": np.nan,
            "MR_SE": np.nan,
            "MR_PVALUE": np.nan,
            "N_INSTRUMENTS": 0,
            "METHOD": "ivw",
            "MR_STATUS": "insufficient_data",
        }

    weights_arr = np.array(weights)
    estimates_arr = np.array(estimates)
    ivw_estimate = np.sum(weights_arr * estimates_arr) / np.sum(weights_arr)
    ivw_se = math.sqrt(1 / np.sum(weights_arr))

    from scipy.stats import norm

    z = abs(ivw_estimate / ivw_se) if ivw_se > 0 else 0
    mr_pvalue = 2 * norm.sf(z)

    return {
        "MR_ESTIMATE": float(ivw_estimate),
        "MR_SE": float(ivw_se),
        "MR_PVALUE": float(mr_pvalue),
        "N_INSTRUMENTS": len(estimates),
        "METHOD": "ivw",
        "MR_STATUS": "ok",
    }


def run_mr_analysis(candidate_genes: pd.DataFrame) -> pd.DataFrame:
    """Run simplified MR for each candidate gene.

    Args:
        candidate_genes: DataFrame with GENE_SYMBOL and LEAD_SNP columns.

    Returns:
        DataFrame with MR results per gene.
    """
    results = []
    for _, row in candidate_genes.iterrows():
        gene = row["GENE_SYMBOL"]
        lead_snp = row["LEAD_SNP"]

        instruments = _get_instruments(lead_snp)

        if not instruments:
            results.append(
                {
                    "GENE_SYMBOL": gene,
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "insufficient_data",
                }
            )
            continue

        snp_ids = [inst["rsid"] for inst in instruments]
        outcomes = _get_outcome_associations(snp_ids)

        if not outcomes:
            results.append(
                {
                    "GENE_SYMBOL": gene,
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "insufficient_data",
                }
            )
            continue

        if len(instruments) == 1:
            inst = instruments[0]
            rsid = inst["rsid"]
            if rsid in outcomes:
                out = outcomes[rsid]
                mr = _wald_ratio(
                    float(inst["beta"]),
                    float(inst["se"]),
                    float(out["beta"]),
                    float(out["se"]),
                )
            else:
                mr = {
                    "MR_ESTIMATE": np.nan,
                    "MR_SE": np.nan,
                    "MR_PVALUE": np.nan,
                    "N_INSTRUMENTS": 0,
                    "METHOD": "none",
                    "MR_STATUS": "insufficient_data",
                }
        else:
            mr = _ivw(instruments, outcomes)

        mr["GENE_SYMBOL"] = gene
        results.append(mr)

    return pd.DataFrame(results)[
        ["GENE_SYMBOL", "MR_ESTIMATE", "MR_SE", "MR_PVALUE", "N_INSTRUMENTS", "METHOD", "MR_STATUS"]
    ].reset_index(drop=True)
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_mr_analysis.py -v`
Expected: 3 passed

**Step 5: Lint and commit**

```bash
uv run ruff check src/gwas_explorer/mr_analysis.py tests/test_mr_analysis.py
uv run ruff format src/gwas_explorer/mr_analysis.py tests/test_mr_analysis.py
git add src/gwas_explorer/mr_analysis.py tests/test_mr_analysis.py
git commit -m "feat: add simplified Mendelian Randomization analysis"
```

---

### Task 7: Prioritization Module

**Files:**
- Create: `src/gwas_explorer/prioritize.py`
- Create: `tests/test_prioritize.py`

**Step 1: Write the failing tests**

`tests/test_prioritize.py`:
```python
"""Tests for target prioritization and scoring module."""

import numpy as np
import pandas as pd

from gwas_explorer.prioritize import prioritize_targets


@pytest.fixture
def full_pipeline_data() -> pd.DataFrame:
    """Combined data simulating full pipeline output for prioritization."""
    return pd.DataFrame(
        {
            "GENE_SYMBOL": ["KCNJ11", "TCF7L2", "PPARG", "ADRB2", "HLA-A"],
            "ENSEMBL_ID": [
                "ENSG00000187486",
                "ENSG00000148737",
                "ENSG00000132170",
                "ENSG00000169252",
                "ENSG00000206503",
            ],
            "LEAD_SNP": ["rs5219", "rs7903146", "rs1801282", "rs222222", "rs111111"],
            "P_VALUE": [5e-15, 1e-30, 2e-10, 1e-9, 1e-20],
            "OR_BETA": [1.2, 1.4, 1.15, 1.1, 1.5],
            "N_SUPPORTING_SNPS": [1, 2, 1, 1, 1],
            "IS_DRUGGABLE": [True, False, True, True, False],
            "TRACTABILITY_SM": [True, False, True, True, False],
            "TRACTABILITY_AB": [False, False, False, False, False],
            "HAS_EXISTING_DRUG": [True, False, True, True, False],
            "DRUGS": [
                [{"name": "GLIBENCLAMIDE", "indication": "type 2 diabetes mellitus", "phase": 4, "mechanism": "", "status": "Approved"}],
                [],
                [{"name": "PIOGLITAZONE", "indication": "type 2 diabetes mellitus", "phase": 4, "mechanism": "", "status": "Approved"}],
                [{"name": "SALBUTAMOL", "indication": "asthma", "phase": 4, "mechanism": "", "status": "Approved"}],
                [],
            ],
            "MR_ESTIMATE": [0.5, np.nan, 0.3, 0.1, np.nan],
            "MR_SE": [0.1, np.nan, 0.08, 0.2, np.nan],
            "MR_PVALUE": [1e-6, np.nan, 0.001, 0.6, np.nan],
            "N_INSTRUMENTS": [1, 0, 1, 1, 0],
            "METHOD": ["wald_ratio", "none", "wald_ratio", "wald_ratio", "none"],
            "MR_STATUS": ["ok", "insufficient_data", "ok", "ok", "error"],
        }
    )


import pytest


def test_validated_target(full_pipeline_data: pd.DataFrame) -> None:
    """KCNJ11 has T2D drug → should be 'validated'."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "KCNJ11"]
    assert row.iloc[0]["TARGET_CATEGORY"] == "validated"


def test_repurposing_candidate(full_pipeline_data: pd.DataFrame) -> None:
    """ADRB2 has drug but not for T2D → should be 'repurposing_candidate'."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "ADRB2"]
    assert row.iloc[0]["TARGET_CATEGORY"] == "repurposing_candidate"


def test_novel_challenging(full_pipeline_data: pd.DataFrame) -> None:
    """TCF7L2 has no drugs and is not druggable → should be 'novel_challenging'."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "TCF7L2"]
    assert row.iloc[0]["TARGET_CATEGORY"] == "novel_challenging"


def test_mr_causality_flag(full_pipeline_data: pd.DataFrame) -> None:
    """KCNJ11 MR p < 0.05 → MR_SUPPORTS_CAUSALITY should be True."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "KCNJ11"]
    assert row.iloc[0]["MR_SUPPORTS_CAUSALITY"] is True


def test_mr_no_causality_flag(full_pipeline_data: pd.DataFrame) -> None:
    """ADRB2 MR p > 0.05 → MR_SUPPORTS_CAUSALITY should be False."""
    result = prioritize_targets(full_pipeline_data)
    row = result[result["GENE_SYMBOL"] == "ADRB2"]
    assert row.iloc[0]["MR_SUPPORTS_CAUSALITY"] is False


def test_output_sorted_by_category_and_pvalue(full_pipeline_data: pd.DataFrame) -> None:
    """Output should be sorted: validated first, then by p-value within category."""
    result = prioritize_targets(full_pipeline_data)
    categories = result["TARGET_CATEGORY"].tolist()
    # validated should come before repurposing_candidate before novel
    category_order = ["validated", "repurposing_candidate", "novel_druggable", "novel_challenging"]
    cat_indices = [category_order.index(c) for c in categories]
    assert cat_indices == sorted(cat_indices)
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_prioritize.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write implementation**

`src/gwas_explorer/prioritize.py`:
```python
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
    """Score and categorize drug target candidates.

    Args:
        combined_data: Merged DataFrame with gene info, druggability, and MR results.

    Returns:
        Prioritized DataFrame sorted by category and evidence strength.
    """
    df = combined_data.copy()

    # Derive flags
    df["HAS_GWAS_SIGNAL"] = True
    df["DRUG_INDICATION_MATCHES_T2D"] = df["DRUGS"].apply(
        lambda drugs: _matches_t2d_indication(drugs) if isinstance(drugs, list) else False
    )
    df["MR_SUPPORTS_CAUSALITY"] = df["MR_PVALUE"].apply(
        lambda p: bool(p < MR_PVALUE_THRESHOLD) if pd.notna(p) else False
    )

    # Assign category
    df["TARGET_CATEGORY"] = df.apply(_assign_category, axis=1)

    # Sort: by category order, then by p-value within category
    df["_cat_rank"] = df["TARGET_CATEGORY"].map({c: i for i, c in enumerate(CATEGORY_ORDER)})
    df = df.sort_values(["_cat_rank", "P_VALUE"]).drop(columns=["_cat_rank"])

    return df.reset_index(drop=True)
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_prioritize.py -v`
Expected: 6 passed

**Step 5: Lint and commit**

```bash
uv run ruff check src/gwas_explorer/prioritize.py tests/test_prioritize.py
uv run ruff format src/gwas_explorer/prioritize.py tests/test_prioritize.py
git add src/gwas_explorer/prioritize.py tests/test_prioritize.py
git commit -m "feat: add target prioritization and categorization"
```

---

### Task 8: Notebooks

**Files:**
- Create: `notebooks/01_extract_gwas_hits.ipynb`
- Create: `notebooks/02_map_variants_to_genes.ipynb`
- Create: `notebooks/03_druggability_lookup.ipynb`
- Create: `notebooks/04_mendelian_randomization.ipynb`
- Create: `notebooks/05_prioritize_targets.ipynb`

Each notebook follows the same pattern:

1. **Imports** — `from gwas_explorer import ...`
2. **Run step** — call the module function
3. **Inspect** — display DataFrame head, shape, describe
4. **Visualize** — relevant plots (e.g., Manhattan-style p-value plot, bar chart of druggability categories, volcano plot of MR results)
5. **Save** — write parquet to `data/processed/` or `data/results/`

**Notebook 01: Extract GWAS Hits**
```python
# Cell 1: Imports
from gwas_explorer.download import download_gwas_catalog
from gwas_explorer.filter import filter_t2d_associations
import pandas as pd

# Cell 2: Download
catalog_path = download_gwas_catalog()
print(f"Downloaded to: {catalog_path}")

# Cell 3: Load and filter
raw = pd.read_csv(catalog_path, sep="\t", low_memory=False)
print(f"Raw associations: {len(raw):,}")

t2d_hits = filter_t2d_associations(raw)
print(f"T2D significant hits: {len(t2d_hits)}")
t2d_hits.head(10)

# Cell 4: Visualize p-value distribution
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(range(len(t2d_hits)), -np.log10(t2d_hits["P_VALUE"].sort_values()), alpha=0.6)
ax.axhline(-np.log10(5e-8), color="red", linestyle="--", label="p = 5e-8")
ax.set_xlabel("SNP index")
ax.set_ylabel("-log10(p-value)")
ax.set_title("T2D GWAS Hits")
ax.legend()
plt.tight_layout()

# Cell 5: Save
from gwas_explorer.config import PROCESSED_DIR
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
t2d_hits.to_parquet(PROCESSED_DIR / "t2d_significant_hits.parquet", index=False)
print(f"Saved {len(t2d_hits)} hits")
```

**Notebook 02: Map Variants to Genes**
```python
# Cell 1: Imports
from gwas_explorer.gene_mapping import map_variants_to_genes
from gwas_explorer.config import PROCESSED_DIR
import pandas as pd

# Cell 2: Load and map
t2d_hits = pd.read_parquet(PROCESSED_DIR / "t2d_significant_hits.parquet")
candidate_genes = map_variants_to_genes(t2d_hits)
print(f"Candidate genes: {len(candidate_genes)}")
candidate_genes.head(10)

# Cell 3: Bar chart of supporting SNPs per gene (top 20)
import matplotlib.pyplot as plt

top = candidate_genes.nlargest(20, "N_SUPPORTING_SNPS")
fig, ax = plt.subplots(figsize=(10, 6))
ax.barh(top["GENE_SYMBOL"], top["N_SUPPORTING_SNPS"])
ax.set_xlabel("Number of supporting SNPs")
ax.set_title("Top 20 T2D Candidate Genes by Supporting Evidence")
plt.tight_layout()

# Cell 4: Save
candidate_genes.to_parquet(PROCESSED_DIR / "t2d_candidate_genes.parquet", index=False)
print(f"Saved {len(candidate_genes)} candidate genes")
```

**Notebook 03: Druggability Lookup**
```python
# Cell 1: Imports
from gwas_explorer.druggability import query_druggability
from gwas_explorer.config import PROCESSED_DIR
import pandas as pd

# Cell 2: Load and query
candidate_genes = pd.read_parquet(PROCESSED_DIR / "t2d_candidate_genes.parquet")
druggability = query_druggability(candidate_genes)
print(f"Druggable: {druggability['IS_DRUGGABLE'].sum()} / {len(druggability)}")
print(f"With existing drugs: {druggability['HAS_EXISTING_DRUG'].sum()}")
druggability.head(10)

# Cell 3: Pie chart of druggability
import matplotlib.pyplot as plt

counts = druggability["IS_DRUGGABLE"].value_counts()
fig, ax = plt.subplots()
ax.pie(counts, labels=["Not druggable", "Druggable"], autopct="%1.0f%%")
ax.set_title("Druggability of T2D Candidate Genes")

# Cell 4: Save
druggability.to_parquet(PROCESSED_DIR / "t2d_druggability.parquet", index=False)
print(f"Saved druggability data for {len(druggability)} genes")
```

**Notebook 04: Mendelian Randomization**
```python
# Cell 1: Imports
from gwas_explorer.mr_analysis import run_mr_analysis
from gwas_explorer.config import PROCESSED_DIR
import pandas as pd

# Cell 2: Load and run MR
candidate_genes = pd.read_parquet(PROCESSED_DIR / "t2d_candidate_genes.parquet")
mr_results = run_mr_analysis(candidate_genes)
status_counts = mr_results["MR_STATUS"].value_counts()
print(f"MR results:\n{status_counts}")
mr_results.head(10)

# Cell 3: Volcano plot of MR results
import matplotlib.pyplot as plt
import numpy as np

ok = mr_results[mr_results["MR_STATUS"] == "ok"]
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(ok["MR_ESTIMATE"], -np.log10(ok["MR_PVALUE"]), alpha=0.6)
ax.axhline(-np.log10(0.05), color="red", linestyle="--", label="p = 0.05")
ax.set_xlabel("MR Estimate (causal effect)")
ax.set_ylabel("-log10(p-value)")
ax.set_title("Mendelian Randomization Results")
ax.legend()
plt.tight_layout()

# Cell 4: Save
mr_results.to_parquet(PROCESSED_DIR / "t2d_mr_results.parquet", index=False)
print(f"Saved MR results for {len(mr_results)} genes")
```

**Notebook 05: Prioritize Targets**
```python
# Cell 1: Imports
from gwas_explorer.prioritize import prioritize_targets
from gwas_explorer.config import PROCESSED_DIR, RESULTS_DIR
import pandas as pd

# Cell 2: Load and merge all data
druggability = pd.read_parquet(PROCESSED_DIR / "t2d_druggability.parquet")
mr_results = pd.read_parquet(PROCESSED_DIR / "t2d_mr_results.parquet")
combined = druggability.merge(mr_results, on="GENE_SYMBOL", how="left")

# Cell 3: Prioritize
prioritized = prioritize_targets(combined)
print("Target categories:")
print(prioritized["TARGET_CATEGORY"].value_counts())
prioritized.head(20)

# Cell 4: Summary visualization
import matplotlib.pyplot as plt

cat_counts = prioritized["TARGET_CATEGORY"].value_counts()
colors = {"validated": "#2ecc71", "repurposing_candidate": "#3498db", "novel_druggable": "#f39c12", "novel_challenging": "#e74c3c"}
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(cat_counts.index, cat_counts.values, color=[colors.get(c, "gray") for c in cat_counts.index])
ax.set_ylabel("Number of genes")
ax.set_title("T2D Drug Target Prioritization")
plt.xticks(rotation=15)
plt.tight_layout()

# Cell 5: Save results
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
prioritized.to_parquet(RESULTS_DIR / "t2d_target_prioritization.parquet", index=False)
# Also save CSV (exclude DRUGS list column for readability)
csv_cols = [c for c in prioritized.columns if c != "DRUGS"]
prioritized[csv_cols].to_csv(RESULTS_DIR / "t2d_target_prioritization.csv", index=False)
print(f"Saved {len(prioritized)} prioritized targets")

# Cell 6: Show top novel druggable targets
novel = prioritized[prioritized["TARGET_CATEGORY"] == "novel_druggable"]
print(f"\n=== Top Novel Druggable Targets ({len(novel)}) ===")
novel[["GENE_SYMBOL", "P_VALUE", "MR_SUPPORTS_CAUSALITY", "TRACTABILITY_SM"]].head(15)
```

**Step 1: Create all notebooks**

Use `uv run jupyter` or write them as .ipynb JSON files.

**Step 2: Verify notebooks load**

Run: `uv run jupyter nbconvert --to script notebooks/01_extract_gwas_hits.ipynb --stdout | head`
Expected: converted Python output

**Step 3: Commit**

```bash
git add notebooks/
git commit -m "feat: add orchestration notebooks for full T2D pipeline"
```

---

### Task 9: Update CLAUDE.md and main.py

**Files:**
- Modify: `CLAUDE.md`
- Modify: `main.py`

**Step 1: Update CLAUDE.md** with new dev commands, architecture, and module descriptions.

**Step 2: Update main.py** to serve as a CLI entry point that runs the full pipeline:

```python
"""Run the full T2D drug target identification pipeline."""

from gwas_explorer.config import PROCESSED_DIR, RESULTS_DIR
from gwas_explorer.download import download_gwas_catalog
from gwas_explorer.druggability import query_druggability
from gwas_explorer.filter import filter_t2d_associations
from gwas_explorer.gene_mapping import map_variants_to_genes
from gwas_explorer.mr_analysis import run_mr_analysis
from gwas_explorer.prioritize import prioritize_targets

import pandas as pd


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
```

**Step 3: Run full test suite**

Run: `uv run pytest -v`
Expected: all tests pass

**Step 4: Lint everything**

Run: `uv run ruff check src/ tests/ main.py && uv run ruff format src/ tests/ main.py`

**Step 5: Commit**

```bash
git add CLAUDE.md main.py
git commit -m "feat: add CLI entry point and update project docs"
```

---

### Task 10: Final Verification

**Step 1: Run full test suite with coverage**

Run: `uv run pytest -v --tb=short`
Expected: all tests pass

**Step 2: Lint and type check**

Run:
```bash
uv run ruff check src/ tests/ main.py
uv run pyright src/
```
Expected: no errors

**Step 3: Final commit (if any fixes needed)**

```bash
git add -A
git commit -m "chore: final cleanup and verification"
```
