# eQTL-Based Two-Sample MR Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the broken instrument lookup in MR analysis with genuine two-sample MR using eQTLGen blood eQTL data as the exposure and T2D GWAS as the outcome.

**Architecture:** Use the OpenGWAS `/tophits` endpoint to retrieve genome-wide significant eQTL instruments for each gene via its Ensembl ID (study ID pattern: `eqtl-a-{ENSEMBL_ID}`). Then look up those SNPs in the T2D outcome study. The existing IVW/Wald ratio math is correct and stays unchanged. `_run_mr_for_gene` gains the Ensembl ID argument and orchestrates the new flow.

**Tech Stack:** Python, pandas, requests, OpenGWAS REST API, pytest + responses

---

### Task 1: Add eQTLGen config constant

**Files:**
- Modify: `src/gwas_explorer/config.py:30-32`

**Step 1: Add the constant**

Add after line 32 (`T2D_OUTCOME_ID`):

```python
EQTLGEN_STUDY_PREFIX = "eqtl-a-"  # eQTLGen blood cis-eQTL studies on OpenGWAS
```

**Step 2: Commit**

```bash
git add src/gwas_explorer/config.py
git commit -m "feat(config): add eQTLGen study prefix constant"
```

---

### Task 2: Replace `_get_instruments` with eQTL tophits lookup

**Files:**
- Modify: `src/gwas_explorer/mr_analysis.py`

**Step 1: Write the failing test**

In `tests/test_mr_analysis.py`, add a new test that exercises the new eQTL tophits flow. The key change: the instrument lookup now calls `/tophits` (not `/associations`), and the outcome lookup calls `/associations` with the T2D study ID.

```python
@responses.activate
def test_eqtl_instruments_from_tophits(sample_candidate_genes: pd.DataFrame) -> None:
    """Should fetch eQTL instruments via /tophits and outcome via /associations."""
    eqtl_instruments = [
        {"rsid": "rs10885396", "beta": -0.285, "se": 0.01, "p": 8.3e-132,
         "ea": "A", "nea": "G", "eaf": 0.4, "id": "eqtl-a-ENSG00000187486",
         "trait": "ENSG00000187486"},
    ]
    outcome_associations = [
        {"rsid": "rs10885396", "beta": 0.02, "se": 0.008, "p": 0.01,
         "ea": "A", "nea": "G", "eaf": 0.4, "id": "ebi-a-GCST006867"},
    ]

    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=eqtl_instruments,
        status=200,
    )
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/associations",
        json=outcome_associations,
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "KCNJ11"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["MR_STATUS"] == "ok"
    assert row["METHOD"] == "wald_ratio"
    assert row["N_INSTRUMENTS"] == 1
    # Wald ratio = 0.02 / -0.285 ≈ -0.0702
    assert abs(row["MR_ESTIMATE"] - (0.02 / -0.285)) < 0.01
```

**Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_mr_analysis.py::test_eqtl_instruments_from_tophits -v`
Expected: FAIL (test doesn't exist yet — add it, then it fails because code still uses old path)

**Step 3: Rewrite `_get_instruments` and update `_run_mr_for_gene`**

Replace `_get_instruments` (lines 41-60) with:

```python
def _get_eqtl_instruments(ensembl_id: str) -> list[dict] | None:
    """Retrieve genome-wide significant cis-eQTL instruments from eQTLGen via OpenGWAS.

    Uses the /tophits endpoint to get the strongest eQTL SNPs for the gene.
    Returns a list of instruments, or None if the API request failed.
    """
    eqtl_study_id = f"{EQTLGEN_STUDY_PREFIX}{ensembl_id}"
    try:
        response = _session.post(
            f"{OPENGWAS_API_URL}/tophits",
            json={"id": [eqtl_study_id]},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        if not isinstance(data, list):
            return None
        return [d for d in data if float(d.get("p", 1)) < MR_INSTRUMENT_PVALUE]
    except Exception:
        logger.warning("Failed to get eQTL instruments for %s (%s)", ensembl_id, eqtl_study_id)
        return None
```

Update the import in `mr_analysis.py` to include `EQTLGEN_STUDY_PREFIX`:

```python
from gwas_explorer.config import (
    EQTLGEN_STUDY_PREFIX,
    MR_INSTRUMENT_PVALUE,
    OPENGWAS_API_URL,
    OPENGWAS_TOKEN,
    T2D_OUTCOME_ID,
)
```

Update `_run_mr_for_gene` signature and body (lines 194-230):

```python
def _run_mr_for_gene(gene: str, ensembl_id: str, lead_snp: str) -> dict:
    """Run two-sample MR for a single gene using eQTL exposure and T2D outcome."""
    instruments = _get_eqtl_instruments(ensembl_id)

    if instruments is None:
        return {"GENE_SYMBOL": gene, **_INSUFFICIENT, "MR_STATUS": "error"}

    if not instruments:
        return {"GENE_SYMBOL": gene, **_INSUFFICIENT, "MR_STATUS": "no_eqtl_instruments"}

    snp_ids = [inst["rsid"] for inst in instruments]
    outcomes = _get_outcome_associations(snp_ids)

    if not outcomes:
        return {"GENE_SYMBOL": gene, **_INSUFFICIENT}

    if len(instruments) == 1:
        inst = instruments[0]
        rsid = inst["rsid"]
        if rsid in outcomes:
            out = _harmonize_alleles(inst, outcomes[rsid])
            if out is None:
                mr = dict(_INSUFFICIENT)
            else:
                mr = _wald_ratio(
                    float(inst["beta"]),
                    float(inst["se"]),
                    float(out["beta"]),
                    float(out["se"]),
                )
        else:
            mr = dict(_INSUFFICIENT)
    else:
        mr = _ivw(instruments, outcomes)

    mr["GENE_SYMBOL"] = gene
    return mr
```

Update `run_mr_analysis` to pass `ensembl_id` (lines 233-255):

```python
def run_mr_analysis(candidate_genes: pd.DataFrame, max_workers: int = 4) -> pd.DataFrame:
    """Run two-sample MR (eQTL exposure → T2D outcome) for each candidate gene.

    Args:
        candidate_genes: DataFrame with GENE_SYMBOL, ENSEMBL_ID, and LEAD_SNP columns.
        max_workers: Number of concurrent API requests (default 4).

    Returns:
        DataFrame with MR results per gene.
    """
    if candidate_genes.empty:
        return pd.DataFrame(
            columns=["GENE_SYMBOL", "MR_ESTIMATE", "MR_SE", "MR_PVALUE",
                      "N_INSTRUMENTS", "METHOD", "MR_STATUS"]
        )

    genes = list(zip(
        candidate_genes["GENE_SYMBOL"],
        candidate_genes["ENSEMBL_ID"],
        candidate_genes["LEAD_SNP"],
    ))
    executor = RateLimitedExecutor(max_workers=max_workers, min_delay=0.1)
    results = executor.map(_run_mr_for_gene, genes)

    return pd.DataFrame(results)[
        ["GENE_SYMBOL", "MR_ESTIMATE", "MR_SE", "MR_PVALUE", "N_INSTRUMENTS", "METHOD", "MR_STATUS"]
    ].reset_index(drop=True)
```

**Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_mr_analysis.py::test_eqtl_instruments_from_tophits -v`
Expected: PASS

**Step 5: Commit**

```bash
git add src/gwas_explorer/mr_analysis.py tests/test_mr_analysis.py
git commit -m "feat(mr): use eQTLGen tophits for two-sample MR instruments"
```

---

### Task 3: Update existing tests for the new API flow

**Files:**
- Modify: `tests/test_mr_analysis.py`

The existing tests use `_make_opengwas_callback` which routes calls to `/associations` based on call count. The new flow hits `/tophits` for instruments and `/associations` for outcomes, so we need separate mocks for each endpoint.

**Step 1: Replace the callback helper and update all existing tests**

Delete `_make_opengwas_callback` entirely. Replace all existing `@responses.activate` tests to use:
- `responses.add(responses.POST, f"{OPENGWAS_API_URL}/tophits", json=instruments, status=200)` for instrument lookup
- `responses.add(responses.POST, f"{OPENGWAS_API_URL}/associations", json=outcomes, status=200)` for outcome lookup

Key changes per test:

- `test_wald_ratio_calculation`: mock `/tophits` returning 1 instrument, `/associations` returning 1 outcome
- `test_ivw_calculation`: mock `/tophits` returning 2 instruments, `/associations` returning 2 outcomes
- `test_insufficient_instruments`: mock `/tophits` returning `[]`
- `test_handles_api_error`: mock `/tophits` returning 500
- `test_wald_ratio_zero_beta`: mock `/tophits` returning 1 instrument with beta=0, `/associations` returning 1 outcome
- `test_allele_harmonization_flipped`: mock `/tophits` returning 1 instrument, `/associations` returning 1 outcome with swapped alleles
- `test_empty_candidates_returns_empty`: unchanged (no API calls)

**Step 2: Run all tests**

Run: `uv run pytest tests/test_mr_analysis.py -v`
Expected: ALL PASS

**Step 3: Commit**

```bash
git add tests/test_mr_analysis.py
git commit -m "test(mr): update all MR tests for eQTL tophits flow"
```

---

### Task 4: Add test for `no_eqtl_instruments` status

**Files:**
- Modify: `tests/test_mr_analysis.py`

**Step 1: Write the test**

```python
@responses.activate
def test_no_eqtl_data_for_gene(sample_candidate_genes: pd.DataFrame) -> None:
    """Genes without eQTLGen data should get no_eqtl_instruments status."""
    responses.add(
        responses.POST,
        f"{OPENGWAS_API_URL}/tophits",
        json=[],
        status=200,
    )

    genes = sample_candidate_genes[sample_candidate_genes["GENE_SYMBOL"] == "HLA-A"].copy()
    result = run_mr_analysis(genes)

    assert len(result) == 1
    assert result.iloc[0]["MR_STATUS"] == "no_eqtl_instruments"
```

**Step 2: Run test**

Run: `uv run pytest tests/test_mr_analysis.py::test_no_eqtl_data_for_gene -v`
Expected: PASS

**Step 3: Commit**

```bash
git add tests/test_mr_analysis.py
git commit -m "test(mr): add test for genes without eQTL instruments"
```

---

### Task 5: Update notebook for the new flow

**Files:**
- Modify: `notebooks/04_mendelian_randomization.ipynb`

**Step 1: Update notebook cell to show eQTL-based status breakdown**

The notebook code should work as-is since `run_mr_analysis` still takes `candidate_genes` and returns the same columns. But update the markdown header and add a note about the new eQTL approach. Also update the status display to include the new `no_eqtl_instruments` status.

**Step 2: Commit**

```bash
git add notebooks/04_mendelian_randomization.ipynb
git commit -m "docs(notebook): update MR notebook for eQTL two-sample approach"
```

---

### Task 6: Run full test suite and lint

**Step 1: Run all tests**

Run: `uv run pytest -v`
Expected: ALL PASS

**Step 2: Run linter**

Run: `uv run ruff check src/ tests/`
Expected: No errors

**Step 3: Run formatter**

Run: `uv run ruff format src/ tests/`

**Step 4: Final commit if formatting changed**

```bash
git add -u
git commit -m "style: format after eQTL MR refactor"
```
