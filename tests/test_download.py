"""Tests for GWAS Catalog bulk download module."""

from pathlib import Path
from unittest.mock import patch

import responses

from gwas_explorer.config import GWAS_CATALOG_FTP_URL
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
    import os
    import time

    existing = tmp_path / "gwas_catalog_associations.tsv"
    existing.write_text("old data")
    # Set mtime to 60 days ago
    old_time = time.time() - (60 * 86400)
    os.utime(existing, (old_time, old_time))

    fresh_tsv = "col1\tcol2\nfresh\tdata\n"
    responses.add(responses.GET, GWAS_CATALOG_FTP_URL, body=fresh_tsv, status=200)

    with patch("gwas_explorer.download.RAW_DIR", tmp_path):
        path = download_gwas_catalog(max_age_days=30)

    assert "fresh" in path.read_text()
