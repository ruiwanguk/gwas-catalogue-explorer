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
