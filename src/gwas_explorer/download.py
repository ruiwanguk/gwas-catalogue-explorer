"""Download GWAS Catalog bulk association data."""

import logging
import tempfile
import time
import zipfile
from pathlib import Path

from gwas_explorer.config import (
    CACHE_STALENESS_DAYS,
    GWAS_ASSOCIATIONS_FILENAME,
    GWAS_CATALOG_FTP_URL,
    RAW_DIR,
)
from gwas_explorer.http_utils import create_session

logger = logging.getLogger(__name__)

_session = create_session()


def download_gwas_catalog(max_age_days: int = CACHE_STALENESS_DAYS) -> Path:
    """Download GWAS Catalog associations ZIP and extract the TSV.

    Streams to a temporary file to avoid loading the entire ZIP into memory.

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
            logger.info("Using cached file (%d days old)", int(age_days))
            return filepath

    logger.info("Downloading GWAS Catalog from %s", GWAS_CATALOG_FTP_URL)
    response = _session.get(GWAS_CATALOG_FTP_URL, stream=True, timeout=300)
    response.raise_for_status()

    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        for chunk in response.iter_content(chunk_size=65536):
            tmp.write(chunk)

    try:
        with zipfile.ZipFile(tmp_path) as zf:
            tsv_names = [n for n in zf.namelist() if n.endswith(".tsv")]
            if not tsv_names:
                raise RuntimeError(f"No TSV file found in archive: {zf.namelist()}")
            with zf.open(tsv_names[0]) as src, open(filepath, "wb") as dst:
                while True:
                    chunk = src.read(65536)
                    if not chunk:
                        break
                    dst.write(chunk)
    finally:
        tmp_path.unlink(missing_ok=True)

    logger.info("Extracted to %s", filepath)
    return filepath
