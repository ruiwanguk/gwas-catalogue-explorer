"""Pipeline configuration: endpoints, thresholds, and paths."""

import os
from pathlib import Path

from dotenv import load_dotenv

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Load .env from project root
load_dotenv(PROJECT_ROOT / ".env")
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"

# GWAS Catalog
GWAS_CATALOG_FTP_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"
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
EQTLGEN_STUDY_PREFIX = "eqtl-a-"  # eQTLGen blood cis-eQTL studies on OpenGWAS

# MR thresholds
MR_PVALUE_THRESHOLD = 0.05
MR_INSTRUMENT_PVALUE = 5e-8

# API tokens
OPENGWAS_TOKEN = os.environ.get("OPENGWAS_TOKEN")

# Download settings
CACHE_STALENESS_DAYS = 30
