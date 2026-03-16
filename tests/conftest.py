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
