# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GWAS Catalogue exploration project using Python 3.14. Focuses on querying, analyzing, and visualizing data from the NHGRI-EBI GWAS Catalogue.

## Development Commands

This project uses `uv` as the package manager.

```bash
# Install dependencies
uv sync

# Run Jupyter notebooks
uv run jupyter notebook

# Add a new dependency
uv add <package-name>
```

## Project Structure

- `notebooks/` - Jupyter notebooks for data exploration and analysis
- `data/` - Downloaded datasets and cached API responses
- `main.py` - Application entry point
- `pyproject.toml` - Project configuration and dependencies

## Key Resources

- GWAS Catalogue REST API: https://www.ebi.ac.uk/gwas/rest/api
- GWAS Catalogue website: https://www.ebi.ac.uk/gwas/
