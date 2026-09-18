# ProteoTrace

**An auditable, research-use-only workflow for proteomics quality control, biomarker analysis, and evidence-based reporting.**

ProteoTrace is a portfolio project designed for translational proteomics research. It demonstrates a safe architecture in which deterministic Python functions perform all calculations, while a future language-model layer may only explain validated outputs and cite sources. It is not a diagnostic or clinical decision-support system.

## What it does

1. Validates protein-abundance and clinical-metadata files.
2. Measures missingness and flags outlier samples.
3. Tests whether batch is associated with outcome.
4. Compares responder and non-responder protein abundance with transparent Welch t-tests and Benjamini-Hochberg false-discovery-rate adjustment.
5. Exports a reproducible CSV report with quality-control findings and ranked candidate proteins.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m src.report --proteins data/synthetic_protein_matrix.csv --metadata data/synthetic_clinical_metadata.csv --output outputs
pytest
```

To open the local research-demo interface:

```bash
streamlit run app.py
```

## Data safety

The included files are entirely synthetic. Do not upload identifiable patient data, unpublished clinical data, or data that lacks ethical and institutional approval. See [data governance](governance/data_privacy.md).

See [synthetic demonstration results](docs/example_results.md) for the expected output and limitations.

## Repository structure

```text
data/        Synthetic demonstration data and dictionary
src/         Reproducible analysis functions
tests/       Automated checks
governance/  Model and data-governance documentation
outputs/     Generated reports, excluded from version control
```

## Portfolio roadmap

The first release delivers the deterministic analysis core and a local Streamlit interface. Planned work: pathway enrichment, literature retrieval with source traceability, a constrained Gemini-powered report writer, and evaluation on pre-registered synthetic benchmark scenarios.
Auditable AI workflow for proteomics quality control, biomarker analysis, and evidence-based reporting.
