import hashlib
from pathlib import Path

import pandas as pd

from src.qc import validate_inputs
from src.report import build_report
from src.statistics import benjamini_hochberg, differential_by_response


ROOT = Path(__file__).resolve().parents[1]


def test_synthetic_pipeline_recovers_expected_signal():
    proteins = pd.read_csv(ROOT / "data/synthetic_protein_matrix.csv")
    metadata = pd.read_csv(ROOT / "data/synthetic_clinical_metadata.csv")
    assert validate_inputs(proteins, metadata) == []
    results = differential_by_response(proteins, metadata)
    assert results.iloc[0]["protein"] in {"IL6", "CRP"}
    assert results["q_value"].between(0, 1).all()
    assert {"n_group_1", "n_group_2"}.issubset(results.columns)


def test_validation_reports_duplicate_sample_identifiers():
    proteins = pd.DataFrame({"sample_id": ["A", "A"], "P1": [1.0, 2.0]})
    metadata = pd.DataFrame({"sample_id": ["A", "B"], "response": ["x", "y"]})
    assert "protein matrix contains duplicated sample_id values" in validate_inputs(proteins, metadata)


def test_benjamini_hochberg_preserves_order_and_monotonicity_when_ranked():
    adjusted = benjamini_hochberg(pd.Series([0.04, 0.001, 0.02]))
    assert adjusted.tolist() == [0.04, 0.003, 0.03]


def test_run_manifest_records_input_fingerprints_and_tool_versions(tmp_path):
    output_dir = tmp_path / "run"
    build_report(
        str(ROOT / "data/synthetic_protein_matrix.csv"),
        str(ROOT / "data/synthetic_clinical_metadata.csv"),
        str(output_dir),
    )
    manifest = (output_dir / "run_manifest.txt").read_text()
    protein_hash = hashlib.sha256((ROOT / "data/synthetic_protein_matrix.csv").read_bytes()).hexdigest()
    metadata_hash = hashlib.sha256((ROOT / "data/synthetic_clinical_metadata.csv").read_bytes()).hexdigest()
    assert f"proteins_sha256={protein_hash}" in manifest
    assert f"metadata_sha256={metadata_hash}" in manifest
    assert "python_version=" in manifest
    assert "pandas_version=" in manifest
    assert "numpy_version=" in manifest
    assert "scipy_version=" in manifest
