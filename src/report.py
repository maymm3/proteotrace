"""Command-line runner that exports traceable ProteoTrace outputs."""

from __future__ import annotations

import argparse
import hashlib
import platform
from pathlib import Path

import numpy
import pandas as pd
import scipy

from src.qc import batch_response_table, missingness_by_sample, outlier_samples, validate_inputs
from src.statistics import differential_by_response


def _sha256(path: str) -> str:
    """Return a streaming SHA-256 digest for an input file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_report(proteins_path: str, metadata_path: str, output_dir: str) -> None:
    proteins = pd.read_csv(proteins_path)
    metadata = pd.read_csv(metadata_path)
    findings = validate_inputs(proteins, metadata)
    if findings:
        raise ValueError("Input validation failed: " + "; ".join(findings))
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    missingness_by_sample(proteins).to_csv(out / "sample_missingness.csv", index=False)
    outlier_samples(proteins).to_csv(out / "sample_outliers.csv", index=False)
    batch_response_table(metadata).to_csv(out / "batch_by_response.csv")
    differential_by_response(proteins, metadata).to_csv(out / "differential_proteins.csv", index=False)
    (out / "run_manifest.txt").write_text(
        "ProteoTrace deterministic analysis run\n"
        f"python_version={platform.python_version()}\n"
        f"pandas_version={pd.__version__}\n"
        f"numpy_version={numpy.__version__}\n"
        f"scipy_version={scipy.__version__}\n"
        f"proteins={proteins_path}\nmetadata={metadata_path}\n"
        f"proteins_sha256={_sha256(proteins_path)}\n"
        f"metadata_sha256={_sha256(metadata_path)}\n"
        "Research use only. Results require expert scientific review.\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the ProteoTrace synthetic demonstration.")
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    build_report(args.proteins, args.metadata, args.output)
