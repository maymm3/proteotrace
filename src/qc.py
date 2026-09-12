"""Deterministic quality-control utilities for ProteoTrace."""

from __future__ import annotations

import pandas as pd


def validate_inputs(proteins: pd.DataFrame, metadata: pd.DataFrame) -> list[str]:
    """Return clear validation findings without altering input data."""
    findings: list[str] = []
    for name, frame in (("protein matrix", proteins), ("metadata", metadata)):
        if "sample_id" not in frame.columns:
            findings.append(f"{name} is missing required column: sample_id")
        elif frame["sample_id"].duplicated().any():
            findings.append(f"{name} contains duplicated sample_id values")
    if "sample_id" in proteins and "sample_id" in metadata:
        unmatched = set(proteins.sample_id).symmetric_difference(metadata.sample_id)
        if unmatched:
            findings.append(f"sample identifiers differ between files: {len(unmatched)} unmatched")
    return findings


def missingness_by_sample(proteins: pd.DataFrame) -> pd.DataFrame:
    values = proteins.drop(columns="sample_id")
    return pd.DataFrame({
        "sample_id": proteins.sample_id,
        "missing_fraction": values.isna().mean(axis=1),
    })


def outlier_samples(proteins: pd.DataFrame, z_threshold: float = 3.0) -> pd.DataFrame:
    """Flag samples with extreme mean abundance, using a documented robust rule."""
    values = proteins.drop(columns="sample_id")
    sample_means = values.mean(axis=1)
    center = sample_means.median()
    mad = (sample_means - center).abs().median()
    robust_z = pd.Series(0.0, index=sample_means.index) if mad == 0 else 0.6745 * (sample_means - center) / mad
    return pd.DataFrame({
        "sample_id": proteins.sample_id,
        "mean_abundance": sample_means,
        "robust_z": robust_z,
        "is_outlier": robust_z.abs() > z_threshold,
    })


def batch_response_table(metadata: pd.DataFrame) -> pd.DataFrame:
    if not {"batch", "response"}.issubset(metadata.columns):
        raise ValueError("metadata must contain batch and response columns")
    return pd.crosstab(metadata["batch"], metadata["response"])
