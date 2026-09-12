"""Transparent group-comparison functions for ProteoTrace."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def differential_by_response(proteins: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Run Welch tests for every protein and return a ranked, reproducible table."""
    merged = metadata[["sample_id", "response"]].merge(proteins, on="sample_id", validate="one_to_one")
    groups = merged["response"].dropna().unique().tolist()
    if len(groups) != 2:
        raise ValueError("response must contain exactly two groups for this demonstration")
    first, second = groups
    rows = []
    for protein in proteins.columns.drop("sample_id"):
        x = merged.loc[merged.response == first, protein].dropna()
        y = merged.loc[merged.response == second, protein].dropna()
        statistic, p_value = ttest_ind(x, y, equal_var=False)
        rows.append({
            "protein": protein,
            "group_1": first,
            "group_2": second,
            "mean_group_1": x.mean(),
            "mean_group_2": y.mean(),
            "mean_difference": x.mean() - y.mean(),
            "p_value": p_value,
        })
    result = pd.DataFrame(rows)
    result["abs_mean_difference"] = result.mean_difference.abs()
    return result.sort_values(["p_value", "abs_mean_difference"], ascending=[True, False]).reset_index(drop=True)
