# Synthetic data dictionary

All values are synthetic log2-normalised protein-abundance values created solely for demonstration.

| Field | Meaning |
|---|---|
| `sample_id` | Synthetic, non-identifying sample identifier |
| `response` | Synthetic study-group label |
| `batch` | Synthetic analytical batch label |
| `timepoint` | Synthetic collection timepoint |
| Protein columns | Synthetic abundance values for named proteins |

The synthetic signal intentionally makes IL6 and CRP higher in the responder group. This enables transparent testing of the differential-analysis workflow.
