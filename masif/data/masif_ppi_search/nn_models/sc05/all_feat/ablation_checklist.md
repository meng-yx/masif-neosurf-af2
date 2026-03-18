# Negative Sampling Ablation Checklist

## Runs
- `baseline_1to1`: `custom_params_baseline_1to1`
- `mixed_3to1`: `custom_params_mixed_3to1`
- `mixed_5to1`: `custom_params_mixed_5to1`

## Procedure
1. Rebuild cache with selected params.
2. Train MaSIF-ppi-search with same selected params.
3. Compare metrics from `model_data/log.txt`.

## Report metrics
- Validation ROC AUC
- Test ROC AUC
- Mean validation positive score
- Mean validation negative score
- Mean test positive score
- Mean test negative score
- Training stability (loss trend, gradient norm trend)

## Acceptance guideline
- Start with `mixed_3to1`.
- Promote to `mixed_5to1` only if validation/test AUC improves and training stays stable.
