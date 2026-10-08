# CBIO012 results (results/results.json, src/run.py, results/run.log). Replication-style benchmark; hypothesis NOT supported.
Data: O'Neil/Merck screen (DeepSynergy labels), Hetionet CbG targets. Only 24 of 38 drugs map to targets, giving 6,084 of 23,052 rows (156 drug pairs, 39 cell lines, synergy rate 26%). 14 unmapped drugs (mostly investigational Merck compounds) are excluded. SN-38 uses irinotecan as a prodrug proxy. No cell-line expression or chemistry features.
AUROC, HistGradientBoosting, identity-only vs identity+gene-relationship features (mean over 3 seeds x 5 folds):
| Split | identity only | + gene features | gain |
|---|---|---|---|
| random rows (D2) | 0.823 | 0.836 | +0.014 |
| cold drug-pair (D3, gate >= +0.01) | 0.762 | 0.726 | -0.037 NOT met |
| cold cell line (D4) | 0.791 | 0.804 | +0.013 |
| cold drug (D5, 24 drugs; gain in 4 of 24) | 0.749 | 0.706 | -0.043 |
- D2 baseline: logistic regression on identity only, cold-pair 0.746.
- D6 ablation (cold-pair gain): overlap -0.035, network links -0.024, pathway -0.021, SVD -0.010. No feature group helps.
- D7 shuffled gene features (cold-pair): gain -0.069, below the real-feature -0.037. Gate (<= +0.01) met, but real features are only less harmful than shuffled ones, not helpful.
- D8 per-cell-line AUROC (cold-pair) 0.52-0.84, mean 0.69.
- D9 threshold 5/10/20: gain -0.028, -0.032, -0.028 (consistent).
- D10 regression Pearson r (cold-pair): identity 0.404, with gene 0.399.
Reading: target-gene relationship features did not help predict synergy for unseen drug pairs or unseen drugs. Small gains on random rows and cold cell lines (+0.013) are the kind that appear when pair identity can be partly recovered from the features, so they are not evidence for the gene-synergy idea. Limits: 24 drugs, hand-limited target coverage (Hetionet drug targets are sparse and old), no expression features. 0 novel findings.

## Audit (second pass) - the negative is solid for added value; gene features alone carry signal that identity already captures
Prereg AUDIT_PREREG.md committed before running (one syntax fix to audit.py before any result existed); src/audit.py, results/audit.json. Labels SHA256 re-verified. Cold-pair, seeds 0-2 averaged out-of-fold predictions, pair-level cluster bootstrap (200 resamples of 156 pairs).
- Q1 HGB gain (identity+gene minus identity): -0.030, 95% interval [-0.057, -0.007] (AUROC 0.769 -> 0.739). Entirely below +0.01.
- Q2 logistic regression with the 6 non-SVD gene features (overlap, links, pathway): gain -0.006, interval [-0.012, -0.002] (0.758 -> 0.753). Entirely below +0.01. Removing the SVD block makes the harm tiny but does not make it positive.
- Q3 gene features alone (no identity), HGB: AUROC 0.665, interval [0.612, 0.718], above 0.5.
Pre-specified reading: negative SOLID (both intervals below +0.01). Q3 says the gene-relationship features are not empty: alone they predict synergy above chance. But they add nothing on top of drug identity, and several of them (summed/product target counts) are drug-level quantities, so Q3 cannot be read as gene-relationship signal; it is consistent with them acting as proxies for which drugs tend to be synergistic. Hypothesis still NOT supported in this benchmark (24 drugs, sparse old Hetionet targets, no expression features); this does not rule out the idea with better data.
