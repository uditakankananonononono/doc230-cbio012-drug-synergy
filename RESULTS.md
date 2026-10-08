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
