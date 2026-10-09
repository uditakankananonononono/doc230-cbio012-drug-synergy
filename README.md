# doc230-cbio012-drug-synergy
Audited: negative confirmed (bootstrap gain intervals below +0.01 for HGB and LR); gene features alone AUROC 0.665 but redundant with drug identity (RESULTS.md Audit).
Preregistered benchmark: do drug-target gene relationship features help predict drug-combination synergy (O Neil/Merck, Hetionet targets)? Answer here: no (cold-pair gate missed). See PREREG.md, RESULTS.md, src/run.py, data_manifest.

## Reproduce / readability audit (2026-10-09)
Data: sources and SHA256 are in data_manifest/. Scripts read fixed paths under /tmp (download the files there first, check the SHA256, then run python3 src/run.py and python3 src/audit.py). Not independently re-run from a clean machine as part of this readability check; the numbers in RESULTS.md come from the original runs on 2 CPU / 2 GB. The O'Neil labels file (labels.csv, DeepSynergy) goes in /tmp/syn/ and Hetionet nodes.tsv and edges.sif.gz in /tmp/het/.
