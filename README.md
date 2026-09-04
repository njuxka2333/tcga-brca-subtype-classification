# TCGA-BRCA Molecular Subtype Classification

Classifies TCGA-BRCA breast tumors into 3-class IHC-surrogate molecular subtypes (Basal-like / HER2-positive / Luminal) from RNA-seq expression, using a leakage-safe nested cross-validation pipeline, followed by GO/KEGG biological interpretation and a shortcut-gene ablation.

## Files

- `TCGA_BRCA_LuminalBasal_pipeline.ipynb` — main analysis notebook: QC and preprocessing, nested cross-validated model comparison (LASSO, Random Forest, LightGBM, SVM), GO/KEGG enrichment, biological interpretation, and an AI usage appendix.
- `data/build_transcriptome_background.py` — builds the transcriptome-wide background gene list (from all TCGA-BRCA Primary Tumor RNA-seq samples on GDC) used as the statistical background for GO/KEGG enrichment in the notebook.

## Requirements

Python 3 with `pandas`, `numpy`, `scikit-learn`, `lightgbm`, `gseapy`, `matplotlib`, `requests`.
